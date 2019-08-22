


function DQ(mat::AMat, d::AVec, q::AFloat; J::Int64=1200, sd::Float64=0.335)

    # Generate uniform 3D unit vectors
    n = length(d);
    x = [(2*j - 1 - J)/J for j = 1:J];
    theta, phi = acos.(x), sqrt(pi*J)*asin.(x);
    qmat = q * [sin.(theta) .* cos.(phi) sin.(theta) .* sin.(phi) cos.(theta)];

    # Calculate I(q)
    qr = mat * qmat';
    A = cos.(qr);
    B = sin.(qr);
    sv = ones(n);

    # the density d in the input is the "excessive density" above the solvents
    d = d + sd * sv;
    solute = [d' * A; -d' * B];
    solvent = sd * [sv' * A; -sv' * B];
    subtracted = solute .- solvent;
    D = mean(sum((subtracted).^2, dims=1));

    return D;

end




function DenSWAXS(mat::AMat, d::AVec, q::AVec; J::Int64=1200, sden::Float64=0.335)

    # check the dimension matching
    n = size(mat, 1);
    length(d) == n ? nothing : error("Dimensions mismatch betewwn grid points and densities ...");

    # Calculate electron density of the solvent
    spacing = norm(mat[1, :] .- mat[2, :]);
    solventden = sden * spacing ^ 3;

    # Calculate the intensity using DQ function
    intensity = pmap(x -> DQ(mat, d, x; sd=solventden, J=J), q, distributed=true);

    return intensity;
end


# Overload DenSWAXS for convience of back calculation
function DenSWAXS(eden::mrc, support::mrc, q::AVec; support_cutoff::Float64=0.75, J::Int64=1200, sded::Float64=0.335)

    eden.nx != support.nx ? error("DenSWAXS: The ::mrc shapes do not match!") : nothing;

    # Calculate electron density of the solvent
    spacing = norm(eden.gridposition[1, :] .- eden.gridposition[2, :]);
    solventden = sden * spacing ^ 3;

    # get the density within the molecular support
    idx = findall(reshape(support.data, (support.nx*support.ny*support.nz, )) .> support_cutoff);
    mat = eden.gridposition[idx, :];
    d = reshape(eden.data, (eden.nx*eden.ny*eden.nz, ))[idx];

    # Calculate the intensity using DQ function
    intensity = pmap(x -> DQ(mat, d, x; sd=solventden, J=J), q, distributed=true);

    return intensity;
end
