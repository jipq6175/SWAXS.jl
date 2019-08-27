

function _Orie(q::AFloat, J::Signed)
    x = [(2*j - 1 - J)/J for j = 1:J];
    theta, phi = acos.(x), sqrt(pi*J)*asin.(x);
    qmat = q * [sin.(theta) .* cos.(phi) sin.(theta) .* sin.(phi) cos.(theta)];
    return qmat;
end



function _DQ(mat::AMat, d::AVec, q::AFloat; J::Int64=1200, sd::Float64=0.335)

    n = length(d);
    qmat = _Orie(q, J);

    qr = mat * qmat';
    A = cos.(qr);
    B = sin.(qr);
    sv = ones(n);

    d = d + sd * sv;
    solute = [d' * A; -d' * B];
    solvent = sd * [sv' * A; -sv' * B];
    subtracted = solute .- solvent;
    D = mean(sum((subtracted).^2, dims=1));
    return D;
end




function DenSWAXS(m::mrc, q::AVec; density_cutoff::Float64=1e-3, J::Int64=1200, sden::Float64=0.335)

    spacing = norm(m.gridposition[1, :] .- m.gridposition[2, :]);
    solventden = sden * spacing ^ 3;

    d = reshape(m.data, (m.nx*m.ny*m.nz, ));
    idx = findall((d .>= density_cutoff) .| (d .<= -density_cutoff));
    mat = m.gridposition[idx, :];

    intensity = pmap(x -> _DQ(mat, d[idx], x; sd=solventden, J=J), q, distributed=true);
    return intensity;
end
