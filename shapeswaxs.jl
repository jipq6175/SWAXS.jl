

function DQ(mat::AMat, q::AFloat; J::Int64=1500)

    # Generate uniform 3d unit vectors
    x = [(2*j - 1 - J)/J for j = 1:J];
    theta, phi = acos.(x), sqrt(pi*J)*asin.(x);
    qmat = q .*[sin.(theta) .* cos.(phi) sin.(theta) .* sin.(phi) cos.(theta)];

    # Calculate I(q)
    qr = mat * qmat';
    sys = [sum(cos.(qr), 1); sum(sin.(qr), 1)];
    D = mean(sum(sys.^2, 1));
    return D;
end



function ShapeWAXS(mat::AMat, q::AVec; J::Int64=1800)



    # Include self everywhere
    @eval @everywhere JJ = $J;

    ####################################
    # Creat a function for pmap later
    @everywhere function DQ(mat::Matrix{T}, q::T; J::Int64=JJ) where T<:Real

        # Generate uniform 3d unit vectors
        x = [(2*j - 1 - J)/J for j = 1:J];
        theta, phi = acos.(x), sqrt(pi*J)*asin.(x);
        qmat = q .*[sin.(theta) .* cos.(phi) sin.(theta) .* sin.(phi) cos.(theta)];

        # Calculate I(q)
        qr = mat * qmat';
        sys = [sum(cos.(qr), 1); sum(sin.(qr), 1)];
        D = mean(sum(sys.^2, 1));
        return D;
    end
    ###################################

    # Compute the scattering profile using pmap
    @time intensity = pmap(x -> DQ(mat, x), q);

    return intensity;

end
