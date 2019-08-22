
function ShapeSWAXS(v::Voxel, d::AFloat, q::AVec; J::Int64=1800)

    mat = v.coordata;
    intensity = pmap(x -> _DQ(mat, d * ones(size(mat, 1)), x; sd=0.00, J=J), q, distributed=true);

    return intensity;
end
