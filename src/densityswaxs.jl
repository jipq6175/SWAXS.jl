mutable struct mrc
    filename::String;
    nx::Int32;
    ny::Int32;
    nz::Int32;
    mode::Int32;
    nxstart::Int32;
    nystart::Int32;
    nzstart::Int32;
    mx::Int32;
    my::Int32;
    mz::Int32;
    cella::Array{Float32, 1};
    cellb::Array{Float32, 1};
    mapc::Int32;
    mapr::Int32;
    maps::Int32;
    dmin::Float32;
    dmax::Float32;
    dmean::Float32;
    ispg::Int32;
    nsymbt::Int32;
    extra::Array{UInt8, 1};
    origin::Array{Float32, 1};
    map::Array{UInt8, 1};
    machst::Array{UInt8, 1};
    rms::Float32;
    nlabl::Int32;
    label::Array{UInt8, 1};
    data::Array{Float32, 3};
    gridposition::Array{Float32, 2};
end


function mrc_init()
    zi = Int32(0);
    zf = Float32(0.0);
    zu = UInt8(0);
    m = mrc("filename",
            ntuple(x->zi, 10)...,
            [zf],
            [zf],
            ntuple(x->zi, 3)...,
            ntuple(x->zf, 3)...,
            ntuple(x->zi, 2)...,
            [zu],
            [zf],
            [zu],
            [zu],
            zf,
            zi,
            [zu],
            Float32.(zeros(1,1,1)),
            Float32.(zeros(1,3)));
    return m;
end


function mrc_reader(filename::String)

    f = open(filename);
    words = read(f);
    close(f);

    m = mrc(filename,
            ntuple(x -> reinterpret(Int32, words[1:40])[x], 10)...,
            Array{Float32, 1}(reinterpret(Float32, words[41: 52])),
            Array{Float32, 1}(reinterpret(Float32, words[53: 64])),
            ntuple(x -> reinterpret(Int32, words[65:76])[x], 3)...,
            ntuple(x -> reinterpret(Float32, words[77:88])[x], 3)...,
            ntuple(x -> reinterpret(Int32, words[89:96])[x], 2)...,
            words[97:196],
            Array{Float32, 1}(reinterpret(Float32, words[197: 208])),
            words[209:212],
            words[213:216],
            reinterpret(Float32, words[217:220])[1],
            reinterpret(Int32, words[221:224])[1],
            words[225:1024],
            Float32.(zeros(1,1,1)),
            Float32.(zeros(1,3)));

    length(words) - 1024 - 4 * prod([m.nx; m.ny; m.nz]) == 0 ? nothing : @error("Dimensions mismatch.");

    m.data = reshape(Array{Float32, 1}(reinterpret(Float32, words[1025:end])), (m.nx, m.ny, m.nz));

    # need to correct the grid positions
    m.gridposition = fill(Float32(0.0), (m.nx*m.ny*m.nz, 3));
    length(unique(m.cella)) != 1 ? @warn("The cell is not isotropic. ") : nothing;
    grid_spacing = sign(m.origin[1]) * unique(m.cella)[1] / m.nx;
    id = 1;
    for i = 1:m.nx, j = 1:m.ny, k = 1:m.nz
        m.gridposition[id, :] = grid_spacing .* ([i j k] .- 1.0);
        id = id + 1;
    end

    return m;
end



function mrc_print(m::mrc)
    println("===============  MRC Header Summary  ===============");
    println("(nx, ny, nz) = ($(m.nx), $(m.ny), $(m.nz)).");
    println("(mx, my, mz) = ($(m.mx), $(m.my), $(m.mz)).");
    println("(nxstart, nystart, nzstart) = ($(m.nxstart), $(m.nystart), $(m.nzstart)).");
    println("Cell dimensions = ($(m.cella[1]), $(m.cella[2]), $(m.cella[3])) (Å).");
    println("Cell angles = ($(m.cellb[1]), $(m.cellb[2]), $(m.cellb[3])) (degrees).");
    println("Axis correspondance = ($(m.mapc), $(m.mapr), $(m.maps)).");
    println("(dmin, dmax, dmean) = ($(m.dmin), $(m.dmax), $(m.dmean)) (e), RMS = $(m.rms).");
    println("Origin = ($(m.origin[1]), $(m.origin[2]), $(m.origin[3])).");
    println("====================================================");
    return nothing;
end




function mrc_writer(filename::String, m::mrc)

    n = m.nx * m.ny * m.nz;
    words = Array{UInt8, 1}(undef, 0);

    # push the headers
    append!(words, reinterpret(UInt8, [m.nx; m.ny; m.nz; m.mode; m.nxstart; m.nystart; m.nzstart; m.mx; m.my; m.mz]) |> Vector);
    append!(words, reinterpret(UInt8, [m.cella; m.cellb]) |> Vector);
    append!(words, reinterpret(UInt8, [m.mapc; m.mapr; m.maps]) |> Vector);
    append!(words, reinterpret(UInt8, [m.dmax; m.dmin; m.dmean]) |> Vector);
    append!(words, reinterpret(UInt8, [m.ispg; m.nsymbt]) |> Vector);
    append!(words, m.extra);
    append!(words, reinterpret(UInt8, m.origin) |> Vector);
    append!(words, [m.map; m.machst]);
    append!(words, reinterpret(UInt8, [m.rms]) |> Vector);
    append!(words, reinterpret(UInt8, [m.nlabl]) |> Vector);
    append!(words, m.label);

    # push the data
    append!(words, reinterpret(UInt8, reshape(m.data, (n,))) |> Vector);

    # write the data into file
    length(words) == 4n + 1024 ? nothing : @warn("Length mismatch!");
    f = open(filename, "w");
    write(f, words);
    close(f);
    @info("MRC: $filename written successfully!!");

    return nothing;
end




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
