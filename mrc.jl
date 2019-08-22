# The mrc/ccp4 module for DENSS volumetric electron densities

__precompile__()


# The mrc module dealing with the 3D electron densities
module MRC

# Required packages within the module
using PyCall, Statistics

# setting alias of type
const AFloat = AbstractFloat;

# Export all the functions defined
export mrc
# export mrc_init
# export mrc_copy
export mrc_reader, mrc_print, mrc_writer
export mrc_zoom, mrc_zoom!, mrc_gaussian, mrc_gaussian!, mrc_scale, mrc_scale!
export mrc_chop, mrc_chop!, mrc_assign, mrc_assign!
export mrc_connect!, mrc_connect, mrc_greyconnect!, mrc_greyconnect
export mrc_convolve!, mrc_convolve


##############################################################################################################
 #     # ######   #####      #####  ####### ######  #     #  #####  #######
 ##   ## #     # #     #    #     #    #    #     # #     # #     #    #
 # # # # #     # #          #          #    #     # #     # #          #
 #  #  # ######  #           #####     #    ######  #     # #          #
 #     # #   #   #                #    #    #   #   #     # #          #
 #     # #    #  #     #    #     #    #    #    #  #     # #     #    #
 #     # #     #  #####      #####     #    #     #  #####   #####     #
##############################################################################################################

"""
    m::mrc

The mutable `mrc` struct for storing the 3D electron densities information from .mrc (ccp4) file whose data structure is all in 32 bits, i.e. `Float32` or `Int32`.

Summary\\
≡≡≡≡≡≡≡≡≡

    mutable struct mrc <: Any

Fields\\
≡≡≡≡≡≡≡≡

    filename     :: String
    nx           :: Int32
    ny           :: Int32
    nz           :: Int32
    mode         :: Int32
    nxstart      :: Int32
    nystart      :: Int32
    nzstart      :: Int32
    mx           :: Int32
    my           :: Int32
    mz           :: Int32
    cella        :: Array{Float32,1}
    cellb        :: Array{Float32,1}
    mapc         :: Int32
    mapr         :: Int32
    maps         :: Int32
    dmin         :: Float32
    dmax         :: Float32
    dmean        :: Float32
    ispg         :: Int32
    nsymbt       :: Int32
    extra        :: Array{UInt8,1}
    origin       :: Array{Float32,1}
    map          :: Array{UInt8,1}
    machst       :: Array{UInt8,1}
    rms          :: Float32
    nlabl        :: Int32
    label        :: Array{UInt8,1}
    data         :: Array{Float32,3}
    gridposition :: Array{Float32,2}
"""
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




##############################################################################################################
  #####  ####### #     #  #####  ####### ######  #     #  #####  ####### ####### ######
 #     # #     # ##    # #     #    #    #     # #     # #     #    #    #     # #     #
 #       #     # # #   # #          #    #     # #     # #          #    #     # #     #
 #       #     # #  #  #  #####     #    ######  #     # #          #    #     # ######
 #       #     # #   # #       #    #    #   #   #     # #          #    #     # #   #
 #     # #     # #    ## #     #    #    #    #  #     # #     #    #    #     # #    #
  #####  ####### #     #  #####     #    #     #  #####   #####     #    ####### #     #
##############################################################################################################

# The NULL mrc constructor
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


"""
    m = mrc_reader(filename::String)

Reads in the `m::mrc` struct from a .mrc file.
"""
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




##############################################################################################################
#     # ####### ### #        #####
#     #    #     #  #       #     #
#     #    #     #  #       #
#     #    #     #  #        #####
#     #    #     #  #             #
#     #    #     #  #       #     #
 #####     #    ### #######  #####
 #############################################################################################################


"""
    mrc_print(m)

Prints out the summary of `mrc` struct.
"""
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



"""
    mrc_writer(filename::String, m::mrc)

Writes the `m::mrc` struct into .mrc file to be visualized by PyMOL or other programs.
"""
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



# Perform a copy of the mrc struct
function mrc_copy(m::mrc)
    m2 = mrc_init();
    fields = fieldnames(typeof(m2));
    for i = 1:length(fields)
        setfield!(m2, fields[i], deepcopy(getfield(m, fields[i])));
        # This has to use deepcopy because the getfield refer to the memory of the field.
        # we need deepcopy to make it independent from the original one
    end

    return m2;
end




##############################################################################################################
 ####### ####### ####### #     #
      #  #     # #     # ##   ##
     #   #     # #     # # # # #
    #    #     # #     # #  #  #
   #     #     # #     # #     #
  #      #     # #     # #     #
 ####### ####### ####### #     #
##############################################################################################################

# zoom in or zoom out, corresponding to the SWAXS resolution and oversampling
# This function changes the original mrc struct and that's why there is a "!"
"""
    mrc_zoom!(m::mrc, zoom::Float64; order::Int64=4, mode::String="constant", verbose::Bool=true)

zooms in (`zoom > 1.0`) or zooms out (`zoom < 1.0`) from the current `m::mrc`. The `!` implies that the `m::mrc` will be changed. \\
zoom > 1.0: zooming in, corresponding to voxel splitting into finer 3D resolution.\\
zoom < 1.0: zooming out, corresponding to voxel merging into a coarser 3D resolution.\\
zoom = 1.0: the `m::mrc` is unchanged, resolution remains.\\
For practical considerations, `0.8 < zoom < 1.75` is a good range for slightly SWAXS differences. \\
It's suggested that `zoom = 1.25` or `1.5` for electron density split. \\
"""
function mrc_zoom!(m::mrc, zoom::AFloat; order::Int64=4, mode::String="constant", verbose::Bool=true, normalize::Bool=true)

    ndimage = pyimport("scipy.ndimage");
    zoomed = ndimage.zoom(m.data, zoom, order=order, mode=mode);
    (nx, ny, nz) = Int32.(size(zoomed));
    m.nxstart = sign(m.nxstart) * Int32(floor((nx-1)/2));
    m.nystart = sign(m.nystart) * Int32(floor((ny-1)/2));
    m.nzstart = sign(m.nzstart) * Int32(floor((nz-1)/2));
    m.nx, m.ny, m.nz, m.mx, m.my, m.mz = nx, ny, nz, nx, ny, nz;
    xunit, yunit, zunit = m.cella[1]/nx, m.cella[2]/ny, m.cella[3]/nz;
    m.origin = Float32.([m.nxstart*xunit; m.nystart*yunit; m.nzstart*zunit]);
    normalize ? scale = sum(m.data) / sum(zoomed) : scale = 1.0;
    m.dmin = Float32(scale * minimum(zoomed));
    m.dmax = Float32(scale * maximum(zoomed));
    m.dmean = Float32(scale * mean(zoomed));
    m.rms = Float32(scale * std(zoomed));
    m.data = Float32.(scale * zoomed);

    # deal with grid positions
    m.gridposition = fill(Float32(0.0), (m.nx*m.ny*m.nz, 3));
    length(unique(m.cella)) != 1 ? @warn("The cell is not isotropic. ") : nothing;
    grid_spacing = sign(m.origin[1]) * unique(m.cella)[1] / m.nx;
    id = 1;
    for i = 1:m.nx, j = 1:m.ny, k = 1:m.nz
        m.gridposition[id, :] = grid_spacing .* ([i j k] .- 1.0);
        id = id + 1;
    end

    verbose ? @info("MRC: The input ::mrc was modified (zoomed) successfully ...") : nothing;

end


# zoom in or zoom out, corresponding to the SWAXS resolution and oversampling
# This function outputs a zoomed mrc struct
"""
    m2 = mrc_zoom(m::mrc, zoom::Float64; order::Int64=4, mode::String="constant")

Zooms in or zooms out of `mrc` and returns the zoomed `mrc` leaving the input intact. For more detailed explanation, type `?mrc_zoom!`.
"""
function mrc_zoom(m::mrc, zoom::AFloat; order::Int64=4, mode::String="constant", normalize::Bool=true)

    m2 = mrc_copy(m);
    mrc_zoom!(m2, zoom; order=order, mode=mode, verbose=false, normalize=normalize);

    @info("MRC: Operation completed (zoomed) successfully ...");
    return m2;

end




##############################################################################################################
  #####     #    #     #  #####   #####  ###    #    #     #
 #     #   # #   #     # #     # #     #  #    # #   ##    #
 #        #   #  #     # #       #        #   #   #  # #   #
 #  #### #     # #     #  #####   #####   #  #     # #  #  #
 #     # ####### #     #       #       #  #  ####### #   # #
 #     # #     # #     # #     # #     #  #  #     # #    ##
  #####  #     #  #####   #####   #####  ### #     # #     #
##############################################################################################################

"""
    mrc_gaussian!(m::mrc, sigma::Float64; order::Int64=0, mode::String="constant", verbose::Bool=true)

Blurs the `m::mrc` using Gaussian Filter with parameter `sigma::Float64`. The size of `mrc` is unchanged but the values are smeared. The `!` implies that the input `m::mrc` is changed. \\
For practical considerations, `sigma` should be between `0.25` and `0.8` for 3D electron densities. \\
It is suggested that `0.25 < sigma < 0.5` for slightly SWAXS differences. \\
"""
function mrc_gaussian!(m::mrc, sigma::AFloat; order::Int64=0, mode::String="constant", verbose::Bool=true)

    ndimage = pyimport("scipy.ndimage");
    gaussianed = ndimage.gaussian_filter(m.data, sigma, order=order, mode=mode);
    scale = sum(m.data) / sum(gaussianed);
    m.dmin = Float32(scale * minimum(gaussianed));
    m.dmax = Float32(scale * maximum(gaussianed));
    m.dmean = Float32(scale * mean(gaussianed));
    m.rms = Float32(scale * std(gaussianed));
    m.data = Float32.(scale * gaussianed);

    verbose ? @info("MRC: The input ::mrc was modified (gaussianed) successfully ...") : nothing;

end


"""
    m2 = mrc_gaussian(m::mrc, sigma::Float64; order::Int64=0, mode::String="constant")

Blurs the `m::mrc` using Gaussian Filter with parameter `sigma::Float64` and returns the blurred `m2::mrc`. For detailed information, type `?mrc_gaussian!`.\\
"""
function mrc_gaussian(m::mrc, sigma::AFloat; order::Int64=0, mode::String="constant")

    m2 = mrc_copy(m);
    mrc_gaussian!(m2, sigma; order=order, mode=mode, verbose=false);

    @info("MRC: Operation completed (gaussianed) successfully ...");
    return m2;

end




##############################################################################################################
  #####   #####     #    #       #######
 #     # #     #   # #   #       #
 #       #        #   #  #       #
  #####  #       #     # #       #####
       # #       ####### #       #
 #     # #     # #     # #       #
  #####   #####  #     # ####### #######
##############################################################################################################

# Scale the mrc struct: scale the overall electron density without changing its shape
function mrc_scale!(m::mrc, scale::AFloat; verbose::Bool=true)

    m.dmin = Float32(scale * m.dmin);
    m.dmax = Float32(scale * m.dmax);
    m.dmean = Float32(scale * m.dmean);
    m.rms = Float32(scale * m.rms);
    m.data = Float32.(scale * m.data);
    verbose ? @info("MRC: The input ::mrc was modified (scaled) successfully ...") : nothing;

end


function mrc_scale(m::mrc, scale::AFloat)

    m2 = mrc_copy(m);
    mrc_scale!(m2, scale; verbose=false);

    @info("MRC: Operation completed (scaled) successfully ...");
    return m2;

end




##############################################################################################################
  #####  #     # ####### ######
 #     # #     # #     # #     #
 #       #     # #     # #     #
 #       ####### #     # ######
 #       #     # #     # #
 #     # #     # #     # #
  #####  #     # ####### #
##############################################################################################################

# chop out the voxels that have lover than threshold electron density
function mrc_chop!(m::mrc, threshold::AFloat; verbose::Bool=true)

    idx = findall(m.data .> threshold);
    cut = maximum(abs.(getindex.(idx, [1 2 3]) .+ [m.nxstart m.nystart m.nzstart])) + 2;

    if (cut < abs(m.nxstart)) && (cut < abs(m.nystart)) && (cut < abs(m.nzstart))
        # this is when everything in data is included.
        support = fill(false, size(m.data));
        support[abs(m.nxstart)-cut:abs(m.nxstart)+cut, abs(m.nystart)-cut:abs(m.nystart)+cut, abs(m.nzstart)-cut:abs(m.nzstart)+cut] = fill(true, (2*cut+1, 2*cut+1, 2*cut+1));

        m.data = m.data[abs(m.nxstart)-cut:abs(m.nxstart)+cut, abs(m.nystart)-cut:abs(m.nystart)+cut, abs(m.nzstart)-cut:abs(m.nzstart)+cut];

        spacing = m.origin[1] / m.nxstart;
        m.cella = Float32.(spacing * [2*cut+1; 2*cut+1; 2*cut+1]);
        m.nxstart = Int32(sign(m.nxstart) * (cut+1));
        m.nystart = Int32(sign(m.nystart) * (cut+1));
        m.nzstart = Int32(sign(m.nzstart) * (cut+1));
        m.gridposition = m.gridposition[reshape(support, (m.nx*m.ny*m.nz, )), :];
        m.nx, m.ny, m.nz, m.mx, m.my, m.mz = 2*cut+1, 2*cut+1, 2*cut+1, 2*cut+1, 2*cut+1, 2*cut+1;
        m.dmin = Float32(minimum(m.data));
        m.dmean = Float32(mean(m.data));
        m.rms = Float32(std(m.data));
        verbose ? @info("MRC: The input ::mrc was modified (chopped) successfully ...") : nothing;
    else
        verbose ? @info("MRC: The input ::mrc was not modified (chopped) at this threshold = $threshold ...") : nothing;
    end

end



function mrc_chop(m::mrc, threshold::AFloat)

    m2 = mrc_copy(m);
    mrc_chop!(m2, threshold; verbose=false);
    size(m2.data) == size(m.data) ? @info("MRC: The input ::mrc was not modified (chopped) at this threshold = $threshold ...") : @info("MRC: Operation completed (chopped) successfully ...");

    return m2;
end



# chop out the voxels using cut::Int64
function mrc_chop!(m::mrc, cut::Signed; verbose::Bool=true)
    if (cut < abs(m.nxstart)) && (cut < abs(m.nystart)) && (cut < abs(m.nzstart))
        # this is when everything in data is included.
        support = fill(false, size(m.data));
        support[abs(m.nxstart)-cut:abs(m.nxstart)+cut, abs(m.nystart)-cut:abs(m.nystart)+cut, abs(m.nzstart)-cut:abs(m.nzstart)+cut] = fill(true, (2*cut+1, 2*cut+1, 2*cut+1));

        m.data = m.data[abs(m.nxstart)-cut:abs(m.nxstart)+cut, abs(m.nystart)-cut:abs(m.nystart)+cut, abs(m.nzstart)-cut:abs(m.nzstart)+cut];

        spacing = m.origin[1] / m.nxstart;
        m.cella = Float32.(spacing * [2*cut+1; 2*cut+1; 2*cut+1]);
        m.nxstart = Int32(sign(m.nxstart) * (cut+1));
        m.nystart = Int32(sign(m.nystart) * (cut+1));
        m.nzstart = Int32(sign(m.nzstart) * (cut+1));
        m.gridposition = m.gridposition[reshape(support, (m.nx*m.ny*m.nz, )), :];
        m.nx, m.ny, m.nz, m.mx, m.my, m.mz = 2*cut+1, 2*cut+1, 2*cut+1, 2*cut+1, 2*cut+1, 2*cut+1;
        m.dmin = Float32(minimum(m.data));
        m.dmean = Float32(mean(m.data));
        m.rms = Float32(std(m.data));
        verbose ? @info("MRC: The input ::mrc was modified (chopped) successfully ...") : nothing;
    else
        verbose ? @info("MRC: The input ::mrc was not modified (chopped) using this cut = $cut ...") : nothing;
    end
end



function mrc_chop(m::mrc, cut::Signed)

    m2 = mrc_copy(m);
    mrc_chop!(m2, cut; verbose=false);
    size(m2.data) == size(m.data) ? @info("MRC: The input ::mrc was not modified (chopped) at this cut = $cut ...") : @info("MRC: Operation completed (chopped) successfully ...");

    return m2;
end




##############################################################################################################
    #     #####   #####  ###  #####  #     #
   # #   #     # #     #  #  #     # ##    #
  #   #  #       #        #  #       # #   #
 #     #  #####   #####   #  #  #### #  #  #
 #######       #       #  #  #     # #   # #
 #     # #     # #     #  #  #     # #    ##
 #     #  #####   #####  ###  #####  #     #
##############################################################################################################

function mrc_assign!(m::mrc, idx::Array{<:Signed, 1}, dd::Array{<:AFloat, 1}; verbose::Bool=true)

    length(idx) != length(dd) ? error("MRC: Diemnsions of indeces and densities do not match !!") : nothing;


    m.data[idx] = Float32.(dd[:]);
    m.dmax = Float32(maximum(m.data));
    m.dmin = Float32(minimum(m.data));
    m.dmean = Float32(mean(m.data));
    m.rms = Float32(std(m.data));

    verbose ? @info("MRC: The input ::mrc was modified (assigned) successfully ...") : nothing;

end



function mrc_assign(m::mrc, idx::Array{<:Signed, 1}, dd::Array{<:AFloat, 1})

    m2 = mrc_copy(m);
    mrc_assign!(m2, idx, dd; verbose=false);
    @info("MRC: Operation completed (assigned) successfully ...");

    return m2;

end


##############################################################################################################
  #####  ####### #     # #     # #######  #####  #######
 #     # #     # ##    # ##    # #       #     #    #
 #       #     # # #   # # #   # #       #          #
 #       #     # #  #  # #  #  # #####   #          #
 #       #     # #   # # #   # # #       #          #
 #     # #     # #    ## #    ## #       #     #    #
  #####  ####### #     # #     # #######  #####     #
##############################################################################################################


function mrc_connect!(m::mrc, connect_cutoff::AFloat, lower_bound::AFloat; verbose::Bool=true);

    d = Float64.(reshape(m.data, (m.nx*m.ny*m.nz, )));
    ndimage = pyimport("scipy.ndimage");
    binary3d = m.data .> connect_cutoff;
    # filter out unconnected voxels
    open_binary3d = ndimage.binary_opening(binary3d);
    close_binary3d = ndimage.binary_closing(open_binary3d);
    keepidx = findall(reshape(close_binary3d, (m.nx*m.ny*m.nz, )));
    discardidx = findall(reshape(.!close_binary3d, (m.nx*m.ny*m.nz, )));
    mrc_assign!(m, keepidx, d[keepidx]; verbose=false);
    mrc_assign!(m, discardidx, fill(lower_bound, length(discardidx)); verbose=false);
    scale = sum(d) / sum(m.data);
    mrc_scale!(m, scale; verbose=false);

    verbose ? @info("MRC: The input ::mrc was modified (connected) successfully ...") : nothing;

    return nothing;

end



function mrc_connect(m::mrc, connect_cutoff::AFloat, lower_bound::AFloat)

    m2 = mrc_copy(m);
    mrc_connect!(m2, connect_cutoff, lower_bound; verbose=false);
    @info("MRC: Operation completed (connected) successfully ...");

    return m2;
end




##############################################################################################################
  #####  ######  ####### #     #        #####  ####### #     # #     # #######  #####  #######
 #     # #     # #        #   #        #     # #     # ##    # ##    # #       #     #    #
 #       #     # #         # #         #       #     # # #   # # #   # #       #          #
 #  #### ######  #####      #    ##### #       #     # #  #  # #  #  # #####   #          #
 #     # #   #   #          #          #       #     # #   # # #   # # #       #          #
 #     # #    #  #          #          #     # #     # #    ## #    ## #       #     #    #
  #####  #     # #######    #           #####  ####### #     # #     # #######  #####     #
##############################################################################################################

function mrc_greyconnect!(m::mrc, size::Signed; verbose=true)

    ndimage = pyimport("scipy.ndimage");
    #grey = ndimage.grey_opening(m.data, size=(size, size, size));
    grey = ndimage.grey_closing(m.data, size=(size, size, size));

    scale = sum(m.data) / sum(grey);
    m.dmin = Float32(scale * minimum(grey));
    m.dmax = Float32(scale * maximum(grey));
    m.dmean = Float32(scale * mean(grey));
    m.rms = Float32(scale * std(grey));
    m.data = Float32.(scale * grey);

    verbose ? @info("MRC: The input ::mrc was modified (grey_connected) successfully ...") : nothing;

    return nothing;
end




function mrc_greyconnect(m::mrc, size::Signed)

    m2 = mrc_copy(m);
    mrc_greyconnect!(m2, size; verbose=false);
    @info("MRC: Operation completed (grey connected) successfully ...");

    return m2;
end



##############################################################################################################
  #####  ####### #     # #     # ####### #       #     # #######
 #     # #     # ##    # #     # #     # #       #     # #
 #       #     # # #   # #     # #     # #       #     # #
 #       #     # #  #  # #     # #     # #       #     # #####
 #       #     # #   # #  #   #  #     # #        #   #  #
 #     # #     # #    ##   # #   #     # #         # #   #
  #####  ####### #     #    #    ####### #######    #    #######
##############################################################################################################

function mrc_convolve!(m::mrc, size::AFloat; mode::String="nearest", verbose::Bool=true)

    ndimage = pyimport("scipy.ndimage");
    spacing = m.cella[1] / m.nx;
    filter_size = ceil(size / spacing);
    zoom = spacing * filter_size / size;
    mrc_zoom!(m, zoom; verbose=false, normalize=true);
    n = Int64(filter_size);
    conv = ndimage.convolve(m.data, ones(n, n, n), mode=mode);

    m.dmin = Float32(minimum(conv));
    m.dmax = Float32(maximum(conv));
    m.dmean = Float32(mean(conv));
    m.rms = Float32(std(conv));
    m.data = Float32.(conv);

    verbose ? @info("MRC: The input ::mrc was modified (convolved) successfully ...") : nothing;

    return nothing;
end


function mrc_convolve(m::mrc, size::AFloat; mode::String="nearest")

    m2 = mrc_copy(m);
    mrc_convolve!(m, size; mode=mode, verbose=false);
    @info("MRC: Operation completed (convolved) successfully ...");

    return m2;
end











# more functionalities to come



end
