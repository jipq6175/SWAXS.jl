# Voxel
mutable struct Voxel
    version::Int64;
    dim::Int64;
    translate::Vector{Float64};
    scale::Float64;
    rawdata::Vector{Bool};
    coordata::Matrix{Float64};
end


function readvox(filename::String)

    # initialize the voxel struct
    vox = Voxel(1, 1, zeros(3), 0.0, fill(true, 2), zeros(2,2));

    # open and read in the headers
    f = open(filename, "r");

    # first line
    line = readline(f);
    startswith(line, "#binvox")? vox.version = parse(Int64, split(line)[2]): error("Unrecognized .binvox file with line 1: $line.");

    # second line
    line = readline(f);
    if startswith(line, "dim")
        p = parse.(Int64, split(line)[2:4]);
        length(unique(p)) == 1? vox.dim = p[1]: error("Dimensions are different on x, y, z with $p.");
    else
        error("Unrecognized .binvox file with line 2: $line.");
    end

    # Third line
    line = readline(f);
    startswith(line, "translate")? vox.translate = parse.(Float64, split(line)[2:4]): error("Unrecognized .binvox file with line 3: $line.");

    # Fourth line
    line = readline(f);
    startswith(line, "scale")? vox.scale = parse(Float64, split(line)[2]): error("Unrecognized .binvox file with line 4: $line.");

    # Fifth line and ready to read data
    line = readline(f);
    line != "data"? error("Unrecognized .binvox file with line 5: $line."): nothing; #info("Headers all good.. Ready to read in binary data ...");

    b = UInt8[];
    readbytes!(f, b, Inf);
    close(f);

    n = Int(length(b)/2);

    b = Int.(reshape(b, 2, n));

    m = 1;
    vox.rawdata = fill(false, vox.dim^3);

    for i = 1:n
        k = b[2, i];
        vox.rawdata[m:m+k-1] = fill(Bool(b[1, i]), k);
        m = m + k;
    end

    grid3d = zeros(vox.dim^3, 3);

    m = 1;
    for i = 1:vox.dim
        for j = 1:vox.dim
            for k = 1:vox.dim
                grid3d[m, :] = 2*[k-1 j-1 i-1];
                m = m + 1;
            end
        end
    end

    vox.coordata = grid3d[vox.rawdata, :];
    return vox;
end



# Voxelize the .obj files or .off files within all the files in a directory
function voxelizer(dir::String)

    filelist = Vector{String}(0);

    for (root, dirs, files) in walkdir(dir)
        info("Getting files in $root ...");
        for file in files
            push!(filelist, joinpath(root, file));
        end
    end

    filelist = filelist[endswith.(filelist, ".obj") .| endswith.(filelist, ".off")];

    info("Starting to voxelize...");
    n = length(filelist);
    for i = 1:n

        command = `binvox -d 51 -cb $(filelist[i])`;
        try
            tic();
            run(pipeline(command, stdout="out.log", append=true));
            t = toq();
            @printf("%.4f%% (%06d / %06d) \n%s Success, with elapsed time = %.3f s.\n", 100*i/n, i, n, filelist[i], t);
            # println("($i/$n): $(filelist[i]) Success, with elapsed time = $t s.");
        catch err
            warn(@sprintf("%.4f%% (%06d / %06d) \n%s failed. Skipped!!", 100*i/n, i, n, filelist[i]));
        end

    end

end
