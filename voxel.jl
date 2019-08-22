# Voxel
mutable struct Voxel
    version::Int64;
    dim::Int64;
    translate::Array{AFloat, 1};
    scale::AFloat;
    rawdata::Array{Bool, 1};
    coordata::Array{AFloat, 2};
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
