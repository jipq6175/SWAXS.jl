
using Distributions, LinearAlgebra, Printf

function random_water(file::String, s::Float64; output::String="output")

    fin = open(file, "r");
    lines = readlines(fin);
    close(fin);
    nlines = length(lines);

    # Do a quick scan
    nwaters = 0;
    waterstart = 1;
    for i = 1:nlines

        substr = split(lines[i]);
        if length(substr) < 4
            continue;
        else
            (substr[3] == "OW") && (substr[4] == "SOL") ? nwaters += 1 : nothing;
            #nwaters == 1 ? waterstart = i : nothing;
        end
    end
    # println(waterstart);
    # pre-define normal distribution and
    Σ = s^2 * Matrix{Float64}(I, 3, 3);
    gaussian3 = MvNormal(Σ);
    displacement = round.(rand(gaussian3, nwaters), digits=3);

    for i = 1:nwaters

        targetline = waterstart + 3*(i - 1);

        for j = 0:2

            coor = parse.(Float64, split(lines[targetline + j][31:56]));
            coor = coor + displacement[:, i];
            lines[targetline + j] = @sprintf "%s %7.3f %7.3f %7.3f %s" lines[targetline + j][1:30] coor[1] coor[2] coor[3] lines[targetline + j][57:end];
        end

    end

    fout = open(output * ".pdb", "w");
    for i = 1:nlines
        write(fout, lines[i]);
        write(fout, "\n");
    end
    close(fout);

    return nothing;
end
