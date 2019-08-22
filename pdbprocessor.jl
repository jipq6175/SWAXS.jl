

"""
    m = pdbprocessor(filename[, atomnames[, waters=true[, ions=true]]])

Process the atomic coordinates `filename` with desired inclusion of `atomnames` and options to include `waters` or `ions`.

`m[:, 1]` is the atomic scattering factors at Energy = 9.91555KeV and `m[:, 2:4]` are the `(x, y, z)` coordinates of the atom.

"""
function pdbprocessor(filename::String, atomnames::Vector{String}=["H", "C", "N", "O", "NA", "MG", "P", "CL", "K"]; waters::Bool=true, ions::Bool=true)

    # atomic scattering factors, at Energy = 9.91555KeV
    const ASF = Dict([("H", 0.999980), ("C", 6.01286), ("N", 7.02218), ("O", 8.03637), ("NA",11.0999), ("MG", 12.1329), ("P", 15.2382), ("CL", 17.3078), ("K+", 19.3598)]);

    # Some residues to worry about
    const IONS = ["CIP"; "K+"; "NA+"];
    const SOL = ["SOL"; "WAT"; "HOH"];

    # Read pdb file
    f = open(filename, "r");
    lines = readlines(f);
    close(f);

    atoms = lines[startswith.(lines, "ATOM")];
    n = length(atoms);

    # trying to locate the coordinates
    sp = split(atoms[1]);
    index = collect(1:length(sp));
    id = index[typeof.(parse.(sp)) .== Float64][1];

    # forming the debye matrix
    mat = Matrix{Float64}(0, 4);
    for i = 1:n
        sp = convert.(String, split(atoms[i]));
        # Check waters
        if (sp[4] in IONS)
            ions ? mat = vcat(mat, [ASF[sp[3]] parse(sp[id]) parse(sp[id+1]) parse(sp[id+2])]) : continue;
        elseif (sp[4] in SOL)
            waters ? mat = vcat(mat, [ASF[atom_identify(sp[3])] parse(sp[id]) parse(sp[id+1]) parse(sp[id+2])]) : continue;
        else
            mat = vcat(mat, [ASF[atom_identify(sp[3])] parse(sp[id]) parse(sp[id+1]) parse(sp[id+2])]);
        end
    end
    return mat;
end


function atom_identify(atom::String)
    const ATOM = ["H", "C", "N", "O", "NA", "MG", "P", "CL", "K"];
    at = split(atom, "");
    id = find([x in ATOM for x in at])[1];
    return at[id];
end
