
######################################
## pdb processor ver 2.0 to extract ##
######################################
"""
    atomid, mat = SimplyPDB(filename::String, atomnames::Vector{String}; waters=true, ions=true)

Process the atomic coordinates `filename` with desired inclusion of `atomnames` and options to include `waters` or `ions`.

`atomid` is the `Vector{String}` that inlcude every elements in the `filename`.
`mat` is the `(x, y, z)` coordinates of the corresponding atom.

# Benchmarks

For a pdb file that contains 7318 atoms, including solvents and ions.
```julia-repl
julia> @time SimplyPDB(solutefile)
1.843963 seconds (821.26 k allocations: 647.427 MiB, 49.35% gc time)
(>String[7318], 7318x3 Array{Float64,2}:)
```

"""
function SimplyPDB(filename::AS, atomnames::Array{String, 1}=["H", "C", "N", "O", "NA", "MG", "P", "CL", "K"]; waters::Bool=true, ions::Bool=true)

    # Some residues to worry about
    IONS = ["CIP"; "K+"; "NA+"; "MG+"];
    SOL = ["SOL"; "WAT"; "HOH"];

    # Read pdb file
    f = open(filename, "r");
    lines = readlines(f);
    close(f);

    atoms = lines[startswith.(lines, "ATOM")];
    n = length(atoms);

    # trying to locate the coordinates
    sp = split(atoms[1]);
    index = collect(1:length(sp));
    id = index[typeof.(parse.(Float64, sp)) .== Float64][1];

    # forming the debye matrix
    mat = Array{Float64, 2}(undef, 0, 3);
    atomid = Array{String, 1}(undef, 0);
    for i = 1:n
        sp = convert.(String, split(atoms[i]));
        # Check waters
        if (sp[4] in IONS)
            if ions
                push!(atomid, sp[3]);
                mat = vcat(mat, parse.(sp[id:id+2])');
            end
        elseif (sp[4] in SOL)
            if waters
                push!(atomid, "SOL-" * atom_identify(sp[3]));
                mat = vcat(mat, parse.(sp[id:id+2])');
            end
        else
            push!(atomid, atom_identify(sp[3]));
            # Avoid the "RX N" residues generated from the
            typeof(parse(sp[id])) != Int64 ? mat = vcat(mat, parse.(sp[id:id+2])') : mat = vcat(mat, parse.(sp[id+1:id+3])');
        end
    end
    return atomid, mat;
end


# Identify atom names
"""
    atom = atom_identify(atomname);

Converting a confusing `atomname` in a PDB file to recognized `atom`.

# Examples

```julia-repl
julia> atom_identify("OP2")
"o"
julia> atom_identify("2H5'")
"H"

```

"""
function atom_identify(atom::String)
    ATOM = ["H", "C", "N", "O", "NA", "MG", "P", "CL", "K"];
    at = split(atom, "");
    id = find([x in ATOM for x in at])[1];
    return convert(String, at[id]);
end




###############################################################
## convert atomid to form factor given a scattering vector q ##
###############################################################
"""
    atomaff = AtomAFF(atomid::Vector{String}, q::Float64; coef_file);

Compute the atomic form factors of `atomid` at specific `q::Float64` value. This computation has to be carried out at every `q`. The calculation is implemented using the following formula:
``
f(q) = c + \\sum_{k=1}^4 a_k\\exp\\left[ - b_k\\left( \\frac{q}{4\\pi} \\right)^2\\right]
``
The coefficients ``a_k, b_k and c`` are tabulated in `coef_file="aff.dat"` for common elements in nucleic acids.

The oxygen and hydrogen in water molecules have electro-negativity corrections as follows.

Oxygens
``
f'(O) = f(O)\\left[ 1 + 0.12\\exp\\left( -\\frac{1}{2} \\left( \\frac{q}{\\delta} \\right)^2 \\right) \\right]
``

Hydrogens:
``
f'(H) = f(H)\\left[ 1 - 0.48\\exp\\left( -\\frac{1}{2} \\left( \\frac{q}{\\delta} \\right)^2 \\right) \\right]
``

where ``\\delta = 2.2`` is the published value.

# Examples

```julia-repl
julia> atomid, mat = SimplyPDB(solutefile);
julia> @time AtomAFF(atomid, 0.3)
0.008292 seconds (17.01 k allocations: 468.625 KiB)
> Float64[7318]
```

# Notes
1. Since everytime the function `AtomAFF()` is executed, it reads in the file `"aff.dat"`, at a fixed `q`, it should run for only one time.

"""
function AtomAFF(atomid::Vector{String}, q::T; coef_file::String="aff.dat") where T<:Real

    # Load the a, b, c constants for atoms and ions form factor
    # Using the data: aff.dat
    f = open(coef_file, "r");
    lines = readlines(f);
    close(f);

    lines = lines[.~startswith.(lines, "##")];
    n = length(lines);
    ASF = Dict("String" => -1.0);

    s = q/(4*pi);

    for i = 1:n
        strings = split(lines[i], "\t");
        atom = strings[1];
        coefs = parse.(strings[3:11]);
        a, b, c = coefs[[1; 3; 5; 7]], coefs[[2; 4; 6; 8]], coefs[9];
        fq = c + a' * exp.(-b * s^2);
        merge!(ASF, Dict(atom => fq));
    end

    # Deal with the water O and H
    fo = ASF["O"] * (1. + 0.12 * exp(-0.5*(q/2.2)^2));
    fh = ASF["H"] * (1. - 0.48 * exp(-0.5*(q/2.2)^2));
    merge!(ASF, Dict("SOL-O" => fo, "SOL-H" => fh));

    # Phosphate screening
    # ASF["P"] = 0.25*ASF["P"];

    # convert all the atoms to the ASFs
    m = length(atomid);
    atomaff = Vector{Float64}(m);
    [atomaff[i] = ASF[atomid[i]] for i = 1:m];

    return atomaff;
end


# # Calculate 3D averaged scattering amplitude at one fixed q::Float64
function _DQ(sysA::Tuple{Array{String,1},Array{Float64,2}}, sysB::Tuple{Array{String,1},Array{Float64,2}}, q::T; J::Int64=1500) where T<:Real

    # Prior procedure
    # Given one q value, we only need to use AtomAFF() twice, for solute & solvant
    atomidA, pdbmatA = sysA;
    atomidB, pdbmatB = sysB;
    atomaffA = AtomAFF(atomidA, q);
    atomaffB = AtomAFF(atomidB, q);

    # Generate J q vectors
    x = [(2*j - 1 - J)/J for j = 1:J];
    theta, phi = acos.(x), sqrt(pi*J)*asin.(x);
    qmat = q * [sin.(theta) .* cos.(phi) sin.(theta) .* sin.(phi) cos.(theta)];

    # Calculate A, B, D
    qrA = pdbmatA * qmat';
    qrB = pdbmatB * qmat';
    A = [atomaffA' * cos.(qrA); -atomaffA' * sin.(qrA)];
    B = [atomaffB' * cos.(qrB); -atomaffB' * sin.(qrB)];

    # Form factor difference across J q vectors
    d = A .- B;

    # Find the average of the |A-B|^2
    D = mean(sum(d.^2, dims=1));

    return D;
end


######################################################################
## Evaluating Scattering Profile using Parallel computation         ##
## Note that here we only using one single PDB files for evaluation ##
######################################################################
"""
    intensity = SWAXS(sysA, sysB, q::Vector{Float64}; J=1500, npr=7)

Accurately compute the solution X-ray scattering profile of ONE macromolecular state, `sysA` and the corresponding solvent conditions, `sysB` using the Park et al. formalism.

`sysA` and `sysB` are the direct outputs of `SimplyPDB()` on a `"solutefile.pdb"` and `"solventfile.pdb"`.

`J` is the number of uniformly distributed ``q_j`` vectors in 3D, spherically.

`npr` is the number of parallel workers to evaluate `intensity` independently.

Park et al. formalism is as follows:

Fourier transform of system A and B:
``
A(\\mathbf{q}) = \\sum_{i=1}^{N_A} f_i(\\mathbf{q}) \\exp\\left( -i\\mathbf{q}\\mathbf{r} \\right)
``

``
B(\\mathbf{q}) = \\sum_{i=1}^{N_B} f_i(\\mathbf{q}) \\exp\\left( -i\\mathbf{q}\\mathbf{r} \\right)
``

The theoretical intensity ``I(\\mathbf{q})`` at one specific ``\\mathbf{q}`` is
``
D(\\mathbf{q}) = \\left| A(\\mathbf{q}) - B(\\mathbf{q}) \\right|^2
``

The experimentally measured intensity is the average of ``D(\\mathbf{q})`` on the solid angle
``
I(q) = \\frac{1}{4\\pi}\\int D(q)\\sin{\\theta}d\\theta d\\phi \\approx \\frac{1}{J} \\sum_{j=1}^J D(\\mathbf{q_j})
``

# Benchmarks

```julia-repl
julia> q = collect(0.0:0.01:0.95);
julia> @time SWAXS(SimplyPDB(solutefile), SimplyPDB(solventfile), q)
29.056434 seconds (2.45 M allocations: 1.316 GiB, 1.57% gc time)
> Float64[96]
```

"""
function PDBSWAXS(sysA::Tuple{Array{String,1},Array{Float64,2}}, sysB::Tuple{Array{String,1},Array{Float64,2}}, q::Vector{T}; J::Int64=1500, npr::Int64=7) where T<:Real

    # Setting up the parallel environment
    # info("=====   Setting up parallel environment   =====");
    if nprocs() != npr
        nprocs() > npr ? rmprocs(workers()[end-(nprocs()-npr)+1:end]) : addprocs(npr-nprocs());
    else
        nothing;
    end

    # Compute the scattering profile using pmap
    intensity = pmap(x -> _DQ(sysA, sysB, x; J=J), q, distributed=true);

    return intensity;

end



# ######################################################################
# ## Evaluating Scattering Profile using Parallel computation         ##
# ## Use an ENSEMBLE of structures to compute the scattering profiles ##
# ######################################################################
# """
#     intensity = EnSWAXS(soluteseries, solventseries, q::Vector{Float64}; J=1500, npr=7)
#
# Accurately compute the solution X-ray scattering profile of multiple macromolecular states, `soluteseries` and the corresponding solvent conditions, `solventseries` using the Park et al. formalism.
#
# `soluteseries::Vector{String}` is a vector of solute pdbs and `solventseries::Vector{String}` is a vector of solvent pdbs.
#
# `q`, `J` and `npr` follow the same convention of the function `SWAXS()`.
#
# The formalism is exactly the same. The following is how `EnSWAXS()` does the ensemble average.
#
# ``
# I(q) = \\langle \\langle D(\\mathbf{q}) \\rangle_E \\rangle_\\Omega = \\frac{1}{4\\pi}\\int \\left[ \\frac{1}{E}\\sum_{k=1}^E D_k(\\mathbf{q}) \\right]\\sin{\\theta}d\\theta d\\phi
# ``
#
# ``
# I(q) = \\frac{1}{E}\\sum_{k=1}^E \\frac{1}{4\\pi} \\int D_k(\\mathbf{q})\\sin{\\theta}d\\theta d\\phi = \\langle \\langle D(\\mathbf{q}) \\rangle_\\Omega \\rangle_E
# ``
#
# So essentially, the ensemble average is the average of all scattering profiles in the ensemble.
# ``
# I(q) = \\frac{1}{E} \\sum_{k=1}^E I_k(q)
# ``
#
# # Benchmarks
#
# ```julia-repl
# julia> length(soluteseries), length(solventseries)
# (12, 12)
# julia> q = collect(0.0:0.01:0.95);
# julia> @time EnSWAXS(soluteseries, solventseries, q)
# INFO: =====   Deal with the series of PDB files   =====
# INFO: =====   Calculating the Scattering Intensity   =====
# elapsed time: 327.045976352
# INFO: =====   Wrapping up the computations   =====
# > Float64[96]
# ```
#
# # Notes
# 1. At this point, since the `EnSWAXS()` is the average of all scattering profiles from an ensemble, it might be a little bit computationally efficient to use `SWAXS()` and average all the profiles. But the time difference is about 1 second per structure.
#
# """
# function EnSWAXS(soluteseries::Vector{String}, solventseries::Vector{String}, q::Vector{T}; J::Int64=1500, npr::Int64=7) where T<:Real
#
#     # Setting up the parallel environment
#     # info("=====   Setting up parallel environment   =====");
#     if nprocs() != npr
#         nprocs() > npr ? rmprocs(workers()[end-(nprocs()-npr)+1:end]) : addprocs(npr-nprocs());
#     else
#         nothing;
#     end
#
#     # Include self everywhere
#     @everywhere include("scattering_form_factor.jl");
#     @eval @everywhere JJ = $J;
#
#     # Deal wih all the files and coordinates
#     info("=====   Deal with the series of PDB files   =====");
#     sysA = map(x -> SimplyPDB(x), soluteseries);
#     sysB = map(x -> SimplyPDB(x), solventseries);
#     N = length(sysA);
#
#     @eval @everywhere sysA_atomids = $([sysA[i][1] for i = 1:N]);
#     @eval @everywhere sysA_pdbmats = $([sysA[i][2] for i = 1:N]);
#     @eval @everywhere sysB_atomids = $([sysB[i][1] for i = 1:N]);
#     @eval @everywhere sysB_pdbmats = $([sysB[i][2] for i = 1:N]);
#
#
#     # Given one value of q
#     # Create a function for pmap
#     @everywhere function EnDQ(sysA_atomids::Vector{Array{String,1}}, sysA_pdbmats::Vector{Array{T,2}}, sysB_atomids::Vector{Array{String,1}}, sysB_pdbmats::Vector{Array{T,2}}, q::T; J::Int64=JJ) where T<:Real
#
#         # number of ensembles
#         n = length(sysA_atomids);
#         D = zeros(J, n);
#
#         # Convert the atom ids to atom affs
#         sysA_atomaffs = AtomAFF.(sysA_atomids, q);
#         sysB_atomaffs = AtomAFF.(sysB_atomids, q);
#
#         # Create J evenly spaced q vectors
#         x = [(2*j - 1 - J)/J for j = 1:J];
#         theta, phi = acos.(x), sqrt(pi*J)*asin.(x);
#         qmat = q .*[sin.(theta) .* cos.(phi) sin.(theta) .* sin.(phi) cos.(theta)];
#
#         A = zeros(2, J, n);
#         B = zeros(2, J, n);
#
#         for k = 1:n
#
#             # Calculate A, B, D
#             # print(sysA_pdbmats[k]);
#             qrA = sysA_pdbmats[k] * qmat';
#             qrB = sysB_pdbmats[k] * qmat';
#             A = [sysA_atomaffs[k]' * cos.(qrA); -sysA_atomaffs[k]' * sin.(qrA)];
#             B = [sysB_atomaffs[k]' * cos.(qrB); -sysB_atomaffs[k]' * sin.(qrB)];
#
#             # Form factor difference across J q vectors
#             # Note that |<Z>|^2 ~= <|Z|^2> need to account for that
#             d = A - B;
#
#             # Find the average of the |A-B|^2
#             D[:, k] = transpose(sum(d.^2, 1));
#         end
#         # AB = A .- B;
#         # D = mean(sum(A.^2, 1), 3)[1, :, 1] - mean(sum(B.^2, 1), 3)[1, :, 1] + 2 * sum( mean(B, 3) .* mean(AB, 3), 1)[1, :, 1];
#
#         intensity = mean(D);
#         return intensity;
#     end
#
#     info("=====   Calculating the Scattering Intensity   =====");
#     tic();
#     intensity = pmap(x -> EnDQ(sysA_atomids, sysA_pdbmats, sysB_atomids, sysB_pdbmats, x), q);
#     toc();
#
#     info("=====   Wrapping up the computations   =====");
#     return intensity;
#
# end
