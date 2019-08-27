# SWAXS.jl

Compute the solution x-ray scattering profiles using orientational average

Possible compatible file formats:
1. Single .pdb file
2. A pair of .pdb files: solute.pdb and solvent.pdb (used for buffer subtraction and excluded volume)
3. 3D electron density map; .mrc or .map. Buffer subtraction will be implemented.


## Usage

```

usage: swaxs.jl [-d DENSITY] [-p PDB] [-t SOLUTE] [-v SOLVENT]
                [-s SOLVENT_DENSITY] [-c DENSITY_CUTOFF] [-J J]
                [-n NPR] [-h] qmin qspacing qmax

positional arguments:
  qmin                  Starting q value, in (1/A) (type: Float64)
  qmax                  Ending q value, in (1/A) (type: Float64)
  qspacing              q grid spacing, in (1/A) (type: Float64)

optional arguments:
  -d, --density DENSITY Electron density map file: .mrc or .map
  -p, --pdb PDB         Single file containing atomic coordinates: .pdb
  -t, --solute SOLUTE   The solute.pdb file
  -v, --solvent SOLVENT The solvent.pdb file
  -s, --solvent_density SOLVENT_DENSITY
                        The bulk solvent electron density in e/A^3 (type: Float64, default: 0.335)
  -c, --density_cutoff DENSITY_CUTOFF
                        The electron density cutoff to exclude near-zero voxels (type: Float64, default: 0.001)
  -J, --J J             Number of orientations to be averaged (type: Int64, default: 1200)
  -n, --npr NPR         Number of parallel workers for computation (type: Int64, default: 8)
  -h, --help            show this help message and exit

```
