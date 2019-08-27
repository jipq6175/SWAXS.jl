# SWAXS.jl

Compute the solution x-ray scattering profiles using orientational average

Possible compatible file formats:
1. Single `.pdb` file
2. A pair of `.pdb` files: `solute.pdb` and `solvent.pdb` (used for buffer subtraction and excluded volume)
3. 3D electron density map; `.mrc` or `.map`. Buffer subtraction will be implemented.


## Usage

Run `julia swaxs.jl --help`.

```

usage: swaxs.jl [-d DENSITY] [-p PDB] [-t SOLUTE] [-v SOLVENT]
                [-s SOLVENT_DENSITY] [-c DENSITY_CUTOFF] [-J J]
                [-n NPR] [-h] qmin qspacing qmax

positional arguments:
  qmin                  Starting q value, in (1/A) (type: Float64)
  qspacing              q grid spacing, in (1/A) (type: Float64)
  qmax                  Ending q value, in (1/A) (type: Float64)

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


## Examples

1. Calculating from single `.pdb` file containing atomic-coordinates: `rna.pdb`

   The `rna.pdb` file contains about 10000 atoms.
   Run

```
julia -O 3 --color=yes swaxs.jl --pdb rna.pdb -J 1500 0.0 0.01 1.0 -o test
```

  The output should look like

```
[ Info: SWAXS: Setting up parallel workers ...
[ Info: Please wait ...
[ Info: Computing SWAXS (J=1500) using single PDB file: rna.pdb.
[ Info: SWAXS program completed successfully: elapsed time = 20.375546201 seconds with 8 cores.
```

   And the swaxs profile from `q = 0.0` to `q = 1.0` with spacing `0.01` is saved as `test.dat`.


2. Calculating from a pair of `.pdb` files to account for buffer subtraction: `solute.pdb` and `solvent.pdb`

   Each of the `.pdb` file contains about 9000 atoms.
   Run
```
julia -O 3 --color=yes swaxs.jl --solute solute.pdb --solvent solvent.pdb -J 1500 0.0 0.01 1.0 -o test
```

  The output should look like the following

```
[ Info: SWAXS: Setting up parallel workers ...
[ Info: Please wait ...
[ Info: Computing SWAXS (J=1500) using solute: solute.pdb and solvent: solvent.pdb.
[ Info: SWAXS program completed successfully: elapsed time = 31.219024699 seconds with 8 cores.
```

  And the swaxs profile from `q = 0.0` to `q = 1.0` with spacing `0.01` is saved as `test.dat`.



3. Calculating from CCP4 3D electron density map `.mrc` or `.map` files: `den.mrc`

   The `den.mrc` contains 96x96x96 volumetric data and is considered to be the Excessive Density
   Run

```
julia -O 3 --color=yes swaxs.jl --density den.mrc -s 0.335 -c 0.001 -J 1500 0.0 0.01 1.0 -o test
```

   The output should look like

```
[ Info: SWAXS: Setting up parallel workers ...
[ Info: Please wait ...
[ Info: Computing SWAXS (J=1500) using electron density file: den.mrc, with sden=0.335 cutoff=0.001.
[ Info: SWAXS program completed successfully: elapsed time = 15.9209419 seconds with 8 cores.
```

   And the swaxs profile from `q = 0.0` to `q = 1.0` with spacing `0.01` is saved as `test.dat`.



## Note

1. The options `--pdb`, `--density` and `--solute --solvent` cannot not be specified at the same time.
2. If high-throughput computation is required, one should bypassing the commandline because it sets up parallel workers everytime. To avoid that, set up your parallel workers and call `@everywhere include("swaxs_module.jl")` and `@everywhere using .SWAXS` in your julia script.
3. If your structures contain too many atoms or voxels, I suggest you to take a look at GPUSWAXS. 
4. Not sure if I will maintain this or not.
