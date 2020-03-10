# SWAXS.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jipq6175.github.io/SWAXS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jipq6175.github.io/SWAXS.jl/dev)
[![Build Status](https://travis-ci.com/jipq6175/SWAXS.jl.svg?branch=master)](https://travis-ci.com/jipq6175/SWAXS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jipq6175/SWAXS.jl?svg=true)](https://ci.appveyor.com/project/jipq6175/SWAXS-jl)
[![Codecov](https://codecov.io/gh/jipq6175/SWAXS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jipq6175/SWAXS.jl)
[![Coveralls](https://coveralls.io/repos/github/jipq6175/SWAXS.jl/badge.svg?branch=master)](https://coveralls.io/github/jipq6175/SWAXS.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/jipq6175/SWAXS.jl.svg)](https://cirrus-ci.com/github/jipq6175/SWAXS.jl)


Compute the solution x-ray scattering profiles using orientational average

Possible compatible file formats:
1. Single `.pdb` file
2. A pair of `.pdb` files: `solute.pdb` and `solvent.pdb` (used for buffer subtraction and excluded volume)
3. 3D electron density map: `.mrc` or `.map`. Buffer subtraction is implemented.
4. 3D shape with dummy voxels: `.binvox`. Buffer subtraction is implemented.


## Required Packages

1. ArgParse.jl
2. DelimitedFiles.jl
3. LinearAlgebra.jl
4. Statistics.jl
5. Reexport.jl


## Usage

Run `julia swaxs.jl --help`.

```

usage: swaxs.jl [-D DENSITY] [-P PDB] [-T SOLUTE] [-V SOLVENT]
                [-B BINVOX] [-s SOLVENT_DENSITY] [-v VOXEL_DENSITY]
                [-c DENSITY_CUTOFF] [-J J] [-n NPR] [-o OUTPUT] [-h]
                qmin qspacing qmax

positional arguments:
  qmin                  Starting q value, in (1/A) (type: Float64)
  qspacing              q grid spacing, in (1/A) (type: Float64)
  qmax                  Ending q value, in (1/A) (type: Float64)

optional arguments:
  -D, --density DENSITY
                        Electron density map file: .mrc or .map
  -P, --pdb PDB         Single file containing atomic coordinates: .pdb
  -T, --solute SOLUTE   The solute.pdb file
  -V, --solvent SOLVENT
                        The solvent.pdb file
  -B, --binvox BINVOX   The shape file of dummy voxels: .binvox
  -s, --solvent_density SOLVENT_DENSITY
                        The bulk solvent electron density in e/A^3 (type: Float64, default: 0.335)
  -v, --voxel_density VOXEL_DENSITY
                        The averaged electron density on dummy voxels in e/A^3 (type: Float64, default: 0.5)
  -c, --density_cutoff DENSITY_CUTOFF
                        The electron density cutoff to exclude near-zero voxels (type: Float64, default: 0.001)
  -J, --J J             Number of orientations to be averaged (type: Int64, default: 1200)
  -n, --npr NPR         Number of parallel workers for computation (type: Int64, default: 8)
  -o, --output OUTPUT   Output file prefix for saving the .dat file (default: "output")
  -h, --help            show this help message and exit

```


## Examples

### 1. Calculating from single `.pdb` file containing atomic-coordinates: `rna.pdb`

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


### 2. Calculating from a pair of `.pdb` files to account for buffer subtraction: `solute.pdb` and `solvent.pdb`

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



### 3. Calculating from CCP4 3D electron density map `.mrc` or `.map` files: `den.mrc`

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




### 4. Calculating from 3D shape with dummy voxels: `cat.binvox`

   The `cat.binvox` contains 51x51x51 volumetric data with a grid size of `2A` and is considered to be the 3D cat shape with Excessive Density `0.5 e/A^3`.
   Run

```
julia -O 3 --color=yes swaxs.jl --binvox cat.binvox -s 0.335 -v 0.5 -J 1500 0.0 0.01 1.0 -o test
```

   The output should look like

```
[ Info: SWAXS: Setting up parallel workers ...
[ Info: Please wait ...
[ Info: Computing SWAXS (J=1500) using shape file: cat.binvox, with voxel_density=0.5, sden=0.335.
[ Info: SWAXS program completed successfully: elapsed time = 13.130944699 seconds with 8 cores.
```

   And the swaxs profile from `q = 0.0` to `q = 1.0` with spacing `0.01` is saved as `test.dat`.



## Notes

1. The options `--pdb`, `--density`, `--binvox` and `--solute --solvent` cannot not be specified at the same time. Otherwise, error will be thrown.
2. If high-throughput computation is required, one should bypassing the commandline because it sets up parallel workers everytime. To avoid that, set up your parallel workers and call `@everywhere include("swaxs_module.jl")` and `@everywhere using .SWAXS` in your julia script.
3. If your structures contain too many atoms or voxels, I suggest you to take a look at GPUSWAXS.
4. The atmo-name processor might be a little bit buggy, I will try to fix some of the parsing issues.
5. Not sure if I will maintain this or not in the future ...
