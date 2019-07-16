# Fast and accurate *ab-initio* Path Integral Molecular Dynamics in Fortran 90

2015-2017 Daniel C. Elton

PIMD-F90 implements path integral molecular dynamics (PIMD) to capture the nuclear quantum effects in water.

# Current status:

'''The code is currently broken! Those interested doing PIMD are encouraged to look at the [i-PI Python package](http://ipi-code.org/) instead.'''

We received a request to remove the TTM / SPC force calculation routines. Those routines were publicly available on [Google code](https://code.google.com/archive/p/nqcbabel/) in the "nqcbabel" repository under /nqcbabel/trunk/extra/ttm-ice and clearly had a [GNU GPL v3](https://www.gnu.org/licenses/gpl.html) license associated with them. We were told this licensing was a mistake however and agreed to comply with the request to remove them.  

The code should still work with SIESTA but needs to be fixed to remove references to the TTM and SPC forcefield routines:
(potential.f90, consts.f90, find_neigh.f90, pot_mod.f90, pot_ttm.f90, smear.f90 math.f90 ewald.f90 nasa.f90 nasa_mod.f90)

TTM3-F ref: George S. Fanourgakis and Sotiris S. Xantheas, J. Chem. Phys. 128, 074506 (2008)  
TTM2-F ref: J. Phys. Chem. A. vol 110, page 4100-4106, (2006)  

## License information

All files written by D. C. Elton include the MIT open source license in their header.

Additional files taken from open source codebases are as follows:

### fsiesta_pipes.f90
This file is taken from the [SIESTA package](https://departments.icmab.es/leem/siesta/) which is licensed under the GNU GPL v3 open source license ([link](https://launchpad.net/siesta)).

## SIESTA
The PIMD code can obtain forces & energies from the density functional theory package [SIESTA](http://departments.icmab.es/leem/siesta/).
It uses unix pipes to communicate with SIESTA. A working flush() routine is required for this (see pflush.f90).


## Features
* Parallelization over the beads with MPI
* Options for a global Nose-Hoover thermostat or Anderson barostat
* virial energy and pressure estimators
* many output options
* Choice of Langevin or Nose-Hoover bead thermostats
* Option to not thermostat the centroid
* Option to perform the PIMD on the monomer only (a subcase of the ring polymer contraction scheme ref: [JCP 129, 024105 (2008)](http://dx.doi.org/10.1063/1.2953308) )
* implementation of custom PES for using the "monomer PIMD" method with DFT
* By changing the adiabaticity parameter (ficticious mass rescaling) the code can also implement Centroid Molecular Dynamics (CMD) and Ring Polymer Molecular Dynamics (RPMD).


## Compilation
currently compiles with openMPI (gcc) mpif90 compiler or Intel mpiifort. However SIESTA functionality doesnt work with mpiifort

other notes: i
compiles with openmpi/gcc/64/1.8.5 but not mvapich/gcc/64/1.2rc1

SIESTA compiles with mvapich but not openmpi


## Input file specification
See the example input file input128.inp for program options

See [manual/PIMD_manual.pdf](manual/PIMD_manual.pdf)


## Acknowledging

If you use this code, please cite our publication:

D. C. Elton, M. Fritz, and M.-V Fernández-Serra. "[Using a monomer potential energy surface to perform approximate path integral molecular dynamics simulation of ab-initio water at near-zero added cost](https://pubs.rsc.org/en/Content/ArticleLanding/2019/CP/C8CP06077K#!divAbstract)", *Phys. Chem. Chem. Phys.*, **21**, 409 (2019) [[arXiv](https://arxiv.org/abs/1803.05740)]

Here is a bibTex entry for it:

@article{Elton2019,  
  doi = {10.1039/c8cp06077k},  
  url = {https://doi.org/10.1039/c8cp06077k},  
  year  = {2019},  
  publisher = {Royal Society of Chemistry ({RSC})},  
  volume = {21},  
  number = {1},  
  pages = {409--417},  
  author = {Daniel C. Elton and Michelle Fritz and Marivi Fern{\'{a}}ndez-Serra},  
  title = {Using a monomer potential energy surface to perform approximate path integral molecular dynamics simulation of ab initio water at near-zero added cost},  
  journal = {Physical Chemistry Chemical Physics}  
}
