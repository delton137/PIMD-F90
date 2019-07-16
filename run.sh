#!/bin/bash -l
##PBS script for running on a cluster with slurm
#SBATCH -p debug
#SBATCH --job-name=testPIMD
#SBATCH --comment=testPIMD
#SBATCH --time=1:00:00
#SBATCH --nodes=1 
##SBATCH --ntasks-per-node=80
#SBATCH --ntasks=80
####SBATCH --exclude=node[00-03]

module purge
module load mpi/openmpi/2.1.0


cd $SLURM_SUBMIT_DIR
./PIMD.x inputmonomer > outmonomer
#mpirun -np 32 ./PIMD.x input32.inp  > out32
