#/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -q long 
#PBS -N test_PIMD
#PBS -j oe

cd $PBS_O_WORKDIR

./PIMD.x inputmonomer.inp > outmonomer.out 
