#!/usr/bin/bash

#PBS -o run_qsub_long.out
#PBS -e run_qsub_long.err
#PBS -l nodes=1:ppn=16
#PBS -q test

eval "$(conda shell.bash hook)"
conda activate pygetm
cd $HOME/BLUE2/medsea
NPROCS=`wc -l < $PBS_NODEFILE`
mpiexec -np $NPROCS python -u $HOME/BLUE2/medsea/medsea.py

