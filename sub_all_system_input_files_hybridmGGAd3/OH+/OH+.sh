#!/bin/bash
#SBATCH --mem-per-cpu=2G             # memory, roughly 2 times %mem defined in the input name.com file
#SBATCH --job-name=OH+
#SBATCH --time=0-03:00           # time (DD-HH:MM)
#SBATCH --account=def-yawang
#SBATCH --ntasks=4               # number of MPI processes
module load StdEnv/2020
mpirun /home/ysy1111/nwchem-origin/nwchem-7.0.2-release/bin/LINUX64/nwchem OH+.nw > OH+.out 2>OH+.log