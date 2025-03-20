#!/bin/bash
#SBATCH --job-name=det_ser       # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)

module purge
./PeleC3d.gnu.ex det.inp