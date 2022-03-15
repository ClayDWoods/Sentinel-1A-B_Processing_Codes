#!/bin/bash

#SBATCH -A ucb62_summit3     # Summit allocation
#SBATCH --partition=shas     # Summit partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=000:05:00           # Max wall time
#SBATCH --nodes=1            # Number of Nodes

#SBATCH --ntasks=1           # Number of tasks per job
#SBATCH --cpus-per-task=4

#SBATCH --job-name=max_curv        # Job submission name
#SBATCH --output=max_curv.%A_%a.out   # Output file name with Job ID

#SBATCH --mail-type=END            # Email user when job finishes
#SBATCH --mail-user=clwo4142@colorado.edu # Email address of user

working_directory=/scratch/summit/clwo4142/SNAP_RUNS/Snap_Runs/Pairs_4_msbasv3/Pair_P114_F76_P150_F510/MSBASv3_Timeseries
cd $working_directory

/projects/clwo4142/software/anaconda/envs/mycustomenv/bin/python3.6 ./Main/Max_curvature_extract.py
