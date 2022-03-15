#!/bin/bash

#SBATCH -A ucb62_summit3     # Summit allocation
#SBATCH --partition=shas     # Summit partition
#SBATCH --qos=normal                 # Summit qos
#SBATCH --time=001:25:00           # Max wall time
#SBATCH --nodes=1            # Number of Nodes

#SBATCH --ntasks=1           # Number of tasks per job
#SBATCH --cpus-per-task=20

#SBATCH --job-name=snap_insar_array        # Job submission name
#SBATCH --output=snap_insar_array.%A_%a.out   # Output file name with Job ID

#SBATCH --mail-type=END            # Email user when job finishes
#SBATCH --mail-user=clwo4142@colorado.edu # Email address of user

#SBATCH --array=1-9

module purge
module load singularity

working_directory=/scratch/summit/clwo4142/SNAP_RUNS/SNAP_CONTAINER_TEST_1/Path_114_Frame_76_jvm_tests
cd $working_directory

/projects/clwo4142/software/anaconda/envs/mycustomenv/bin/python3.6 /scratch/summit/clwo4142/SNAP_RUNS/SNAP_CONTAINER_TEST_1/Path_114_Frame_76_jvm_tests/Main/Run_Snap_Steps.py $SLURM_ARRAY_TASK_ID


