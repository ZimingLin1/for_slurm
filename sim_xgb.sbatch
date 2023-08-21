#!/bin/bash
#SBATCH --job-name simulation_xgb        # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition medium       # Slurm partition to use
#SBATCH --ntasks 1            # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 1-12:00          # Time limit in D-HH:MM
#SBATCH --mem-per-cpu=100000       # Memory limit for each tasks (in MB)
#SBATCH -o myscript_%j.out      # File to which STDOUT will be written
#SBATCH -e myscript_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL        # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=zmlin@uw.edu # Email to which notifications will be sent
 
module load R
R CMD BATCH CalibProp_simulation_ver.2.3.R
R CMD BATCH CalibProp_simulation_xgb.R