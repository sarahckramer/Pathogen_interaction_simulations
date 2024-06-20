#!/bin/bash -l

# Standard output and error:
#SBATCH -o temp/job-%A_%a.out # Standard output, %A = job ID, %a = job array index
#SBATCH -e temp/job-%A_%a.err # Standard error, %A = job ID, %a = job array index

# Job Name:
#SBATCH -J tm

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=pirikahu@mpiib-berlin.mpg.de

# --- resource specification (which resources for how long) ---
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100000 # memory in MB required by the job

#SBATCH --time=24:00:00 # run time in h:m:s, up to 24h possible
 
# --- start from a clean state and load necessary environment modules ---
module purge
module load R/4.2

# --- run your executable via srun ---
R --no-save --no-restore <seitr_x_seitr_v2.R >Rout/R-tm-${SLURM_ARRAY_TASK_ID}.Rout
