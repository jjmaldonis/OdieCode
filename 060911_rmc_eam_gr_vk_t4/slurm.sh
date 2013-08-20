#!/bin/sh

#SBATCH --job-name=rmc                  # job name
#SBATCH --partition=univ                # default "univ" if not specified
#SBATCH --error=/home/maldonis/060911_rmc_eam_gr_vk_t4/job.%J.err
#SBATCH --output=/home/maldonis/060911_rmc_eam_gr_vk_t4/job.%J.out

#SBATCH --time=0-04:30:00               # run time in days-hh:mm:ss

#SBATCH --ntasks=16                    # required number of CPUs
#SBATCH --nodes=1                       # number of nodes requested
#SBATCH --ntasks-per-node=16            # default 16
#SBATCH --ntasks-per-node=16            # default 16
##SBATCH --cpus-per task=1               # default 1
##SBATCH --mem=16384                     # total RAM in MB, max 64GB  per node
#SBATCH --mem-per-cpu=4000              # RAM in MB (default 4GB, max 8GB)

##SBATCH --export=ALL

echo "JobID = $SLURM_JOB_ID"
echo "Using $SLURM_NNODES nodes"
echo "Number of cores per node: $SLURM_TASKS_PER_NODE"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""

#Now list your executable command (or a string of them). Example:
mpirun rmc $SLURM_JOB_ID

