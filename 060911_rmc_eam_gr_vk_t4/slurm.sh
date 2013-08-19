#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=univ                # default "univ" if not specified
#SBATCH --time=0-04:30:00               # run time in days-hh:mm:ss
##SBATCH --ntasks=32                    # require 32 CPUs (CPUs)
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4000              # RAM in MB (default 4GB, max 8GB)
#SBATCH --error=/home/maldonis/060911_rmc_eam_gr_vk_t4/job.%J.err
#SBATCH --output=/home/maldonis/060911_rmc_eam_gr_vk_t4/job.%J.out

#Now list your executable command (or a string of them). Example:
mpirun /home/maldonis/060911_rmc_eam_gr_vk_t4/rmc

#SBATCH --nodes=2                       # number of nodes requested
#SBATCH --ntasks-per-node=16            # default 16
##SBATCH --nodes=2                       # number of nodes requested
##SBATCH --ntasks-per-node=16            # default 16
##SBATCH --cpus-per task=1               # default 1
##SBATCH --mem=16384                     # total RAM in MB, max 64GB  per node

##SBATCH --export=ALL
