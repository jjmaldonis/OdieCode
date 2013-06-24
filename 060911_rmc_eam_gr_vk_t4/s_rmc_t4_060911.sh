##!/bin/bash
#
#$ -N rmc_060911_rmc_t4
#$ -cwd
#$ -S /bin/bash
#$ -m be
# -M jhwang3@wisc.edu
#$ -e standard_error.err
# -i standard_input.input
#$ -o rmc_060911_t4.out
# 
#$ -V
#$ -l h_rt=48:00:00
#$ -q long
#$ -A TG-DMR100061

ibrun -n $NSLOTS -o 0 $HOME/060911_rmc_eam_gr_vk_t4/rmc_test 

