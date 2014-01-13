#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -m be
#$ -M voyles@engr.wisc.edu
#$ -e alsm_fcc.err
#$ -i alsm_fcc.input
#$ -o alsm_fcc.out
#
/home/voyles/sim/model_fft3d/m3dfft
