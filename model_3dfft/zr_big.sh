#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -m be
#$ -M voyles@engr.wisc.edu
#$ -e zr_big.err
#$ -i zr_big.input
#$ -o zr_big.out
#
/home/voyles/sim/model_fft3d/m3dfft
