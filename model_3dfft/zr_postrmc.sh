#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -m be
#$ -M voyles@engr.wisc.edu
#$ -e zr_postrmc.err
#$ -i zr_postrmc.input
#$ -o zr_postrmc.out
#
/home/voyles/sim/model_fft3d/m3dfft
