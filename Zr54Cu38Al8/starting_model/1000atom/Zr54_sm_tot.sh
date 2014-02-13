##!/bin/bash
#
#$ -cwd
#$ -q all.q
#$ -S /bin/bash
#$ -m be
#$ -M voyles@engr.wisc.edu
#$ -pe orte 4
#$ -e Zr54_sm_tot.err
#$ -o Zr54_sm_tot.out
#$ -i Zr54_sm_tot.in
#$ -V
#
/share/apps/openmpi_intel_20130712/bin/mpiexec -n $NSLOTS /home/voyles/bin/lmp_odie

