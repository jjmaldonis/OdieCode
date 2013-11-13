#!/bin/bash
#$ -N md_Emin
#$ -cwd
#$ -S /bin/bash
#$ -q all.q
#$ -pe orte 1
#$ -e energy_min_RT.err
#$ -o energy_min_RT.out
#$ -i energy_min_RT.in
#
/opt/openmpi/bin/mpirun -np $NSLOTS /home/jjmaldonis/bin/lmp_odie
