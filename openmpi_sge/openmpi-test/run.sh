#!/bin/bash
#$ -q testqueue.q
#$ -V
#$ -pe orte 128
#$ -cwd
#$ -e trun.err
#$ -o trun.out
/share/apps/openmpi_intel_20130712/bin/mpiexec -np $NSLOTS /home/jjmaldonis/OdieCode/openmpi_sge/openmpi-test/mpi-ring
