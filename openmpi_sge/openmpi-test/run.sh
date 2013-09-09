#!/bin/bash
#$ -q all.q
#$ -V
#$ -pe orte 4
#$ -cwd
#$ -e mpi-ring.err
#$ -o mpi-ring.out
#/share/apps/openmpi_intel_20130712/bin/mpirun -np $NSLOTS /home/jjmaldonis/OdieCode/openmpi_sge/openmpi-test/mpi-ring
#mpirun -np $NSLOTS amplxe-cl -r my_result -collect hotspots -- mpi-ring
mpirun -n $NSLOTS inspxe-cl -r my_result -collect mi1 -- /home/jjmaldonis/OdieCode/openmpi_sge/openmpi-test/mpi-ring

#amplxe-cl –r my_result -collect hotspots
#/share/apps/openmpi_intel_20130712/bin/mpiexec –np $NSLOTS amplxe-cl –r my_result -collect hotspots mpi-ring

#mpirun –n 4 inspxe-cl –r my_result -collect mi1 -- my_app [my_app_ options]
