#!/bin/bash
echo $1
echo $2
#
#$ -cwd
#$ -S /bin/bash
#$ -N vor_on_MD
#$ -q all.q
#$ -V
#$ -o vor_errors_and_outs/$JOB_NAME.$JOB_ID.out
#$ -e vor_errors_and_outs/$JOB_NAME.$JOB_ID.err
#
vor vor_params.in $1 $2
