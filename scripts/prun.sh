#!/bin/sh

echo "qsub -cwd -V -pe MPI $1 -q $2.q ./runme.sh"
qsub -cwd -V -pe MPI "$1" -q "$2".q ./runme.sh