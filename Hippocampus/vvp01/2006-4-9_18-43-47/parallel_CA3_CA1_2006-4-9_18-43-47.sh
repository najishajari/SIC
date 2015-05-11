#!/bin/sh
#
#$ -t 1-1025
# -o /agbs/cluster/naji/Multivariate IGCI/Oscillators/out/
#  Request 8G of RAM
#$ -l h_vmem=8G

JOBDIR=/agbs/cluster/naji/Linear\ Filters/
JOBOUT="$JOBDIR/out"
JOBDATA="$JOBDIR/data"

parallel_migci () {
	dataset="$1"
	/is/ei/naji/Enthought/Canopy_64bit/User/bin/python "/agbs/cluster/naji/Linear Filters/real_data/vvp01/2006-4-9_18-43-47/parallel_CA3_CA1_2006-4-9_18-43-47.py" --t $dataset
}

parallel_migci $SGE_TASK_ID
