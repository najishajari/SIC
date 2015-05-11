#!/bin/sh
#
#$ -t 1-1024
# -o /agbs/cluster/naji/Multivariate IGCI/Oscillators/out/
#  Request 8G of RAM
#$ -l h_vmem=2G

JOBDIR=/agbs/cluster/naji/Linear\ Filters/
JOBOUT="$JOBDIR/out"
JOBDATA="$JOBDIR/data"

parallel_migci () {
	dataset="$1"
	/is/ei/naji/Enthought/Canopy_64bit/User/bin/python "/agbs/cluster/naji/Linear Filters/real_data/ec016/parallel_real_data_test.py" --t $dataset
}

parallel_migci $SGE_TASK_ID
