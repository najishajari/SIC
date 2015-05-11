#!/bin/sh
#
#$ -t 1-15000
# -o /agbs/cluster/naji/Multivariate IGCI/Oscillators/out/
#$ -l h_vmem=0.5G

JOBDIR="/agbs/cluster/naji/Linear Filters/Echo/"
JOBOUT="$JOBDIR/out"
JOBDATA="$JOBDIR/data"

Echo () {
	dataset="$1"
	/is/ei/naji/Enthought/Canopy_64bit/User/bin/python "/agbs/cluster/naji/Linear Filters/Echo/man_echo.py" --t $dataset
}

Echo $SGE_TASK_ID
