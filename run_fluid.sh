#!/bin/sh

. /etc/profile
module avail
module load gcc/12.1.0
echo "*** BASE ***"
perf stat -r 5 ./build/fluid/fluid 1000 ./in/large.fld final.fld
