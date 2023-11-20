. /etc/profile
module avail
module load gcc/12.1.0
echo "*** BASE ***"
perf stat 1 ./build/fluid/fluid 1000 "large.fld" "final.fld"
