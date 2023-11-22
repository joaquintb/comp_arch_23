#!/bin/sh

. /etc/profile
module avail
module load gcc/12.1.0
echo "*** Running Unit Tests ***"
./build/utest/utest
echo "*** Running Functional Tests ***"
./build/ftest/ftest
