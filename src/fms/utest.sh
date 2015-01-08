#!/bin/bash

# testing each executable using input/p1

./ad_fms.x ../../input/p1 2 123456 | grep Final | awk -v p="ad_fms p1 2" -v t=1 -f utest.awk

# testing each executable using input/p9

./ad_fms.x ../../input/p9 2 123456 | grep Final | awk -v p="ad_fms p9 2" -v t=85 -f utest.awk
./ad_fms.x ../../input/p9 3 123456 | grep Final | awk -v p="ad_fms p9 3" -v t=201 -f utest.awk

# EOF


