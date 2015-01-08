#!/bin/bash

# testing each executable using input/p1

./ad_pfm.x ../../input/p1 2 1 1 123456 | grep Final | awk -v p="ad_pfm p1 2 1 1" -v t=1 -f utest.awk
./ad_pfm.x ../../input/p1 2 2 1 123456 | grep Final | awk -v p="ad_pfm p1 2 2 1" -v t=1 -f utest.awk
./ad_pfm.x ../../input/p1 2 3 1 123456 | grep Final | awk -v p="ad_pfm p1 2 3 1" -v t=1 -f utest.awk

# testing each executable using input/p9

./ad_pfm.x ../../input/p9 2 1 1 123456 | grep Final | awk -v p="ad_pfm p9 2 1 1" -v t=75 -f utest.awk
./ad_pfm.x ../../input/p9 3 1 1 123456 | grep Final | awk -v p="ad_pfm p9 3 1 1" -v t=69 -f utest.awk
./ad_pfm.x ../../input/p9 2 2 1 123456 | grep Final | awk -v p="ad_pfm p9 2 2 1" -v t=53 -f utest.awk
./ad_pfm.x ../../input/p9 3 2 1 123456 | grep Final | awk -v p="ad_pfm p9 3 2 1" -v t=75 -f utest.awk
./ad_pfm.x ../../input/p9 2 3 2 123456 | grep Final | awk -v p="ad_pfm p9 2 3 2" -v t=32 -f utest.awk
./ad_pfm.x ../../input/p9 3 3 2 123456 | grep Final | awk -v p="ad_pfm p9 3 3 2" -v t=58 -f utest.awk

# EOF


