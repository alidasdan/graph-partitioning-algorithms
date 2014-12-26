#!/bin/bash

# testing each executable using input/p1

./ad_fms.x ./input/p1 2 123456 | grep Final | awk -v p="ad_fms p1 2" -v t=1 -f utest.awk
./ad_plm.x ./input/p1 2 1 1 123456 | grep Final | awk -v p="ad_plm p1 2 1 1" -v t=1 -f utest.awk
./ad_plm.x ./input/p1 2 2 1 123456 | grep Final | awk -v p="ad_plm p1 2 2 1" -v t=1 -f utest.awk
./ad_plm.x ./input/p1 2 1 2 123456 | grep Final | awk -v p="ad_plm p1 2 1 2" -v t=1 -f utest.awk
./ad_pfm.x ./input/p1 2 1 123456 | grep Final | awk -v p="ad_pfm p1 2 1" -v t=1 -f utest.awk
./ad_pfm.x ./input/p1 2 2 123456 | grep Final | awk -v p="ad_pfm p1 2 2" -v t=1 -f utest.awk
./ad_pfm.x ./input/p1 2 3 123456 | grep Final | awk -v p="ad_pfm p1 2 3" -v t=1 -f utest.awk

# testing each executable using input/p9

./ad_fms.x ./input/p9 2 123456 | grep Final | awk -v p="ad_fms p9 2" -v t=85 -f utest.awk
./ad_fms.x ./input/p9 3 123456 | grep Final | awk -v p="ad_fms p9 3" -v t=201 -f utest.awk
./ad_plm.x ./input/p9 2 1 1 123456 | grep Final | awk -v p="ad_plm p9 2 1 1" -v t=26 -f utest.awk
./ad_plm.x ./input/p9 3 1 1 123456 | grep Final | awk -v p="ad_plm p9 3 1 1" -v t=157 -f utest.awk
./ad_plm.x ./input/p9 2 2 1 123456 | grep Final | awk -v p="ad_plm p9 2 2 1" -v t=35 -f utest.awk
./ad_plm.x ./input/p9 3 2 1 123456 | grep Final | awk -v p="ad_plm p9 3 2 1" -v t=163 -f utest.awk
./ad_plm.x ./input/p9 2 1 2 123456 | grep Final | awk -v p="ad_plm p9 2 1 2" -v t=27 -f utest.awk
./ad_plm.x ./input/p9 3 1 2 123456 | grep Final | awk -v p="ad_plm p9 3 1 2" -v t=175 -f utest.awk
./ad_pfm.x ./input/p9 2 1 123456 | grep Final | awk -v p="ad_pfm p9 2 1" -v t=75 -f utest.awk
./ad_pfm.x ./input/p9 3 1 123456 | grep Final | awk -v p="ad_pfm p9 3 1" -v t=69 -f utest.awk
./ad_pfm.x ./input/p9 2 2 123456 | grep Final | awk -v p="ad_pfm p9 2 2" -v t=53 -f utest.awk
./ad_pfm.x ./input/p9 3 2 123456 | grep Final | awk -v p="ad_pfm p9 3 2" -v t=75 -f utest.awk
./ad_pfm.x ./input/p9 2 3 123456 | grep Final | awk -v p="ad_pfm p9 2 3" -v t=50 -f utest.awk
./ad_pfm.x ./input/p9 3 3 123456 | grep Final | awk -v p="ad_pfm p9 3 3" -v t=65 -f utest.awk

# EOF


