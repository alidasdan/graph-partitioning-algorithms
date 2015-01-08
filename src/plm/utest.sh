#!/bin/bash

# testing each executable using input/p1

./ad_plm.x ../../input/p1 2 1 1 123456 | grep Final | awk -v p="ad_plm p1 2 1 1" -v t=1 -f utest.awk
./ad_plm.x ../../input/p1 2 2 1 123456 | grep Final | awk -v p="ad_plm p1 2 2 1" -v t=1 -f utest.awk
./ad_plm.x ../../input/p1 2 1 2 123456 | grep Final | awk -v p="ad_plm p1 2 1 2" -v t=1 -f utest.awk

# testing each executable using input/p9

./ad_plm.x ../../input/p9 2 1 1 123456 | grep Final | awk -v p="ad_plm p9 2 1 1" -v t=26 -f utest.awk
./ad_plm.x ../../input/p9 3 1 1 123456 | grep Final | awk -v p="ad_plm p9 3 1 1" -v t=157 -f utest.awk
./ad_plm.x ../../input/p9 2 2 1 123456 | grep Final | awk -v p="ad_plm p9 2 2 1" -v t=35 -f utest.awk
./ad_plm.x ../../input/p9 3 2 1 123456 | grep Final | awk -v p="ad_plm p9 3 2 1" -v t=163 -f utest.awk
./ad_plm.x ../../input/p9 2 1 2 123456 | grep Final | awk -v p="ad_plm p9 2 1 2" -v t=27 -f utest.awk
./ad_plm.x ../../input/p9 3 1 2 123456 | grep Final | awk -v p="ad_plm p9 3 1 2" -v t=175 -f utest.awk

# EOF


