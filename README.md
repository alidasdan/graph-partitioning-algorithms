graph-partitioning-algorithms
=============================

This directory contains multi-way partitioning algorithms: FMS
(Fiduccia-Mattheyses-Sanchis), PLM (Partitioning by Locked Moves), PFM
(Partitioning by Free Moves) as detailed in [DaAy97]. 

A short history: The basis for these algorithms go back to the
Kernighan-Lin (KL) algorithm for graph partitioning. The KL algorithm
produces very good partitions but it is slow. The Fiduccia-Mattheyses
(FM) algorithm is not only a faster version of the KL algorithm but it
also generalizes the KL algorithm to run on hypergraphs. Sanchis
generalized the FM algorithm from 2-way partitioning to multi-way
partitioning. Seeing the limitations of the "locking" mechanism, I
devised a way to relax this mechanism, resulting in the Partitioning
by Locked Moves (PLM) algorithm and the Partitioning by Free Moves
(PFM) algorithm. The PLM algorithm still uses locking but it goes
through multiple phases of locking and unlocking within a pass over
the input graph. The PFM algorithm does not use locking at all; it
uses a way to penalize the moves to march towards a local minimum. For
more details, please refer to [DaAy97].

I originally developed this package (written in C) during my MSc study
(around 1991-1993). Before putting this package on github, I converted
my code to ANSI C (c99). 

This package is available on an "as is" basis. I do not say or imply
that it will be useful for whatever you want to do with it. It may
also contain bugs, and I assume no responsibility for any potential
problems associated with its use. You can use this package free of
charge in academic research and teaching. For any commercial use,
contact Ali Dasdan at ali_dasdan@yahoo.com. See the COPYRIGHT section
below.

## HOW TO COMPILE AND BUILD

Type 'make' (or 'gmake') to build all executables. The executables are
all have .x extension: ad_fms.x, ad_plm.x, ad_pfm.x.

With no targets following the make command, the following executables
will be generated:
- 'ad_fms.x'
- 'ad_plm.x'
- 'ad_pfm.x;

# HOW TO RUN

Type the name of one of the executables in your command line to get
the usage information. At minimum, each executable requires the input
graph and the number of parts to partition the graph. PlM and PFM
require additional parameters to create different versions of them,
which trade off runtime for partition quality.

## HOW TO TEST

Type 'make test' to test each executable on the input graphs under the
'input' directory. The result will be a 'pass' or a 'fail'.

## INPUT FILE FORMAT

The input file format is explained below.

```
> cat input/p1
6
7
1 1 0 1
1 1 1 2
1 1 2 0
1 1 3 4
1 1 4 5
1 1 5 3
1 1 0 3
1
1
1
1
1
1

```

The first two lines give the number of vertices (or cells) and the
number of edges (or nets), respectively. Thus, 'p1' has 6 vertices and
7 edges. 

The following 7 lines describe the edges, one edge per line. The first
number is the edge weight, the second number to be ignored, the third
number is the source vertex and the fourth number is the target
vertex. For example, the first edge from vertex 0 to vertex 1 has a
weight of 1.

The last 6 lines describe the vertex weights.

## REFERENCE

Please cite this reference if you use my programs in your research
work.

```
@article{DaAy97,
 author = {Ali Dasdan and C. Aykanat},
 title = {Two Novel Circuit Partitioning Algorithms Using Relaxed Locking},
 journal = {IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD)},
 volume = {16},
 number = {2},
 year = {1997},
 pages = {169-178},
 }
```

## COPYRIGHT

COPYRIGHT C 1991 - Ali Dasdan

Permission to use for non-commercial purposes is granted provided that
proper acknowledgments are given. For a commercial licence, contact
Ali Dasdan at ali_dasdan@yahoo.com.

This software is provided on an "as is" basis, without warranties or
conditions of any kind, either express or implied including, without
limitation, any warranties or conditions of title, non-infringement,
merchantability or fitness for a particular purpose.

## END OF FILE

