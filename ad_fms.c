
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "ad_defs.h"
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_readinput.h"
#include "ad_partition.h"
#include "ad_print.h"
#include "ad_bucketio.h"
#include "ad_lib.h"
#include "ad_lib_fms.h"

/* FOR SANCHIS' VERSION OF MULTI-WAY PARTITIONING */
/* Also mentioned as the SN algorithm */
/* Direct multi-way partitioning.
   Locking is used.
   Cells are moved wrt their gains.
*/

/* definitions */
int nocells;           /* number of cells */
int noparts;           /* number of partitions */
int nonets;            /* number of nets */
int totsize;           /* total net weight of the partition */
int totcellsize;       /* total cell weight of the partition */
int cutsize;           /* cutsize of the partition */
int max_gain;          /* max gain of a cell */
int max_density;       /* max density of a cell */
int max_cweight;       /* max cell weight */
int max_nweight;       /* max net weight */
int bucketsize;        /* max size of a bucket array */
int msize;             /* index to mcells */

int main(int argc, char *argv[])
{
    if (argc < 3) {
        printf("\nUsage: %s InputFileName NoParts [Seed]\n", argv[0]);
        exit(1);
    }  /* if */

    char fname[STR_SIZE];
    sprintf(fname, "%s", argv[1]);

    noparts = atoi(argv[2]);                         

    long seed;
    if (argc > 3) {
        seed = (long) atoi(argv[3]);
    } else {
        seed = -1;
    }
    seed = randomize((long)  seed);
    printf("SEED = %ld fname = %s\n", seed, fname);

    read_graph_size(fname, &nocells, &nonets, noparts);

    /* alloc memory (statically if possible) */
    cells_t            cells[nocells];
    nets_t             nets[nonets];
    corn_t             cnets[2 * nonets];
    ind_t              pop[MAX_POP];             /* population */
    for (int i = 0; i < MAX_POP; ++i) {
        pop[i].chrom = (allele *) calloc(nocells, sizeof(allele));
        pop[i].parts = (parts_t *) calloc(noparts, sizeof(parts_t));
    }
    partb_t            partb[noparts][noparts - 1];  /* partition buckets */
    cells_info_t       cells_info[nocells];
    for (int i = 0; i < nocells; ++i) {
        cells_info[i].mgain = (int *) calloc(noparts, sizeof(int));
        cells_info[i].partb_ptr = (bnode_ptr_t *) calloc(noparts - 1, sizeof(bnode_ptr_t));
        cells_info[i].partb_gain_inx = (int *) calloc(noparts - 1, sizeof(int));
    }
    /* additional information for cells */
    selected_cell_t    scell[1];     /* selected cell type */
    mcells_t           mcells[nocells];  /* array of cells moved */
    parts_info_t       parts_info[noparts]; 
 
    read_graph(fname, nocells, nonets, noparts, &totsize, &totcellsize,
               &max_density, &max_cweight, &max_nweight,
               cells, nets, cnets);

    max_gain = max_density * max_nweight; 
    bucketsize = 2 * max_gain + 1;
#ifdef DEBUG
    printf("%d < %d\n", 2 * max_gain + 1, bucketsize);
#endif

    /* alloc memory (statically if possible) */
    for (int i = 0; i < noparts; ++i) {
        for (int j = 0; j < noparts - 1; ++j) {
            partb[i][j].bnode_ptr = (bnode_ptr_t *) calloc(bucketsize, sizeof(bnode_ptr_t));
        }
    }

    create_partition(nocells, noparts, totcellsize, 
                     cells, &pop[0]);

#ifdef DEBUG
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    init_buckets(noparts, bucketsize, partb);
    cutsize = find_cut_size(nonets, totsize, nets, &pop[0]);

#ifdef DEBUG
    printf("Totalsize = %d Initial cutsize = %d\n", totsize, cutsize);
#endif

    int gain_sum;
    int no_iter = 0;
    int glob_inx = 0;
    do {

        copy_partition(noparts, parts_info, &pop[0]);

        compute_gains(nocells, noparts, pop[0].chrom, 
                      cells, nets, cnets, cells_info);

        create_buckets(nocells, noparts, max_gain, pop[0].chrom, 
                       partb, cells_info);

        msize = 0;

        int nlocked = 0;
        do {

            int move_possible = select_cell(noparts, scell, parts_info, cells, 
                                            partb, cells_info);

            delete_partb_nodes_of_cell(noparts, scell[0].mov_cell_no, 
                                       scell[0].from_part, partb, cells_info);

            /* lock cell */
            cells_info[scell[0].mov_cell_no].locked = True;
            if (move_possible == True) {
                move_cell(mcells, msize, scell);  
                msize++;
                update_gains(noparts, max_gain, scell, pop[0].chrom,
                             cells, nets, cnets,
                             partb, cells_info);
            }   /* if */
            nlocked++;

        } while (nlocked < nocells); 

        int max_mcells_inx;
        gain_sum = find_move_set(mcells, msize, &max_mcells_inx);

#ifdef DEBUG
        printf("gain_sum=%d max_mcells_inx=%d msize = %d\n",
               gain_sum, max_mcells_inx, msize);
#endif

        if (gain_sum > 0) {
            int cut_gain = move_cells(False, mcells, max_mcells_inx, cutsize,
                                      &glob_inx, &pop[0], cells);
            cutsize -= cut_gain;
        }   /* if */
        no_iter++;

#ifdef DEBUG
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", no_iter,
               cutsize, find_cut_size(nonets, totsize, nets, &pop[0]));
#endif

    } while ((gain_sum > 0) && (cutsize > 0) && (no_iter < NO_ITERATIONS));

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", no_iter,
           cutsize, find_cut_size(nonets, totsize, nets, &pop[0]));

    free_nodes(noparts, bucketsize, partb);

#ifdef DEBUG
    printf("Final : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("FF %d %d %d %d\n", i, pop[0].parts[i].pmin_size,
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    /* free memory */
    for (int i = 0; i < MAX_POP; ++i) {
        cfree(pop[i].chrom);
        cfree(pop[i].parts);
    }
    for (int i = 0; i < nocells; ++i) {
        cfree(cells_info[i].mgain);
        cfree(cells_info[i].partb_ptr);
        cfree(cells_info[i].partb_gain_inx);
    }
    for (int i = 0; i < noparts; ++i) {
        for (int j = 0; j < noparts - 1; ++j) {
            cfree(partb[i][j].bnode_ptr);
        }
    }

    return (0);
}   /* main-fms */

/* EOF */
