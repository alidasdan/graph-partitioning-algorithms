
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <malloc.h>
#include "ad_defs.h"
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_readinput.h"
#include "ad_partition.h"
#include "ad_print.h"
#include "ad_bucketio.h"
#include "ad_lib.h"
#include "ad_lib_pfm.h"

/* PARTITIONING BY FREE MOVES */

int main(int argc, char *argv[])
{
    /* definitions */
    int nocells;           /* number of cells */
    int nonets;            /* number of nets */
    int nopins;            /* number of pins */
    int noparts;           /* number of partitions */
    int totsize;           /* total net weight of the partition */
    int totcellsize;       /* total cell weight of the partition */
    int cutsize;           /* cutsize of the partition */
    int max_gain;          /* max gain of a cell */
    int max_density;       /* max density of a cell */
    int max_cweight;       /* max cell weight */
    int max_nweight;       /* max net weight */

    if (argc < 5) {
        printf("\nUsage: %s InputFileName NoParts Version BucketSizeFactor [Seed]\n", argv[0]);
        exit(1);
    }  /* if */

    char fname[STR_SIZE];
    sprintf(fname, "%s", argv[1]);

    noparts = atoi(argv[2]);
    int version = atoi(argv[3]);

    int bucketsize_factor = atoi(argv[4]);

    long seed;
    if (argc > 5) {
        seed = (long) atoi(argv[5]);
    } else {
        seed = (long) -1;
    }
    seed = randomize((long)  seed);
    printf("SEED = %ld fname = %s\n", seed, fname);

    read_graph_size(fname, &nocells, &nonets);
    nopins = 2 * nonets;

    /* determine max_noiter based on pfm version */
    /* pfm1: size=max_cells; 
       pfm2: size=max_cells * max_parts; 
       pfm3: size=max_cells * max_parts^2 */
    int max_noiter = nocells;
    switch (version) {
    case 1 : break;
    case 2 : max_noiter *= noparts; break;
    case 3 : max_noiter *= noparts * noparts; break;
    default : break;
    }
    cells_t *cells = (cells_t *) calloc(nocells, sizeof(cells_t));
    assert(cells != NULL);
    cells_info_t *cells_info = (cells_info_t *) calloc(nocells, sizeof(cells_info_t));
    assert(cells_info != NULL);
    for (int i = 0; i < nocells; i++) {
        cells_info[i].mgain = (int *) calloc(noparts, sizeof(int));
        cells_info[i].partb_ptr = (bnode_ptr_t *) calloc(noparts - 1, sizeof(bnode_ptr_t));
        cells_info[i].partb_gain_inx = (int *) calloc(noparts - 1, sizeof(int));
    }

    nets_t *nets = (nets_t *) calloc(nonets, sizeof(nets_t));
    assert(nets != NULL);

    /* cells of nets */
    corn_t *cnets = (corn_t *) calloc(nopins, sizeof(corn_t));
    assert(cnets != NULL);

    /* partition buckets */
    partb_t partb[noparts][noparts - 1];  
    parts_info_t parts_info[noparts]; 

    /* population (w/ one individual!) */
    ind_t pop[MAX_POP];             
    for (int i = 0; i < MAX_POP; i++) {
        pop[i].chrom = (allele *) calloc(nocells, sizeof(allele));
        pop[i].parts = (parts_t *) calloc(noparts, sizeof(parts_t));
    }

    /* selected cell */
    selected_cell_t scell[1];     
    selected_cell_t prev_scell[1];

    /* moved cells */
    mcells_t *mcells = (mcells_t *) calloc(2 * max_noiter, sizeof(mcells_t));
    assert(mcells != NULL);

    /* temp chrom */
    allele *tchrom = (allele *) calloc(nocells, sizeof(allele));
    assert(tchrom != NULL);

    read_graph(fname, nocells, nonets, noparts, &totsize, &totcellsize,
               &max_density, &max_cweight, &max_nweight,
               cells, nets, cnets);

    /* bucketsize has impact on cutsize and runtime */
    max_gain = max_density * max_nweight;
    int bucketsize = 2 * max_gain + 1;
    if (bucketsize_factor > 0) {
        bucketsize *= bucketsize_factor;
    }

    /* cache to speed up math heavy function evals */
    eval_t *eval = (eval_t *) calloc(2 * max_gain + 1, sizeof(eval_t));
    float K;
    calculate_scale(nocells, noparts, max_gain, &K);
    fill_eval(max_gain, K, eval); 

    create_partition(nocells, noparts, totcellsize, cells, &pop[0]);

#ifdef DEBUG
    printf("Initial : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("II %d %d %d %d\n", i, pop[0].parts[i].pmin_size, 
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif

    /* alloc memory (statically if possible) */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            partb[i][j].bnode_ptr = (bnode_ptr_t *) calloc(bucketsize, sizeof(bnode_ptr_t));
        }
    }

    init_buckets(noparts, bucketsize, partb);
    cutsize = find_cut_size(nonets, totsize, nets, &pop[0]);

#ifdef DEBUG
    printf("BB bucketsize=%d\n", bucketsize);
    printf("Totalsize = %d Initial cutsize = %d\n", totsize, cutsize); 
#endif

    int gain_sum;
    int glob_inx = 0;
    int pass_no = 0;
    do {

        copy_partition(noparts, parts_info, &pop[0]);

        for (int i = 0; i < nocells; i++) {
            tchrom[i] = pop[0].chrom[i];
        }

        compute_gains(nocells, noparts, tchrom,
                      cells, nets, cnets, cells_info);

        create_buckets(bucketsize, nocells, noparts, max_gain, eval, pop[0].chrom,
                       partb, cells_info);

        prev_scell[0].mov_cell_no = -1;

        int msize = 0; /* index to mcells */
        int mov_count = 0; /* count of total moves */
        while (mov_count < max_noiter) {

            select_cell(noparts, scell, parts_info, cells, partb, cells_info);

            move_cell(mcells, msize, scell, tchrom, cells_info);  

            msize++;

            /* insert previous cell */
            create_partb_nodes_of_cell(bucketsize, noparts, max_gain,
                                       prev_scell[0].mov_cell_no,
                                       prev_scell[0].to_part,
                                       eval, partb, cells_info);

            update_gains(bucketsize, noparts, max_gain, eval, scell, tchrom, cells,
                         nets, cnets, partb, cells_info);

            /* currently moved cell becomes previous cell */
            bcopy((void *) scell, (void *) prev_scell, sizeof(scell[0]));

            /* delete current cell */
            delete_partb_nodes_of_cell(noparts, scell[0].mov_cell_no,
                                       scell[0].from_part, partb, cells_info);

            mov_count++;
        }   /* while */

        int max_mcells_inx;
        gain_sum = find_move_set(mcells, msize, &max_mcells_inx);

#ifdef DEBUG
        printf("gain_sum=%d max_mcells_inx=%d msize = %d\n", 
               gain_sum, max_mcells_inx, msize);
#endif

        if (gain_sum > 0) {
            int cut_gain = move_cells(False, mcells, max_mcells_inx, 
                                      cutsize, &glob_inx, &pop[0], cells);
            cutsize -= cut_gain;
        }   /* if */

        pass_no++;

#ifdef DEBUG
        printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", pass_no, 
               cutsize, find_cut_size(nonets, totsize, nets, &pop[0]));
#endif

        free_nodes(noparts, bucketsize, partb);

    }  while ((gain_sum > 0) && (cutsize > 0)); /* while pass_no */
    /* or (pass_no < noparts) */

    printf("pass_no = %d Final cutsize = %d Check cutsize = %d\n", pass_no, 
           cutsize, find_cut_size(nonets, totsize, nets, &pop[0]));

#ifdef DEBUG
    printf("Final : Part_no min_size curr_size max_size\n");
    for (int i = 0; i < noparts; i++) {
        printf("FF %d %d %d %d\n", i, pop[0].parts[i].pmin_size, 
               pop[0].parts[i].pcurr_size, pop[0].parts[i].pmax_size);
    }
#endif
 
    /* free memory for all data structures */
    cfree(cells);
    for (int i = 0; i < nocells; i++) {
        cfree(cells_info[i].mgain);
        cfree(cells_info[i].partb_ptr);
        cfree(cells_info[i].partb_gain_inx);
    }
    cfree(cells_info);

    cfree(nets);

    cfree(cnets);

    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < noparts - 1; ++j) {
            cfree(partb[i][j].bnode_ptr);
        }
    }

    for (int i = 0; i < MAX_POP; i++) {
        cfree(pop[i].chrom);
        cfree(pop[i].parts);
    }

    cfree(mcells);

    cfree(tchrom);

    cfree(eval);

    return (0);
}   /* main-pfm */

/* EOF */
