
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include <stdlib.h>
#include "ad_defs.h"
#include "ad_lib.h"

/* find the neighbor cell of cell_no, on the net net_no */
void find_other_cell(int net_no, 
                     int cell_no,
                     int *other_cell, 
                     int *other_part_no, 
                     int *net_weight,
                     allele chrom[], 
                     nets_t nets[])
{
    /* first find the index of other_cell */
    int cell_inx = 0;
    if (nets[net_no].ncells[0] == cell_no) {
        cell_inx = 1;
    }
    *other_cell = nets[net_no].ncells[cell_inx];
    *other_part_no = chrom[*other_cell];
    *net_weight = nets[net_no].nweight;
}   /* find_other_cell */
 
/* initialize all bucket indices and pointers */
void init_buckets(int noparts, 
                  int bucketsize,
                  partb_t partb[][noparts - 1])
{
    /* init partition bucket indices */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {

            partb[i][j].max_inx = partb[i][j].min_inx = -1;
            partb[i][j].nobuckets = 0;

            /* init partb bucket pointers */
            for (int k = 0; k < bucketsize; k++) {
                partb[i][j].bnode_ptr[k] = NULL;
            }

        }   /* for j */
    }   /* for i */
}   /* init_buckets */
 
/* map part no such that home_part is excluded */
int map_part_no(int dest_part, int home_part)
{
    if (dest_part < home_part) {
        return (dest_part);
    } else if (dest_part > home_part) {
        return (dest_part - 1);
    } else {
        printf("Error: Unexpected inputs\n");
        exit(1);
    }
}   /* map_part_no */
 
/* compute move gain from home_part to dest_part */
int calculate_gain(int cell_no, 
                   int home_part, 
                   int dest_part,
                   cells_info_t cells_info[])
{
    int mov_gain = cells_info[cell_no].mgain[dest_part] - 
        cells_info[cell_no].mgain[home_part];
    return mov_gain;
}   /* calculate_gain */

/* compute gains of all cells and place them into cells_info */
void compute_gains(int nocells, 
                   int noparts,
                   allele tchrom[],
                   cells_t cells[],
                   nets_t nets[],
                   corn_t cnets[],
                   cells_info_t cells_info[])
{
    for (int cell_no = 0; cell_no < nocells; cell_no++) {

        /* initialize cells_info */
        cells_info[cell_no].locked = False;
        cells_info[cell_no].mcount = 0;

        /* find info about cell_no */
        int net_ptr = cells[cell_no].netlist;   /* ptr to nets of cell_no */
        cells[cell_no].cno_inets = 0;       /* # of internal nets of cell_no */

        /* init external & internal costs of moving cell_no to each partition */
        for (int j = 0; j < noparts; j++) {
            cells_info[cell_no].mgain[j] = 0;
        }

        /* for each cell connected to cell_no do */
        for (int j = 0; j < cells[cell_no].cno_nets; j++) {

            /* find the neighbor cell */
            int net_no = cnets[net_ptr].corn_no;

            int other_cell, other_part_no, net_weight;
            find_other_cell(net_no, cell_no,
                            &other_cell, &other_part_no, &net_weight,
                            tchrom, nets);
            cells_info[cell_no].mgain[other_part_no] += net_weight;
            net_ptr++;  

        }   /* for j */

    }   /* for cell_no */
}   /* compute_gains */

/* free all allocated nodes */
void free_nodes(int noparts, 
                int bucketsize,
                partb_t partb[][noparts - 1])
{
    /* delete nodes connected to partb */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {

            for (int k = 0; k < bucketsize; k++) {
                bnode_ptr_t next = partb[i][j].bnode_ptr[k];
                partb[i][j].bnode_ptr[k] = NULL;
                while (next != NULL) {
                    bnode_ptr_t prev = next;
                    next = next->rptr;
                    free(prev);
                }   /* while */
            }   /* for k */

            partb[i][j].max_inx = -1;
            partb[i][j].min_inx = -1;
            partb[i][j].nobuckets = 0;
        }   /* for j */
    }   /* for i */
}   /* free_nodes */

/* count number of bucket nodes */
void number_nodes(int noparts, 
                  int bucketsize, 
                  int *npartb,
                  partb_t partb[][noparts - 1])
{
    *npartb = 0;

    /* count nodes connected to partb */
    for (int i = 0; i < noparts; i++) {
        for (int j = 0; j < (noparts - 1); j++) {
            for (int k = 0; k < bucketsize; k++) {
                bnode_ptr_t next = partb[i][j].bnode_ptr[k];
                while (next != NULL) {
                    next = next->rptr;
                    (*npartb)++;
                }   /* while */
            }   /* for k */
        }   /* for j */
    }   /* for i */
}   /* number_nodes */

/* find set of cells to be actually moved */
int find_move_set(mcells_t mcells[],
                  int msize,
                  int *max_mcells_inx)
{
    int max_gain_sum = 0;
    *max_mcells_inx = -1;
    int gain_sum = 0;
    for (int i = 0; i < msize; i++) {
        gain_sum += mcells[i].mgain;
        if (gain_sum > max_gain_sum) {
            *max_mcells_inx = i;
            max_gain_sum = gain_sum;
        }   /* if */
    }   /* for i */

    return max_gain_sum;
}   /* find_move_set */

/* move cells actually */
int move_cells(int wflag,
               mcells_t mcells[],
               int max_mcells_inx, 
               int cutsize, 
               int *glob_inx,
               ind_t *ind,
               cells_t cells[])
{
    /* used only under wflag == True */
    int tcutsize, fcutsize;
    tcutsize = fcutsize = cutsize;

    int cut_gain = 0;

    for (int i = 0; i <= max_mcells_inx; i++) {

        if (wflag == True) {
            tcutsize -= mcells[i].mgain;
            if (tcutsize < fcutsize) fcutsize = tcutsize;
            printf ("%d %d %d\n", *glob_inx, tcutsize, fcutsize);
            (*glob_inx)++;
        }    /* if wflag */

        ind->chrom[mcells[i].cell_no] = mcells[i].to;
        cut_gain += mcells[i].mgain;

        /* update partition size limits */
        ind->parts[mcells[i].from].pmax_cells--;
        ind->parts[mcells[i].from].pcurr_size -= cells[mcells[i].cell_no].cweight;
        ind->parts[mcells[i].to].pmax_cells++;
        ind->parts[mcells[i].to].pcurr_size += cells[mcells[i].cell_no].cweight;

    }   /* for i */

    return cut_gain;
}   /* find_move_set */

/* finds cut size of a given partition - used for control */
int find_cut_size(int nonets, 
                  int totsize,
                  nets_t nets[],
                  ind_t *ind)
{
    ind->incost = 0;
    for (int i = 0; i < nonets; i++) {
        if (ind->chrom[nets[i].ncells[0]] ==
            ind->chrom[nets[i].ncells[1]])
            ind->incost += nets[i].nweight;
    }   /* for i */

    return (totsize - ind->incost);
}   /* find_cut_size */

/* EOF */
