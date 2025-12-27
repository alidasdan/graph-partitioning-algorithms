#ifndef AD_LIB_PFM_INCLUDED
#define AD_LIB_PFM_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* map a given mov_gain into index to a bucket array */
int map_gain(int bucketsize,
             int mov_gain, 
             int mov_count, 
             int max_gain,
             eval_t eval[]);
 
/* fill all bucket arrays */
void create_buckets(int bucketsize,
                    int nocells,
                    int noparts,
                    int max_gain,
                    eval_t eval[],
                    allele chrom[],
                    partb_t **partb,
                    cells_info_t cells_info[]);

/* select a cell to move */
void select_cell(int noparts,
                 selected_cell_t scell[],
                 parts_info_t *parts_info,
                 cells_t cells[],
                 partb_t **partb,
                 cells_info_t cells_info[]);

/* move selected cell, and save the move in a file */
void move_cell(mcells_t mcells[],
               int msize,
               selected_cell_t scell[],
               allele tchrom[],
               cells_info_t cells_info[]);

/* update gains after a move */
void update_gains(int bucketsize,
                  int noparts,
                  int max_gain,
                  eval_t eval[],
                  selected_cell_t scell[],
                  allele tchrom[],
                  cells_t cells[],
                  nets_t nets[],
                  corn_t cnets[],
                  partb_t **partb,
                  cells_info_t cells_info[]);

/* fill eval array */
void fill_eval(int max_gain,
               float K,
               eval_t eval[]);

void create_partb_nodes_of_cell(int bucketsize,
                                int noparts,
                                int max_gain,
                                int cell_no,
                                int part_no,
                                eval_t eval[],
                                partb_t **partb,
                                cells_info_t cells_info[]);

/* calculate K and scale factor */
void calculate_scale(int nocells, 
                     int noparts, 
                     int max_gain,
                     float *K);

#endif
