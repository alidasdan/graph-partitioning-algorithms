#ifndef AD_LIB_INCLUDED
#define AD_LIB_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* find the neighbor cell of cell_no, on the net net_no */
void find_other_cell(int net_no, 
                     int cell_no,
                     int *other_cell, 
                     int *other_part_no, 
                     int *net_weight,
                     allele chrom[], 
                     nets_t nets[]);
 
/* initialize all bucket indices and pointers */
void init_buckets(int noparts,
                  int bucketsize,
                  partb_t **partb);
 
/* map part no such that home_part is excluded */
int map_part_no(int dest_part, int home_part);
 
/* compute move gain from home_part to dest_part */
int calculate_gain(int cell_no, 
                   int home_part, 
                   int dest_part, 
                   cells_info_t cells_info[]);

/* compute gains of all cells and place them into cells_info */
void compute_gains(int nocells, 
                   int noparts, 
                   allele tchrom[],
                   cells_t cells[], 
                   nets_t nets[], 
                   corn_t cnets[], 
                   cells_info_t cells_info[]);

/* free all allocated nodes */
void free_nodes(int noparts,
                int bucketsize,
                partb_t **partb);

/* count number of bucket nodes */
void number_nodes(int noparts,
                  int bucketsize,
                  int *npartb,
                  partb_t **partb);

/* find set of cells to be actually moved */
int find_move_set(mcells_t mcells[], 
                  int msize, 
                  int *max_mcells_inx);

/* move cells actually */
int move_cells(int wflag, 
               mcells_t mcells[], 
               int max_mcells_inx, 
               int cutsize, 
               int *glob_inx, 
               ind_t *ind, 
               cells_t cells[]);

/* finds cut size of a given partition - used for control */
int find_cut_size(int nonets, 
                  int totsize, 
                  nets_t nets[], 
                  ind_t *ind);

#endif
