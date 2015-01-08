#ifndef AD_PARTITION_INCLUDED
#define AD_PARTITION_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* create initial partition */
int create_partition(int nocells, 
                     int noparts, 
                     int totsize, 
                     cells_t cells[], 
                     ind_t *ind);

/* copy partition properties to parts_info for temporary use */
void copy_partition(int noparts, 
                    parts_info_t parts_info[], 
                    ind_t *ind);

/* read a partition prepared beforehand */
int read_partition(FILE *fp, 
                   char *filename, 
                   int noparts, 
                   ind_t *ind);

/* write a partition */
void write_partition(FILE *fp, 
                     char *filename, 
                     int nocells, 
                     int noparts, 
                     ind_t *ind);

#endif
