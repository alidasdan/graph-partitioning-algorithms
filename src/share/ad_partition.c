
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdlib.h>
#include <string.h>
#include <malloc/malloc.h>
#include <errno.h>
#include "ad_random.h"
#include "ad_fileio.h"
#include "ad_defs.h"
#include "ad_partition.h"

/* create initial partition */
int create_partition(int nocells,
                     int noparts,
                     int totsize,
                     cells_t cells[],
                     ind_t *ind)
{
    /* init current size */
    for (int i = 0; i < noparts; i++) {
        ind->parts[i].pratio = 1.0 / noparts;
        ind->parts[i].pcurr_size = 0;   
        ind->parts[i].pmax_cells = 0;
    }   /* for i */

    /* allocate temporary memory */
    int *tparts = (int *) calloc(noparts, sizeof(int));
    if (tparts == NULL) {
        printf("Error: Unable to allocate memory for tparts.\n");
        exit(1);
    }   /* if */

    /* insert cells */
    for (int i = 0; i < nocells; i++) {

        /* find minimum */
        int min_inx = 0;
        int min_size = ind->parts[0].pcurr_size;
        for (int j = 1; j < noparts; j++) {
            if (min_size > ind->parts[j].pcurr_size) {
                min_inx = j;
                min_size = ind->parts[j].pcurr_size;
            }   /* if */
        }   /* for j */

        /* find a new minimum among the same minimums */
        int tcount = -1; /* tparts is used to randomize partitioning */
        for (int j = 0; j < noparts; j++) {
            if (ind->parts[j].pcurr_size == min_size) {
                tcount++; 
                tparts[tcount] = j;
            }   /* if */
        }   /* for j */
        min_inx = tparts[irandom(0, tcount)];

        /* assign cell i to part[min_inx] */
        ind->chrom[i] = min_inx;
        ind->parts[min_inx].pcurr_size += cells[i].cweight; 
        ind->parts[min_inx].pmax_cells++; 
    }   /* for i */

    float off_ratio = 0.1;
    int max_size = -1;
    for (int i = 0; i < noparts; i++) {
        float part_size = ((float) totsize) * ind->parts[i].pratio;
        ind->parts[i].pmax_size = (int) (part_size * (1.0 + off_ratio) + 1.0);
        if (ind->parts[i].pmax_size > max_size) { 
            max_size = ind->parts[i].pmax_size;
        }
        ind->parts[i].pmin_size = (int) (part_size * (1.0 - off_ratio) - 1.0);
    }   /* for i */ 

    free(tparts);

    return max_size;
}   /* create_partition */

/* copy partition properties to parts_info for temporary use */
void copy_partition(int noparts,
                    parts_info_t parts_info[],
                    ind_t *ind)
{
    for (int i = 0; i < noparts; i++) {
        parts_info[i].pmax_cells = ind->parts[i].pmax_cells;
        parts_info[i].pmin_size = ind->parts[i].pmin_size;
        parts_info[i].pcurr_size = ind->parts[i].pcurr_size;
        parts_info[i].pmax_size = ind->parts[i].pmax_size;
    }   /* for i */
}   /* copy_partition */

/* read a partition prepared beforehand */
int read_partition(FILE *fp,
                   char *filename,
                   int noparts,
                   ind_t *ind)
{
    int max_size = -1; 
    open_file(&fp, filename, "r");

    for (int i = 0; i < noparts; i++) {

        int part_no;
        if (fscanf(fp, "%d", &part_no) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", filename, errno, strerror(errno));
        }

        if (fscanf(fp, "%d%d%d%d",
                   &(ind->parts[part_no].pmax_cells),
                   &(ind->parts[part_no].pmin_size),
                   &(ind->parts[part_no].pcurr_size),
                   &(ind->parts[part_no].pmax_size)) == EOF) {
            printf("Error: Cannot read from %s: errno= %d error= %s\n", filename, errno, strerror(errno));        
        }

        if (ind->parts[part_no].pmax_size > max_size)
            max_size = ind->parts[part_no].pmax_size;

        for (int j = 0; j < ind->parts[part_no].pmax_cells; j++) {

            int cell_no;
            if (fscanf(fp, "%d", &cell_no) == EOF) {
                printf("Error: Cannot read from %s: errno= %d error= %s\n", filename, errno, strerror(errno));
            }
            ind->chrom[cell_no] = part_no; 

        }   /* for j */

    }   /* for i */

    close_file(&fp);

    return max_size;
}   /* read_partition */

/* write a partition */
void write_partition(FILE *fp,
                     char *filename,
                     int nocells,
                     int noparts,
                     ind_t *ind)
{
    open_file(&fp, filename, "w");

    for (int i = 0; i < noparts; i++) {
        fprintf(fp, "%d %d %d %d %d\n", i,
                ind->parts[i].pmax_cells,
                ind->parts[i].pmin_size,
                ind->parts[i].pcurr_size,
                ind->parts[i].pmax_size);
        for (int j = 0; j < nocells; j++)
            if (ind->chrom[j] == i)
                fprintf(fp, "%d ", j);
        fprintf(fp, "\n");
    }   /* for i */

    close_file(&fp);
}   /* read_partition */

/* EOF */
