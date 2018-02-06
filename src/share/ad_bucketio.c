
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>
#include <stdlib.h>
#include <malloc/malloc.h>
#include "ad_defs.h"
#include "ad_bucketio.h"

/* create a partb node */
void create_partb_node(int noparts,
                       int mov_cell_no, 
                       int home_part, 
                       int mapped_dest_part, 
                       int mov_gain_inx,
                       partb_t partb[][noparts - 1], 
                       cells_info_t cells_info[])
{
    bnode_ptr_t tnode_ptr = (bnode_ptr_t) malloc(sizeof(bnode_t));
    if (tnode_ptr == NULL) {
        printf("Error: Cannot allocate memory for tnode_ptr.\n");
        exit(1);
    }   /* if */
    tnode_ptr->cell_no = mov_cell_no;
    insert_partb_node(tnode_ptr, 
                      mapped_dest_part, mov_gain_inx,
                      &partb[home_part][mapped_dest_part], 
                      &cells_info[mov_cell_no]);
}   /* create_partb_node */

/* insert a node into partb */
void insert_partb_node(bnode_ptr_t tnode_ptr,
                       int mapped_dest_part, 
                       int mov_gain_inx,
                       partb_t *partb_ptr,
                       cells_info_t *cells_info_ptr)
{
    /*
      partb_ptr = partb[home_part][mapped_dest_part]
      cells_info_ptr = cells_info[mov_cell_no]
    */
    int bucket_empty = False;
    if (partb_ptr->bnode_ptr[mov_gain_inx] == NULL) {
        bucket_empty = True;
    }
    tnode_ptr->lptr = NULL;
    tnode_ptr->rptr = partb_ptr->bnode_ptr[mov_gain_inx];
    if (! bucket_empty) {
        (tnode_ptr->rptr)->lptr = tnode_ptr;
    }
    partb_ptr->bnode_ptr[mov_gain_inx] = tnode_ptr;

    /* insert into cells_info */
    cells_info_ptr->partb_ptr[mapped_dest_part] = tnode_ptr;
    cells_info_ptr->partb_gain_inx[mapped_dest_part] = mov_gain_inx;

    /* update indices */
    if (partb_ptr->nobuckets == 0) {
        partb_ptr->max_inx = mov_gain_inx;
        partb_ptr->min_inx = mov_gain_inx;
    } else if (mov_gain_inx > partb_ptr->max_inx) {
        partb_ptr->max_inx = mov_gain_inx;
    } else if (mov_gain_inx < partb_ptr->min_inx) {
        partb_ptr->min_inx = mov_gain_inx;
    }
    if (bucket_empty) {
        partb_ptr->nobuckets++;
    }
}   /* insert_partb_node */

/* delete all nodes of a cell from partb */
void delete_partb_nodes_of_cell(int noparts, 
                                int mov_cell_no, 
                                int home_part,
                                partb_t partb[][noparts - 1],
                                cells_info_t cells_info[])
{
    for (int mapped_dest_part = 0; mapped_dest_part < (noparts - 1); mapped_dest_part++) {
        delete_partb_node(True, mapped_dest_part, 
                          &partb[home_part][mapped_dest_part],
                          &cells_info[mov_cell_no]);
    }
}   /* delete_partb_nodes_of_cell */

/* delete a partb node */
bnode_ptr_t delete_partb_node(int deletion_ok, 
                              int mapped_dest_part,
                              partb_t *partb_ptr,
                              cells_info_t *cells_info_ptr)
{
    /*
      cells_info_ptr = cells_info[mov_cell_no]
      partb_ptr = partb[home_part][mapped_dest_part]
    */
    bnode_ptr_t tnode_ptr = cells_info_ptr->partb_ptr[mapped_dest_part];
    int mov_gain_inx = cells_info_ptr->partb_gain_inx[mapped_dest_part];
    if (tnode_ptr != NULL) {
        if (tnode_ptr->lptr == NULL) {  /* if 1st node */
            if (tnode_ptr->rptr != NULL) {  /* make next node the 1st node */
                (tnode_ptr->rptr)->lptr = NULL;
            }
            partb_ptr->bnode_ptr[mov_gain_inx] = tnode_ptr->rptr;
            if (partb_ptr->bnode_ptr[mov_gain_inx] == NULL) {
                partb_ptr->nobuckets--;
                if (partb_ptr->nobuckets == 0) { /* if no more full buckets */
                    partb_ptr->max_inx = -1;
                    partb_ptr->min_inx = -1;
                } else {  /* set max_inx to the next full bucket index */
                    while ((partb_ptr->max_inx > partb_ptr->min_inx) &&
                           (partb_ptr->bnode_ptr[partb_ptr->max_inx] == NULL)) {
                        partb_ptr->max_inx--;
                    }  /* while */
                }   /* else */
            }   /* if bucket is empty */
        } else {
            (tnode_ptr->lptr)->rptr = tnode_ptr->rptr;
            if (tnode_ptr->rptr != NULL) {
                (tnode_ptr->rptr)->lptr = tnode_ptr->lptr;
            }
        }   /* else */

        if (deletion_ok) {
            free(tnode_ptr);
            cells_info_ptr->partb_ptr[mapped_dest_part] = NULL;
        }   /* if */
    }   /* if tnode_ptr is not nil */

    return (tnode_ptr);
}   /* delete_partb_node */

/* EOF */
