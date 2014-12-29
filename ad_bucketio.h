#ifndef AD_BUCKETIO_INCLUDED
#define AD_BUCKETIO_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* create a partb node */
void create_partb_node(int noparts,
                       int mov_cell_no, 
                       int home_part, 
                       int mapped_dest_part, 
                       int mov_gain_inx,
                       partb_t partb[][noparts - 1], 
                       cells_info_t cells_info[]);

/* insert a node into partb */
void insert_partb_node(bnode_ptr_t tnode_ptr,
                       int mapped_dest_part, 
                       int mov_gain_inx,
                       partb_t *partb_ptr,
                       cells_info_t *cells_info_ptr);

/* delete all nodes of a cell from partb */
void delete_partb_nodes_of_cell(int noparts, 
                                int mov_cell_no, 
                                int home_part,
                                partb_t partb[][noparts - 1],
                                cells_info_t cells_info[]);

/* delete a partb node */
bnode_ptr_t delete_partb_node(int deletion_ok, 
                              int mapped_dest_part,
                              partb_t *partb_ptr,
                              cells_info_t *cells_info_ptr);

#endif
