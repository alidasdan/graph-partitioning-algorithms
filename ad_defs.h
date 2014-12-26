#ifndef AD_DEFS_INCLUDED
#define AD_DEFS_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* naming conventions :
   cell = node, net = edge,
*/

#define True              1    /* boolean true */
#define False             0    /* boolean false */

#define MAX_POP            1     /* max numbr of individuals in population */
#define NIL               -1     /* point to nowhere */
#define STR_SIZE          30     /* string size */
#define NO_ITERATIONS     200    /* max number of iterations for KL */
#define EPSILON           0.01   /* min move value of a cell */

/* cells type */
typedef struct cells_st {
    int  cno_nets;           /* total number of nets on the cell */
    int  cno_inets;          /* number of internal nets */
    int  cweight;            /* weight of cell */
    int  netlist;            /* pointer to nets on the cell */
}   cells_t;
 
/* nets type */
typedef struct nets_st {
    int nweight;            /* weight of node */
    int ncells[2];          /* For GP, nno_cells = 2 and */
    /* ncells[0] & ncells[1] are the cells on the net */
}   nets_t;
 
/* cell cornlist type - CorN = Cell or Net */
typedef struct corn_st {
    int corn_no;            /* cell or net number */
}   corn_t;

/* allele definition - for compatibility purposes */
typedef int allele;

/* partition type */
typedef struct parts_st {
    int   pmax_cells;     /* maximum number of cells in partition */
    int   pmax_size;      /* maximum size of part */
    int   pmin_size;      /* minimum size of part */
    int   pcurr_size;     /* current size of part */
    float pratio;         /* pmax_size / totsize */
}   parts_t;

/* individual type */
typedef struct pop_st {
    allele  *chrom;    /* string holding partitions */
    parts_t *parts;    /* partition array */
    int      incost;   /* sum of net weights - cut cost */
}   ind_t;
 
/* temporary partition type */
typedef struct tparts_st {  /* used while creating partition celllists */
    char filled;            /* set if partition is filled with cells */
    int  pcells_inx;        /* an index to pcells */
}   tparts_t;

/* bucket nodes (= bnode) */
typedef struct bnode_st *bnode_ptr_t;

typedef struct bnode_st {
    int         cell_no; /* cell in the node */
    bnode_ptr_t lptr;    /* pointer to left node */
    bnode_ptr_t rptr;    /* pointer to right node */
}   bnode_t;

/* partition bucket (partb) type */
typedef struct partb_st {
    bnode_ptr_t *bnode_ptr;   
    int          max_inx;        /* max index of filled bucket cell */
    int          min_inx;        /* min index of filled bucket cell */
    int          nobuckets;      /* number of filled buckets */
}   partb_t;

/* additional information for cells */
typedef struct cells_info_st {
    int          mcount;          /* count of moves that cell performs */
    int          locked;          /* set if cell is locked */
    int         *mgain;           /* external costs of moving cell_no to all parts */ 
    /* only mgain[part_no of cell_no] is internal cost */
    bnode_ptr_t *partb_ptr;       /* pointer to partb node */
    int         *partb_gain_inx;  /* index to partb bucket */
    /* pointers to part_bucket */
}   cells_info_t; 

/* selected cell structure */
typedef struct selected_cell_st {
    int mov_cell_no;        /* current properties */
    int from_part;
    int to_part;
    int mov_gain;
}   selected_cell_t;

/* moved cells array */
typedef struct mcells_st {
    int cell_no;     /* the cell moved */
    int from;        /* cell's home partition */
    int to;          /* cell's destination partition */
    int mgain;       /* move gain */
}  mcells_t;

/* partition information type */
typedef struct parts_info_st {
    int  pmax_cells;       /* maximum number of cells in partition */
    int  pmax_size;        /* maximum size of part */
    int  pmin_size;        /* minimum size of part */
    int  pcurr_size;       /* current size of part */
}   parts_info_t;

/* precomputed exponential values */
typedef struct eval_st {
    float val;
}   eval_t;

#endif
