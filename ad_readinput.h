#ifndef AD_READINPUT_INCLUDED
#define AD_READINPUT_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

/* initialize cells array */
void init_cells(int nocells, cells_t cells[]);

/* initialize netlist pointers */
void init_netlist(int nonets, 
                  cells_t cells[], 
                  nets_t nets[], 
                  corn_t cnets[]);

/* read input graph size to allocate memory */
void read_graph_size(char fname[],
                     int  *nocells,
                     int  *nonets);

/* read input graph and construct cell, net arrays */
void read_graph(char fname[], 
                int  nocells, 
                int  nonets, 
                int  noparts, 
                int *totsize, 
                int *totcellsize,
                int *max_density, 
                int *max_cweight, 
                int *max_nweight,
                cells_t cells[], 
                nets_t  nets[], 
                corn_t  cnets[]);

#endif
