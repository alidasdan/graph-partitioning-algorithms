#ifndef AD_FILEIO_INCLUDED
#define AD_FILEIO_INCLUDED

/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdio.h>

/* open file fp with filename = fname and mode = mode */
void open_file(FILE **fp, 
               char *fname, 
               char *mode);
 
/* close file fp */
void close_file(FILE **fp);

#endif
