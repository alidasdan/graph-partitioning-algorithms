
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "ad_fileio.h"

/* open file fp with filename = fname and mode = mode */
void open_file(FILE **fp, 
               char *fname, 
               char *mode)
{
    if ((*fp = fopen(fname, mode)) == NULL) {
        printf("Error: File %s can NOT be opened with mode %s: errno= %d error= %s\n", 
               fname, mode, errno, strerror(errno));
        exit(1);
    }
}     /* open_file */
 
/* close file fp */
void close_file(FILE **fp)
{
    if (fclose(*fp) != 0) {
        printf("Error: Cannot close file: errno= %d error= %s\n", errno, strerror(errno));
        exit(1);
    }
}     /* close_file */

/* EOF */
