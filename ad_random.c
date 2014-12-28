
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <time.h>
#include <stdlib.h>
#include "ad_random.h"

/* initializes random number generator with seed or */
/* with any value if seed = -1 */
long randomize(long seed)
{
    long in_seed;

    if (seed == -1) {
        time(&in_seed);  /* init with current time */
    } else { 
        in_seed = seed;   /* init with seed */
    }
    srand48(in_seed);

    return in_seed;
}     /* randomize */

/* generates a float random number in [0, 1) */
float rand01()
{
    return (float) drand48();
}     /* rand01 */

/* generates a int random number in [min, max] */
/* if int output is needed, prms are true min & max */
int irandom(int min, int max)
{
    if (min >= max) {
        return min;
    } else {
        int retval = (int) (rand01() * (float) (max + 1 - min) + (float) min);
        if (retval > max) {
            retval = max;
        }
        return retval;
    }   /* else */
}     /* irandom */ 

/* EOF */
