
/* COPYRIGHT C 1991- Ali Dasdan */ 

#include <time.h>
#include <stdlib.h>
#include "ad_random.h"

#define True              1    /* boolean true */
#define False             0    /* boolean false */

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

/* generates a float random number in [min, max] */
float frandom(float min, float max)
{
    if (min >= max) { 
        return min;
    } else { 
        return (float) (rand01() * (max - min) + (float) min);
    }
}     /* frandom */ 

/* generates a boolean true value with specified probability */
int flip(float prob)
{
    if (prob >= 1.0) {
        return True;
    } else if (prob <= 0.0) {
        return False; 
    } else if (rand01() <= prob) {
        return True;
    } else {
        return False; 
    }
}     /* flip */

/* EOF */
