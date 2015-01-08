CC0 = gcc
CC = $(CC0)
LD = $(CC0)
FLAGS1 = -O3
FLAGS2 = -lm
FLAGS3 = -std=c99 -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
FLAGS = $(FLAGS1) $(FLAGS3)

DEPS_COMMON = ad_bucketio.h ad_fileio.h ad_lib.h ad_partition.h ad_print.h ad_random.h ad_readinput.h
OBJS = $(patsubst %.h,%.o,$(DEPS_COMMON))

all: $(OBJS)

ad_bucketio.o: ad_bucketio.c ad_bucketio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_fileio.o: ad_fileio.c ad_fileio.h 
	$(CC) $(FLAGS) -c -o $@ $<

ad_lib.o: ad_lib.c ad_lib.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_partition.o: ad_partition.c ad_partition.h ad_defs.h ad_fileio.h ad_random.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_print.o: ad_print.c ad_print.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_random.o: ad_random.c
	$(CC) $(FLAGS) -c -o $@ $<

ad_readinput.o: ad_readinput.c ad_readinput.h ad_fileio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

# Testing:
test:
	./utest.sh

# Cleaning:
clean c cl cle clea: 
	rm -f *.o *~ core *.x
cleano: 
	rm -f *.o *~ core

# End of file
