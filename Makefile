CC0 = gcc
CC = $(CC0)
LD = $(CC0)
FLAGS1 = -O3
FLAGS2 = -lm
FLAGS3 = -std=c99 -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
FLAGS = $(FLAGS1) $(FLAGS3)

DEPS_COMMON = ad_bucketio.h ad_fileio.h ad_lib.h ad_partition.h ad_print.h ad_random.h ad_readinput.h 
DEPS_FMS = ad_lib_fms.h ad_defs.h $(DEPS_COMMON)
DEPS_PLM = ad_lib_plm.h ad_defs.h $(DEPS_COMMON)
DEPS_PFM = ad_lib_pfm.h ad_defs.h $(DEPS_COMMON)
OBJS = $(patsubst %.h,%.o,$(DEPS_COMMON))
OBJS_FMS = $(OBJS) ad_lib_fms.o
OBJS_PLM = $(OBJS) ad_lib_plm.o
OBJS_PFM = $(OBJS) ad_lib_pfm.o

all: ad_fms ad_plm ad_pfm

ad_bucketio.o: ad_bucketio.c ad_bucketio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_fileio.o: ad_fileio.c ad_fileio.h 
	$(CC) $(FLAGS) -c -o $@ $<

ad_lib.o: ad_lib.c ad_lib.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_lib_fms.o: ad_lib_fms.c ad_lib_fms.h ad_lib.h ad_bucketio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_lib_plm.o: ad_lib_plm.c ad_lib_plm.h ad_lib.h ad_bucketio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_lib_pfm.o: ad_lib_pfm.c ad_lib_pfm.h ad_lib.h ad_bucketio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $< $(FLAGS2)

ad_partition.o: ad_partition.c ad_partition.h ad_defs.h ad_fileio.h ad_random.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_print.o: ad_print.c ad_print.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_random.o: ad_random.c
	$(CC) $(FLAGS) -c -o $@ $<

ad_readinput.o: ad_readinput.c ad_readinput.h ad_fileio.h ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_fms: ad_fms.c $(DEPS_FMS) $(OBJS_FMS)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_FMS)

ad_plm: ad_plm.c $(DEPS_PLM) $(OBJS_PLM)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_PLM)

ad_pfm: ad_pfm.c $(DEPS_PFM) $(OBJS_PFM)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_PFM) $(FLAGS2)

# Testing:
test:
	./utest.sh

# Cleaning:
clean c cl cle clea: 
	rm -f *.o *~ core *.x
cleano: 
	rm -f *.o *~ core

# End of file
