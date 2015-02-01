CC0 = gcc
CC = $(CC0)
LD = $(CC0)
FLAGS1 = -O3
FLAGS2 = -lm
FLAGS3 = -std=c99 -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
DIR_H = ../share
FLAGS = $(FLAGS1) $(FLAGS3) -I$(DIR_H)

DEPS_SHARE = $(DIR_H)/ad_bucketio.h \
	$(DIR_H)/ad_fileio.h \
	$(DIR_H)/ad_lib.h \
	$(DIR_H)/ad_partition.h \
	$(DIR_H)/ad_print.h \
	$(DIR_H)/ad_random.h \
	$(DIR_H)/ad_readinput.h 
DEPS_PFM = ad_lib_pfm.h $(DIR_H)/ad_defs.h $(DEPS_SHARE)
OBJS = $(patsubst %.h,%.o,$(DEPS_SHARE))
OBJS_PFM = $(OBJS) ad_lib_pfm.o

all: ad_pfm

ad_lib_pfm.o: ad_lib_pfm.c ad_lib_pfm.h $(DIR_H)/ad_lib.h $(DIR_H)/ad_bucketio.h $(DIR_H)/ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $< $(FLAGS2)

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
