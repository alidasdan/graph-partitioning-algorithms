CC0 = gcc
CC = $(CC0)
LD = $(CC0)
FLAGS1 = -O3
FLAGS2 = -lm
FLAGS3 = -std=c99 -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
FLAGS = $(FLAGS1) $(FLAGS3)

DEPS_SHARE = ../share/ad_bucketio.h \
	../share/ad_fileio.h \
	../share/ad_lib.h \
	../share/ad_partition.h \
	../share/ad_print.h \
	../share/ad_random.h \
	../share/ad_readinput.h 
DEPS_FMS = ad_lib_fms.h ../share/ad_defs.h $(DEPS_SHARE)
OBJS = $(patsubst %.h,%.o,$(DEPS_SHARE))
OBJS_FMS = $(OBJS) ad_lib_fms.o

all: ad_fms

ad_lib_fms.o: ad_lib_fms.c ad_lib_fms.h ../share/ad_lib.h ../share/ad_bucketio.h ../share/ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_fms: ad_fms.c $(DEPS_FMS) $(OBJS_FMS)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_FMS)

# Testing:
test:
	./utest.sh

# Cleaning:
clean c cl cle clea: 
	rm -f *.o *~ core *.x
cleano: 
	rm -f *.o *~ core

# End of file
