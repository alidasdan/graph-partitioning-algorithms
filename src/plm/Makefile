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
DEPS_PLM = ad_lib_plm.h ../share/ad_defs.h $(DEPS_SHARE)
OBJS = $(patsubst %.h,%.o,$(DEPS_SHARE))
OBJS_PLM = $(OBJS) ad_lib_plm.o

all: ad_plm

ad_lib_plm.o: ad_lib_plm.c ad_lib_plm.h ../share/ad_lib.h ../share/ad_bucketio.h ../share/ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $<

ad_plm: ad_plm.c $(DEPS_PLM) $(OBJS_PLM)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_PLM)

# Testing:
test:
	./utest.sh

# Cleaning:
clean c cl cle clea: 
	rm -f *.o *~ core *.x
cleano: 
	rm -f *.o *~ core

# End of file
