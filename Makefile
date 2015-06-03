TARGET=bycatch
OBJS=bycatch.o
SAMDIR=../samtools
HTSDIR=../htslib
CFLAGS=-I$(HTSDIR) -I$(SAMDIR) -Wall -O2 -ggdb3 -std=c99
LDFLAGS=-L$(HTSDIR) -L$(SAMDIR)
LIBS=-lbam -lhts -lm

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET) $(OBJS)
