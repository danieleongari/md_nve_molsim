CC          = gcc
CLINKER     = gcc

CFLAGS      =   -Wall -O4 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = force.c integrate.c md.c mdloop.c readdat.c sample_diff.c sample_rdf.c writepdb.c system.c init.c ran_uniform.c
OBJS    = force.o integrate.o md.o mdloop.o readdat.o sample_diff.o sample_rdf.o writepdb.o system.o init.o ran_uniform.o
EXECS   = md_nve

default: md_nve

all: $(EXECS)

md_nve:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o md_nve $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

force.o: system.h
init.o: system.h
integrate.o: system.h
md.o: system.h
mdloop.o: system.h
ran_uniform.o: system.h
readdat.o: system.h
sample_diff.o: system.h
sample_rdf.o: system.h
system.o: system.h
writepdb.o: system.h
