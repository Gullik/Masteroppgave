##
##      MG test_environment Makefile
##      Author: Gullik V. Killie
##      University of Oslo, Norway 
##      2015
##

CC = mpicc
EXEC = MGSOLVER

CFLAGS= -UDEBUG -DSTANDALONE -UCLOCK -DPERIODIC\
## Flags
# -DMPI 		~Must be switched on to use MPI properly


DEPS = const.h nrutil.h funct.h

LIBS = -lm

OBJECTS= main.o \
		shortcut.o \
		input.o \
		dustg.o \
		grid.o \
		./mg/fmg_P.o \
		./mg/nrutil.o \


$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) %.c

clean:
	rm -f *~ \#* $(MGOBJECTS) $(EXEC)