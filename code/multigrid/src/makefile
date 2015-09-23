##
##      MG test_environment Makefile
##      Author: Gullik V. Killie
##      University of Oslo, Norway 
##      2015
##

CC = mpicc
EXEC = MGSOLVER

CFLAGS_TRI=-UDEBUG -DSTANDALONE -UCLOCK\
	   -I$(inclpath) -L/lib/pa1.1 -g -I.

DEPS = const.h nrutil.h funct.h


LIBS = -lm

OBJECTS= main.o \
		shortcut.o \
		input.o \
		dustg.o \
		./mg/fmg_P.o \
		./mg/nrutil.o \




$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

%.o: %.c $(DEPS)
	$(CC) -c %.c



clean:
	rm -f *~ \#* $(MGOBJECTS) $(EXEC)