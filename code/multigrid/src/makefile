##
##      MG test_environment Makefile
##      Author: Gullik V. Killie
##      University of Oslo, Norway 
##      2015
##

CC = gcc
EXEC = MGSOLVER

CFLAGS_TRI=-UDEBUG -DSTANDALONE -UCLOCK\
	   -I$(inclpath) -L/lib/pa1.1 -g -I.
inclpath = ./mg

LIBS = -lm

OBJECTS= main.o \
		./mg/hellofunc.o \
		./mg/fmg_P.o \


$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJECTS) -I.

main.o: main.c const.h
	$(CC) -c main.c

hellofunc.o: ./mg/hellofunc.c
	$(CC) -c ./mg/hellofunc.c

fmg_P.o:	./mg/fmg_P.c const.h
		$(CC) -c $(CFLAGS) ./mg/fmg_P.c


clean:
	rm -f *~ \#* $(MGOBJECTS) $(EXEC)