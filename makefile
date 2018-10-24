# I am a comment, and I want to say that the variable CC will be
# the compiler to use and CFLAGS are some flags.
CC=g++
CFLAGS=-c -g -O3 -fopenmp
PAT=Src


all: MAIN

MAIN: maincode.o memfuncs.o mersenne.o stoc1.o userintf.o solfuncs.o poly34.o
	$(CC) -fopenmp maincode.o memfuncs.o mersenne.o stoc1.o userintf.o solfuncs.o poly34.o -o MAIN

maincode.o: maincode.cpp
	$(CC) $(CFLAGS) maincode.cpp
memfuncs.o : $(PAT)/memfuncs.cpp
	$(CC) $(CFLAGS) $(PAT)/memfuncs.cpp
mersenne.o : $(PAT)/mersenne.cpp
	$(CC) $(CFLAGS) $(PAT)/mersenne.cpp
stoc1.o : $(PAT)/stoc1.cpp
	$(CC) $(CFLAGS) $(PAT)/stoc1.cpp
userintf.o : $(PAT)/userintf.cpp
	$(CC) $(CFLAGS) $(PAT)/userintf.cpp
solfuncs.o : $(PAT)/solfuncs.cpp
	$(CC) $(CFLAGS) $(PAT)/solfuncs.cpp
poly34.o : $(PAT)/poly34.cpp
	$(CC) $(CFLAGS) $(PAT)/poly34.cpp	

clean:
	rm -rf *o MAIN

