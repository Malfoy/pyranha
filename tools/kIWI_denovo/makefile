CC=/usr/bin/g++
#CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread -fmax-errors=3 -flto -pipe
LDFLAGS=-pthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=kIWI

all: $(EXEC)

kIWI:	main.o	kIWI.o
	$(CC) -o $@ $^ $(LDFLAGS)

kIWI.o: 	kIWI.cpp	kIWI.h
	$(CC) -o $@ -c $< $(CFLAGS)

main.o:	main.cpp	kIWI.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
