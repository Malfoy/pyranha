CC=/usr/bin/g++
# CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread -pipe
LDFLAGS=-pthread -pipe -flto


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif

ifeq ($(polly),1)
CFLAGS+= -Xclang -load -Xclang /local/polly/llvm_build/lib/LLVMPolly.so -mllvm -polly -mllvm -polly-vectorizer=stripmine
CC=/local/polly/llvm_build/bin/clang
endif


EXEC=ftl

all: $(EXEC)

ftl: ftl.o  utils.o mapping.o distances.o ssw_cpp.o ssw.o
	$(CC) -o $@ $^ $(LDFLAGS)

ftl.o: main.cpp  utils.h mapping.h
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

distances.o: distances.cpp distances.h
	$(CC) -o $@ -c $< $(CFLAGS)

mapping.o: mapping.cpp utils.h mapping.h distances.h
	$(CC) -o $@ -c $< $(CFLAGS)

ssw_cpp.o: ssw_cpp.cpp ssw_cpp.h ssw.h
	$(CC) -o $@ -c $< $(CFLAGS)

ssw.o: ssw.c ssw.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
