BIN = ../bin
CPP = clang++
CFLAGS = -Wall -std=c++0x -stdlib=libc++ -I/usr/local/include
LIBS = -L/usr/local/lib -lboost_program_options

INPUT = $(wildcard *.cpp)
OBJS  = $(INPUT:.cpp=.o)
TARGET = $(BIN)/bisearch $(BIN)/tm

all: $(TARGET) 

include Makefile.deps

$(BIN)/bisearch: bisearch.o primer.o main.o
	$(CPP) $(CFLAGS) $(LIBS) bisearch.o primer.o main.o -o $@

$(BIN)/tm: tm.o
	$(CPP) $(CFLAGS) $(LIBS) tm.o bisearch.o -o $@

dep: Makefile.deps

Makefile.deps: $(INPUT)
	$(CPP) $(CFLAGS) -MM $< > $@

clean:
	rm *.o

.SUFFIXES: .o .cpp

.cpp.o:
	$(CPP) $(CFLAGS) -O2 -c $< -o $@

.PHONY: all dep clean
