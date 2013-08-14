BIN = .
CPP = clang++
CFLAGS = -Wall -std=c++0x -stdlib=libc++ -I/usr/local/include
LIBS = -L/usr/local/lib

INPUT = $(wildcard *.cpp)
OBJS  = $(INPUT:.cpp=.o)
TARGET = $(BIN)/convert

all: $(TARGET) 

include Makefile.deps

$(BIN)/convert: convert.o
	$(CPP) $(CFLAGS) $(LIBS) convert.o -o $@

dep: Makefile.deps

Makefile.deps: $(INPUT)
	$(CPP) $(CFLAGS) -MM $< > $@

clean:
	rm *.o

.SUFFIXES: .o .cpp

.cpp.o:
	$(CPP) $(CFLAGS) -O2 -c $< -o $@

.PHONY: all dep clean
