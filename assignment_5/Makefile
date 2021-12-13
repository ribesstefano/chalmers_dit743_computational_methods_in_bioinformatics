# Author: Stefano Ribes
CC = gcc
CFLAGS = -g -std=c99 -O3 -fopenmp # -Wall 
CXX = g++
CXXFLAGS = -g -std=c++11 -O3 -fopenmp # -Wall 
LDFLAGS = -lm
DEFINES = -DNO_LOOKUP_TABLE=1

all: detect_steric_clashes.exe

SRC = $(wildcard *.cc)
CSRC = $(wildcard *.c)
COBJS = $(CSRC:.c=.o)

%.o : %.c
	$(CC) -c $(CFLAGS) $(DEFINES) $< -o $@ $(LDFLAGS)

detect_steric_clashes.exe: $(SRC) $(COBJS)
	$(CXX) $(CXXFLAGS) $(DEFINES) -I. -o detect_steric_clashes.exe $(SRC) $(COBJS) $(LDFLAGS)
	# ./detect_steric_clashes.exe 1cdh.pdb 2csn.pdb

clean:
	rm -f *.exe *.o *.ps