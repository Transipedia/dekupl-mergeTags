CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function -std=c99
HEADERS=kstring.h
OBJECTS=$(HEADERS:.h=.o)
LIBS=-lz

all:mergeTags

mergeTags: mergeTags.c $(HEADERS) $(OBJECTS) dna.h
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)
