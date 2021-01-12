#!/bin/bash

TARGETS=lupp

XBLASPATH=../lib/xblas-1.0.248

CROSS_TOOL=
CC_CPP=$(CROSS_TOOL)g++
CC_C=$(CROSS_TOOL)gcc

CFLAGS = -O3 -Wall -g 

INCLUDEPATH = -I./ -I${XBLASPATH}/src -I${MKLROOT}/include 
LIBPATH = -L./ -L/usr/local/lib -L${XBLASPATH} -L${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/ 

LIBS = -lxblas ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/../compiler/lib/intel64/libiomp5.a -lpthread -ldl

all: clean $(TARGETS) 

$(TARGETS):
	$(CC_CPP) $(CFLAGS) $(INCLUDEPATH) $@.cpp -o $@ $(LIBPATH) $(LIBS)	

clean:
	rm -f $(TARGETS)