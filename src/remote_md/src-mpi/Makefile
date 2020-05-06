# Makefile for src-parallel

SHELL = /bin/bash

### CoMD can be built in either double or single precision and with or
### without MPI.  Select desired precision and MPI here.

# double precision (ON/OFF)
DOUBLE_PRECISION = ON
# MPI for parallel (ON/OFF)
DO_MPI = ON


### Set your desired C compiler and any necessary flags.  Note that CoMD
### uses some c99 features.  You can also set flags for optimization and
### specify paths to include files that the compiler can't find on its
### own.  If you need any -L or -l switches to get C standard libraries
### (such as -lm for the math library) put them in C_LIB.
CC = mpicc
CFLAGS = -std=c99
OPTFLAGS = -g -O5
INCLUDES = 
C_LIB = -lm


### If you need to specify include paths, library paths, or link flags
### for MPI, put them here.  Put both -L and -l switches into MPI_LIB.
MPI_LIB =
MPI_INCLUDE =

### A place to specify any other include or library switches your
### platform requires.
margo_dir = /g/g90/peng8/flash/deps
umap_dir  = /g/g90/peng8/flash/forked_umap_remote_region/build5

OTHER_LIB = -L$(margo_dir)/lib -L${umap_dir}/lib -lumap -lmargo
OTHER_INCLUDE = -I$(margo_dir)/include -I${umap_dir}/include




#########################################
### Below here, it is pitch black.  
### You are likely to be eaten by a grue.
##########################################

# clear all suffixes
.SUFFIXES:
# list only those that we use 
.SUFFIXES: .c .o

.PHONY: DEFAULT clean distclean depend

BIN_DIR=../bin

# Check for double precision
ifeq ($(DOUBLE_PRECISION), ON)
CFLAGS += -DDOUBLE
else
CFLAGS += -DSINGLE
endif

# Set executable name and add includes & libraries for MPI if needed.
CoMD_VARIANT = CoMD-serial
ifeq ($(DO_MPI), ON)
CoMD_VARIANT = CoMD-mpi
INCLUDES += ${MPI_INCLUDE}
CFLAGS += -DDO_MPI
LDFLAGS += ${MPI_LIB}
endif
CoMD_EXE = ${BIN_DIR}/${CoMD_VARIANT}

LDFLAGS += ${C_LIB} ${OTHER_LIB}
CFLAGS  += ${OPTFLAGS} ${INCLUDES} ${OTHER_INCLUDE}


SOURCES := $(wildcard *.c)
SOURCES := $(filter-out client.c CoMD.c, $(SOURCES) )
OBJECTS=$(SOURCES:.c=.o)


DEFAULT: client ${CoMD_EXE}

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

${CoMD_EXE}: ${BIN_DIR} CoMD_info.h CoMD.c ${OBJECTS} 
	${CC} ${CFLAGS} -o ${CoMD_EXE} CoMD.c ${OBJECTS} ${LDFLAGS}

client: ${BIN_DIR} CoMD_info.h client.c ${OBJECTS} 
	${CC} ${CFLAGS} -o $@ client.c ${OBJECTS} ${LDFLAGS}

CoMD_info.h: Makefile
	./generate_info_header ${CoMD_VARIANT} "$(CC)" "$(CFLAGS)" "$(LDFLAGS)"

${BIN_DIR}:
	@if [ ! -d ${BIN_DIR} ]; then mkdir -p ${BIN_DIR} ; fi

clean:
	rm -f *.o CoMD_info.h .depend client 

distclean: clean
	rm -f ${CoMD_EXE} .depend.bak
	rm -rf html latex

.depend: $(SOURCES)
	@touch .depend
	@$(MAKE) --no-print-directory depend

depend:
	@echo "Rebuilding dependencies..."
	@$(MAKE) CoMD_info.h
	@makedepend -f .depend -Y. --$(CFLAGS)-- $(SOURCES) 2> /dev/null


-include .depend
