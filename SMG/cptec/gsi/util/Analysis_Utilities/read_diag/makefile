#! /bin/sh -v

# Define the fortran compiler and options

include ../../../configure.gsi

FC   = $(SFC)
FLAGS= $(FFLAGS_DEFAULT)   
INCLD=  -I../../../include
LIBS =  -L../../../src/main -lgsi

all: read_diag_conv read_diag_rad


read_diag_conv: read_diag_conv.o
	${FC} -o read_diag_conv.exe ${FLAGS} read_diag_conv.o ${LIBS} 
	rm read_diag_conv.o

read_diag_rad: read_diag_rad.o
	${FC} -o read_diag_rad.exe ${FLAGS} read_diag_rad.o ${LIBS} 
	rm read_diag_rad.o

.SUFFIXES : .f90  .o

.f90.o :
	${FC} ${FLAGS} ${INCLD} -c $<

clean:
	/bin/rm -f *.o *.exe
