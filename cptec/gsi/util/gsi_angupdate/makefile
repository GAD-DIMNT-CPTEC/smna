#! /bin/sh -v

# Define the fortran compiler and options

include ../../configure.gsi

FC   =  ${DM_F90}       
FFFLAGS=  -O3 $(FFLAGS_DEFAULT) $(LDFLAGS) -DLINUX
INCLD= -I./  -I$(INC_DIR)
LIBS =  -L../../lib -lw3nco_i4r8 ${MYLIBsys} -L../../src/main -lgsi

OBJS = determine_time_levels.o kinds.o read_diag.o

all: gsi_angupdate

gsi_angupdate: ${OBJS} main.o
	${FC} -o gsi_angupdate.exe ${FFFLAGS} main.o ${OBJS} ${LIBS} 

# ------------
#  # Dependencies
# ------------
    MAKE_DEPEND = Makefile.dependency
include $(MAKE_DEPEND)
#
#
.SUFFIXES : .f90 .F90 .o

.f90.o :
	${FC} ${FFFLAGS} ${INCLD} -c $<

.F90.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPP_FLAGS) $< > $*.f90
	$(F90) $(FFFLAGS) ${INCLD} -c $*.f90
	$(RM) $*.f90


clean:
	/bin/rm -f *.o *.exe *.mod
