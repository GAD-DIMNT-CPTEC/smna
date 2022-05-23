# Automatically generated GNU Makefile. Dom Jul 27 11:51:57 BRT 2014
# fgen v0.3 (C) 1997,98 Beroud Jean-Marc

#include ../../../configure.gsi

# simultaneous parallel jobs & load average limit
MAXJOBS = 1
MAXLOAD = 1


# utils
ECHO    = echo
RM      = rm
CP      = cp
MV      = mv
CD      = cd
MKDIR   = mkdir

# dirs
srcdir = ./
LIBDIR = ../../../lib
INC_DIR= ../../../include
FFLAGS_SIGIOBAM=-DuseMPI
# extra flags
SF90 = ftn
FFLAGS = -g -fbacktrace -fconvert=big-endian
# preprocessor, compilers, linker & archiver
AR  = ar

# objectlist file
include Makefile-objs

# target names
LIB = $(LIBDIR)/libbamsigio_i4r4

# Not real file targets
.PHONY: $(MAKEFILE) Makefile-deps Makefile-objs \
        all

# targets
all: $(LIB)


$(LIB): $(OBJS)
	@$(ECHO) ""
	@$(ECHO) "Creating archive $(@F)"
	@$(ECHO) ""
	$(AR) $(ARFLAGS) $@.a $(notdir $(OBJS))
	$(CP) *.mod $(INC_DIR)

# cleanup
clean:
	$(RM) -f $(LIB) $(notdir $(OBJS)) *.[dlMT] *.lst *.mod work.pc* core

# suffixes
.SUFFIXES: .h .F .f .F90 .f90 .c .o

# remove target on error
.DELETE_ON_ERROR:

# implicit rules
# Want full path? Change $(<F) to $< and add -o $(@F) or -o $@
%.o: %.F   ; $(SF90) -c $(FFLAGS) $(FFLAGS_SIGIOBAM) $(<F)
%.o: %.f   ; $(SF90) -c $(FFLAGS) $(FFLAGS_SIGIOBAM) $(<F)
%.o: %.F90 ; $(SF90) -c $(FFLAGS) $(FFLAGS_SIGIOBAM) $(<F)
%.o: %.f90 ; $(SF90) -c $(FFLAGS) $(FFLAGS_SIGIOBAM) $(<F)

#%.o: %.F90
#	$(CPP) $(CPP_FLAGS) $(CPP_F90FLAGS) $*.F90  > $*.fpp
#	$(FC) $(FFLAGS) -c $*.fpp

# if the compiler do no support the F90 extension
#%.o: %.F90
#	$(MV) $(<F) $(*F).c
#	$(FPP) $(FPPFLAGS) $(*F).c > $(*F)-tmp.f90
#	$(FC) -c $(FFLAGS) $(INCDIRS) $(*F)-tmp.f90
#	$(MV) $(*F)-tmp.o $(*F).o
#	$(RM) -f $(*F).c $(*F)-tmp.f90 

# dependencies file
include Makefile-deps
