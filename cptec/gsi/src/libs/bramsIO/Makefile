# Automatically generated GNU Makefile. Dom Jul 27 11:51:57 BRT 2014
# fgen v0.3 (C) 1997,98 Beroud Jean-Marc

include ../../../configure.gsi

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

# preprocessor, compilers, linker & archiver
AR  = ar

# objectlist file
include Makefile-objs

# target names
LIB = $(LIBDIR)/libbramsio_i4r4

# Extra Flags
FFLAGS_BRAMSIO += #$(FFLAGS_DEBUG)

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
%.o: %.F   ; $(SF90) -c $(FFLAGS_BRAMSIO)  $(<F)
%.o: %.f   ; $(SF90) -c $(FFLAGS_BRAMSIO)  $(<F)
#%.o: %.F90 ; $(SF90) -c $(FFLAGS_BRAMSIO)  $(<F)
%.o: %.f90 ; $(SF90) -c $(FFLAGS_BRAMSIO)  $(<F)
%.o: %.c   ; $(SCC)  -c $(CFLAGS2) $(<F)
.F90.o:
	$(CPP) $(CPP_FLAGS) $(CPP_F90FLAGS) $*.F90  > $*.f90
	$(F90) $(FFLAGS) -c $*.f90
	$(RM) $*.f90

# dependencies file
include Makefile-deps
