SHELL=/bin/sh
include ../../../configure.gsi
LIBDIR = ../../../lib/
LIB = $(LIBDIR)/libw3emc_i4r8.a

# The file definitions. This include must occur before targets.
include make.filelist

all:    $(LIB)

$(LIB): $(OBJS)
	$(AR) -ruv $(LIB) $(OBJS)
	$(CP) *.mod $(INC_DIR)

#gblevents.o:
#	$(SFC) -c $(FFLAGS_W3) gblevents.f


# Dependency include file
include make.dependencies

.SUFFIXES :
.SUFFIXES : .f .a .c .o

.f.o:
	$(SFC) -c $(FFLAGS_W3) $<

.c.o:
	$(SCC) -c $(CFLAGS) $<

clean:
	$(RM) $(OBJS) *.mod $(INC_DIR)/args_mod.mod  $(INC_DIR)/gblevn_module.mod  $(INC_DIR)/mersenne_twister.mod $(LIB)

