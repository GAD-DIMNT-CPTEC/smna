SHELL=/bin/sh
include ../../../configure.gsi
LIBDIR = ../../../lib/
LIB = $(LIBDIR)/libsp_i4r8.a
INCMOD  = ../../../include


OBJ = bll2ps.o spanaly.o spdz2uv.o speps.o spfft.o spfft1.o \
spffte.o spfftpt.o spgradq.o spgradx.o spgrady.o splaplac.o \
splegend.o sppad.o spsynth.o sptez.o sptezd.o sptezm.o sptezmd.o sptezmv.o \
sptezv.o sptgpm.o sptgpmd.o sptgpmv.o sptgps.o sptgpsd.o sptgpsv.o sptgpt.o \
sptgptd.o sptgptsd.o sptgptv.o sptgptvd.o sptran.o sptrand.o sptranf.o sptranf0.o \
sptranf1.o sptranfv.o sptranv.o sptrun.o sptrund.o sptrung.o sptrungv.o sptrunl.o \
sptrunm.o sptrunmv.o sptruns.o sptrunsv.o sptrunv.o spuv2dz.o spvar.o spwget.o \
fftpack.o ncpus.o splat.o lapack_gen.o

all:  $(LIB)

$(LIB): $(OBJ)
	$(AR) -ruv $(ARFLAGS_SP) $(LIB) $(OBJ)

.SUFFIXES:
.SUFFIXES: .F .f .o

.F.o:
	$(CPP) $(CPP_FLAGS) $(CPP_F90FLAGS) $*.F  > $*.f
	$(SFC) -c $(FFLAGS_SP) $*.f
	$(RM) $*.f

.f.o: 
	$(SFC) $(FFLAGS_SP) -c $*.f

clean:
	$(RM) $(OBJ) $(LIB)

