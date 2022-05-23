#
#  Author: pkubota $
#  Date: 2007/10/10 20:28:03 $
#  Revision: 1.4 $
#  02/09/2021 - Eduardo Khamis - v2.2.1
#
SHELL=/bin/bash
#
# makefile for the global model
#
# Targets are the suffixes of the system-specific makefiles in
# the makefiles/ directory.
# For example, to build PosGrib for Solaris, give the command
#
#   make solaris
#
### then builds executables in the exec/ directory.
#
# This builds an intermediate library in the util/ directory,
# then builds the nedit and nc executables in the source/ directory.
#
#

MAKE        = make

LV          = 28
#LV          = 64

BUILD_DIR   = $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

MIN_DIR     = `cd ${BUILD_DIR};cd ../;pwd`

SOURCE_DIR  = $(MIN_DIR)/source

LIB_DIR     = $(SOURCE_DIR)/Physics/Radiation/RRTMG/ext/rrtmg_sw

LIB2_DIR    = $(SOURCE_DIR)/Physics/Radiation/RRTMG/ext/rrtmg_lw



all:
	@echo "Please specify target:"
	@echo "(For example, type \"make linux\" for a Linux system.)"
	@(cd  makefiles && ls -C makefile* | sed -e 's/makefile.//g')


.DEFAULT:


	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f ${LIB_DIR}/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         ${LIB_DIR}/makefile.$@; fi)

	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f ${LIB2_DIR}/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         ${LIB2_DIR}/makefile.$@; fi)

	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f ${BUILD_DIR}/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         ${BUILD_DIR}/makefile.$@; fi)

#	@- (cd .; cat ${SOURCE_DIR}/Assimilation/SpecDump.f90.in | sed s/"NLEVS=28"/"NLEVS=$(LV)"/g \
		    > ${SOURCE_DIR}/Assimilation/SpecDump.f90) 


#
# builds the binaries
#

	@- (cd $(LIB_DIR)              ;$(MAKE) -f makefile.$@)
	(cd .;cp -fp ${LIB_DIR}/*.o .)
	(cd .;echo ll.mod > ${LIB_DIR}/ll.mod ;cp -fp ${LIB_DIR}/*.mod .)

	@- (cd $(LIB2_DIR)             ;$(MAKE) -f makefile.$@)
	(cd .;cp -fp ${LIB2_DIR}/*.o .)
	(cd .;echo ll.mod > ${LIB2_DIR}/ll.mod ;cp -fp ${LIB2_DIR}/*.mod .)

	@- (cd $(BUILD_DIR)            ;$(MAKE) -f makefile.$@)


# We need a "dev-all" target that builds the docs plus binaries, but
# that doesn't work since we require the user to specify the target.  More
# thought is needed

clean:
	@echo "$(MAKE)"

	@- (cd $(LIB_DIR)              ;$(MAKE) -f Makefile.common clean)
	@- (cd $(LIB2_DIR)             ;$(MAKE) -f Makefile.common clean)
	@- (cd $(BUILD_DIR)            ;$(MAKE) -f Makefile.common clean)

install:
	@echo "$(MAKE)"

	@- (cd $(BUILD_DIR);$(MAKE) -f Makefile.common install)

