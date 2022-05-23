# Macro file for use with GNU AutoTools for generating config files and Makefiles.
# This file was manually created based on information found via Google and the
# GNU documentation. Although it works for the source code it's included with,
# it may not be 100% correct so additional research/tweaking may be needed.
#
# bin_PROGRAMS is used to define the binary that's to be built from the
# source code. The binary name will also be used in variable names to
# pass values that will only apply to this binary.
bin_PROGRAMS = inctime
# Note: that the following names all start with the binary name defined in
# bin_PROGRAMS.
# <binary_name>_LDADD is used to pass extra parameters at link time, such as
# libraries that need to be linked in.
inctime_LDADD = -lmisc
# <binary_name>_LDFLAGS is used to pass extra parameters at compilation time,
# such as the paths to libraries that are needed that are not in the default
# paths.
# If the user supplies a directory where libpcap is installed, the contents
# of the variable setup in configire.in will be substituted below so the
# library file will be found
inctime_LDFLAGS = -L../libmisc
# <binary_name>_SOURCES takes a space delimited list of the source code files
# that need to be compiled/linked to build the binary.
inctime_SOURCES = inctime.f90
