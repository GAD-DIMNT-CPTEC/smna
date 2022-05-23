# The lib_LIBRARIES macro tells the autotools the name and type of library we
# want to build - in this case _not_ a shared library. A shared library
# would end with '.la' instead of '.a'.
lib_LIBRARIES = libmisc.a
# This specially named variable contains a space delimited list of the source code
# files that must be compiled in order to create the library. The variable name
# is in the format:
# <library_name>_<library_type>_SOURCES
# Where:
# <library_name> = the name of the library specified in lib_LIBRARIES
# <library_type> = either 'a' for non-shared library or 'la' for shared.
libmisc_a_SOURCES =  m_inpak90.F90  m_stdio.f90  m_string.f90  m_time.f90
# Header files that must be installed are specified by the HEADERS family
# of variables.
include_HEADERS = ../include/m_inpak90.mod ../include/m_stdio.mod ../include/m_string.mod ../include/m_time.mod
