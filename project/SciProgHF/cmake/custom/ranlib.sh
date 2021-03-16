# ranlib.sh; hjaaj January 2017 
# Script for darwin (MacOS) to replace "ranlib" in cmake with
#   set(CMAKE_RANLIB     <path_to_this_directory>/ranlib.sh)
# Note that C language probably needs to be enabled before Fortran for this to work.
# The ar.sh script should also be anabled for this to work.
# This silences the annoying and often many messages as the following example:
#/opt/local/bin/ranlib: file: lib/libdalton.a(par_cfg.F90.o) has no symbols

#echo "Input options to ranlib.sh: " $*

ranlib -no_warning_for_no_symbols $*
