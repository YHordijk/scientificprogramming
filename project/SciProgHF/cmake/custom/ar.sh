# ar.sh; hjaaj January 2017 
# Script for darwin (MacOS) to replace "ar" in cmake with
#   set(CMAKE_AR     <path_to_this_directory>/ar.sh)
# Note that C language probably needs to be enabled before Fortran for this to work.
# The ranlib.sh script should also be anabled for this to work.
# This silences the annoying and often many messages as the following example:
#/opt/local/bin/ranlib: file: lib/libdalton.a(par_cfg.F90.o) has no symbols

#echo "Input options to ar.sh: " $*

# note: CMAKE parameter list starts with "qc" as in "qc somelib.a somthing" (i.e. not "-qc" or "-q -c")
#       and "-S qc" means qc is interpreted as library name. Thus no space between -S and $* !!!!!!
#       There is no way in cmake to change the argument list to ar, according to the internet.
#       Thus this "trick" of replacing "ar" with this shell script..
#       The "-S" prevents ar from running an internal ranlib each time a file is added to the library.

ar -S$*
