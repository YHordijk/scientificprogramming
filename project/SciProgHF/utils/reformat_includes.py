#! /usr/bin/env python
##  -*- Mode: Python; py-indent-offset: 8 -*-

##------------------------------------------------------
## reformat_includes.py
##------------------------------------------------------
##
## Written by
## Christoph Jacob <jacob@few.vu.nl>, May 2004
##
##------------------------------------------------------

## * Comments:
## Reformats F77 fixed form source comments (starting with C or * in first column)
## to F90 style comments starting with !
##
## * Continuation
## Reformats F77 fixed form source continuation so that the code can be used in
## both F77 and F90 code.
## This is done by adding a & in column 73 of the line that is continued.
## This & is ignored by the F77 compiler but makes the F90 compiler happy
## (which in turn ignores the & in column 6 of the next line)
##
## We interpret it as continuation only if the is neither a space nor a letter
## in column 6 to prevent trouble with free form source starting in column 6
## If anybody used a letter in column 6 as continuation marker, the code
## is pretty fucked up anyway
##
## * Declarations
## Adds explicit declarations of all used variables, if they are not already
## there

import sys, os, fortran

major, minor, dum, dum , dum = sys.version_info
if not (major >= 2 and minor >= 3) :
	print "you need Python version >= 2.3 to run this script "
	sys.exit(1)

def process_file (filename, options) :

	try:
		f = fortran.ffile (filename, options.ignore_errors)
	except IOError:
		print "\nfile "+filename+" not found"
		return 1
	except fortran.error.FortranException, ex :
		print "Error processing file"+filename+":"
		print str(ex)
		return 1
			
	try:
		f.reformat()

		if options.backup :
			try:
				f.backup(options.backup_suffix, options.overwrite_backup)
				print "saved backup copy of original version"
			except fortran.error.FortranIOError, ex :
				print str(ex)
				print "no backup was saved"

		f.save()

	except fortran.error.FortranException, ex :
		print "Error processing file"+filename+":"
		print str(ex)
		return 1	

	return 0


## ---- MAIN PROGRAM ----

## parse command line options
from optparse import OptionParser

version="%prog, Version 1.0 \nA tool to reformat Fortran 77 include files so that they can be used in F90 programs \n"+\
        "Written by Christoph Jacob <jacob@few.vu.nl>, May 2004\n"

optionparser = OptionParser(usage="%prog [options] filenames\nUse option --help for details", version=version)

optionparser.add_option("-f", "--force", action="store_true", dest="ignore_errors", default=False,
			help="Try to ignore errors if possible")
optionparser.add_option("--no-backup", action="store_false", dest="backup", default=True,
			help="Do not save backups of the original file")
optionparser.add_option("-s", "--suffix", action="store", dest="backup_suffix", type="string", default="orig",
			metavar="SUFFIX", help="suffix for backup copies [default: orig]")
optionparser.add_option("-o", "--overwrite-old-backups", action="store_true", dest="overwrite_backup", default=False,
			help="Overwrite existing backup files")
optionparser.add_option("-e", "--skip-if-backup-exists", action="store_true", dest="skip", default="False",
                        help="Skip all files where a backup already exists ")


(options, args) = optionparser.parse_args()

if len(args) == 0 :
	optionparser.error("you have to give the names of the files to reformat")

## process the files

failed = []

for filename in args :
	if options.backup_suffix[0] == '.' :
		backup_filename = filename + options.backup_suffix
	else :
		backup_filename = filename + '.' + options.backup_suffix
	if not options.skip or not os.path.exists(backup_filename) :

		print "\nProcessing file "+filename

		if not process_file (filename, options)	== 0 :
			failed.append(filename)
			print "this file is left unchanged. Try using --force option to ignore errors "
		else :
			print "Successfully saved reformatted version of file " + filename

if not failed == [] :
	print "\nErrors occured for the following files: "
	for n in failed :
		print n
	print ''
	sys.exit(1)
