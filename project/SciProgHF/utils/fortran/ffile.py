##  -*- Mode: Python; py-indent-offset: 8 -*-

##------------------------------------------------------
##
## Python Package "fortran"
## A library for processing Fortran files 
##
## Written by
## Christoph Jacob <jacob@few.vu.nl>, May 2004
##
##------------------------------------------------------
## ffile.py
##------------------------------------------------------
##
## A class representing a Fortran file
##
## important attributes:
##   - filename
##   - lines         -- The lines of the file (as strings)
##   - program_units -- A list of fortran.program_units
##
## important methods:
##   - constructor fortran.file(filename)
##   - backup()   -- make a backup copy
##   - save()     -- save the (possibly new version) under the origninal name
##   - reformat() -- 1. reformats comments
##                   2. add explicit declarations for all used variables 
##
## error handling:
##   the following exceptions can appear:
##
##   IOError -- when the constructor or save have trouble reading or writing
##              when backup has problems copying
##
##   FortranIOError -- raised in backup when filename.orig already exists
##
##   various FortranExceptions  -- when the constructor, reformat and so on
##                                 encounter problems
##                                 if the constructor is called with the optional argument
##                                 ignore_errors set to true, a lot of them will be ignored
##                                 and only produce warnings
##
##------------------------------------------------------

import os, shutil
import error
from program_unit import *		
                
                                

class ffile :
        """ A class to represent a fortan file and to do some processing
            
        Attributes:
            filename (String)          -- The filename
	    lines (list of Strings)    -- The lines as read from the file (read-only)

	    program_units (list of fortran.program_units)
	                               -- The programm units 

        """

## Constructor
        def __init__ (self, filename, ignore_errors=False, do_prog_units=True):
		self.filename = filename
		
                f = file (filename, 'r')
                self.lines = f.readlines () 
                f.close

		if do_prog_units :
			self.setup_program_units (ignore_errors)

## Make a copy with suffix appended to filename
## Raises FortranIOError if backupfile already exists

        def backup (self, suffix, force=False):
		if not suffix[0] == '.' :
			suffix = '.' + suffix 
                backup_filename = self.filename + suffix
                if not force and os.path.exists (backup_filename):
                        msg = "Backup file "+ backup_filename + " exists"
                        raise error.FortranIOError(msg) 
                else:
			shutil.copyfile (self.filename, backup_filename)

## Print the file to stdout
 
        def str (self):
		
		lines = self.lines
		insert_before = {}
		insert_after = {}
		for unit in self.program_units :
			for f_line in unit.lines :

				if f_line.insert :
					if f_line.insert_after :
						if f_line.insert_at in insert_after :
							insert_after[f_line.insert_at] += f_line.str()
						else :
							insert_after[f_line.insert_at] = f_line.str()
					else :
						if f_line.insert_at in insert_before :
							insert_before[f_line.insert_at] += f_line.str()
						else :
							insert_before[f_line.insert_at] = f_line.str()
				else :
					for i in f_line.numbers :
						lines[i] = ""
					lines[f_line.numbers[0]] = f_line.str()

				## insert empty lines
				if f_line.statement.is_declaration and not f_line.statement.type == "declaration" :
					if f_line.insert and f_line.insert_after:
						insert_after[f_line.insert_at] += "\n"
					elif f_line.insert :
						insert_before[f_line.insert_at] += "\n"
					else :
						next = f_line.numbers[-1] + 1
						if next < len(lines) and not lines[next].strip(' ') == "\n" :
							lines[f_line.numbers[0]] += '\n'

		result = ""
		for i in range(len(lines)) :
			if i in insert_before :
				result = result + insert_before[i]
			result = result + lines[i]
			if i in insert_after :
				result = result + insert_after[i]

		return result
		
## Save the file, overwriting the original version
## use the parsed lines

        def save (self):
                f = file (self.filename, 'w')
                f.write ( self.str() )
                f.close

## Save the file, overwriting the original version
## use the lines as strings

        def save_strings (self):
                f = file (self.filename, 'w')
               	f.writelines ( self.lines )
                f.close

## Get the Fortran part (for files that have a Fortan and a C part)

        def fortran_part (self):
                result = []; in_C_part = False; in_F_part = False; nested = 0;
                for i in range (len(self.lines)) :
                        stripped = self.lines[i].replace(" ", "").replace("(", "");

                        # We are in the C part or Fortran part: handle nested "#if"'s
                        if (in_C_part or in_F_part) and stripped.startswith ("#if") :
                                nested = nested + 1
                        elif (in_C_part or in_F_part) and nested > 0 and stripped.startswith ("#endif") :
                                nested = nested - 1

                        # Beginning of C part by #ifdef __CVERSION__ or something similar
                        if (stripped.startswith ("#ifdefined") and stripped[10:].startswith ("__CVERSION__") )        \
                               or (stripped.startswith ("#ifdef") and stripped[6:].startswith ("__CVERSION__") ) :
                                in_C_part = True;

                        # Beginning of explicit Fortran part 
                        if (stripped.startswith ("#if!defined") and stripped[11:].startswith ("__CVERSION__") )       \
                               or (stripped.startswith ("#ifndef") and stripped[7:].startswith ("__CVERSION__") ) :
                                in_F_part = True;


                        # Handle else branch
                        if nested == 0 and stripped.startswith ("#else") :
                                if in_C_part:
                                        in_C_part = False
                                elif in_F_part:
                                        in_F_part = False
                                        in_C_part = True

                        # if we are not in the C part and this line is no preprocessor statement,
                        # we want to use it!
                        if not in_C_part and not stripped.startswith("#") and not stripped.rstrip('\n') == "":
                                result.append (i)
                        
                return result;
 
	def setup_program_units (self, ignore_errors=False) :

		self.program_units = [ program_unit (ignore_errors) ]
		
		for i in self.fortran_part () :
			# new program units begin with FUNCTION or SUBROUTINE
			# (at least in F77 code)
			#
			# FIXME: add support for F90 features (MODULES etc. )
			if self.lines[i].lstrip(" ").upper().startswith("FUNCTION") or \
			   self.lines[i].lstrip(" ").upper().startswith("SUBROUTINE") :
				if not self.program_units[-1].lines == [] :
                                        self.program_units[-1].find_declarations ()
					self.program_units.append( program_unit(ignore_errors) )
			self.program_units[-1].append_line (i, self.lines[i])

                self.program_units[-1].unit_is_complete ()

	def reformat_comments (self) :
		for i in self.fortran_part () :
			if not self.lines[i].strip(' \n') == "" :
				if self.lines[i][0] in ['C', 'c', '*'] :
					self.lines[i] = '!' + self.lines[i][1:]
				if self.lines[i].strip(' \n') == "!" and i>0 and \
				       not self.lines[i-1].strip(' \n')[0] == '!' :
					self.lines[i] = "\n"   

	def reformat (self) :
		self.reformat_comments()

		for i in range(len(self.program_units)) :
			self.program_units[i].add_explicit_declarations ()

	def add_copyright (self, dirac_copyright) :

		cr_string = ""
		for l in dirac_copyright :
			cr_string = cr_string + l
		self.lines[0] = cr_string + self.lines[0]


	def remove_copyright (self) :
## loop ueber Zeilen, Ende suchen
		end_of_header = 0
		for i in range(len(self.lines)) :
			if not self.lines[i].strip(' \n') == "" :
				if not self.lines[i][0] in ['C', 'c', '*', '!', '#'] :
					## non-comment line -> we have not found it
					break
				if self.lines[i][0] in ['C', 'c', '*', '!'] :
					if not self.lines[i].find("http://dirac.chem.ou.dk") == -1 :
						end_of_header = i+1
						break
                                        if not self.lines[i].find("http://dirac.chem.sdu.dk") == -1 :
                                                end_of_header = i+1
                                                break
                                        if not self.lines[i].find("http://www.kjemi.uio.no") == -1 :
                                                end_of_header = i+1
                                                break

## loeschen, nur Kommentare
		if end_of_header > 0 :
			for i in range(end_of_header) :
				if not self.lines[i].strip(' \n') == "" :
                               		if self.lines[i][0] in ['C', 'c', '*', '!'] :
                                        	self.lines[i] = ""
		






