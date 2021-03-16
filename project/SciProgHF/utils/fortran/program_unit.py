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
## program_unit.py
##------------------------------------------------------
##
## A class 'program_unit' to represent a fortran program unit
## (that is something like a function or subroutine,
##  modules are not supported yet)
##
## important attributes:
##    - lines              -- a list of 'fortran.line's
##    - declared_variables -- a dictionary with the explicitly
##                            declared variables
##
## important methods:
##    - is_declared(name)  -- tells if 'name' is explicitly declared
##    - add_explicit_declarations -- adds explicit declarations for
##                                   all variables that are used
##
## constructing a new program unit works as follows:
##    - call constructor (optional: ignore_errors)
##    - call append_line (number, str) for every line in this unit
##    - call unit_is_complete ()
##
## error handling:
##     The tokenizing and parsing of the lines can raise exceptions.
##     When the optional 'ignore_errors' argument to the constructor
##     is true, these exceptions are caught and the lines that
##     cause errors are just removed (after printing a worning)
##
##------------------------------------------------------

import tokenizer, statement, error
from variable import *
from line import *

class program_unit :
	""" A class that represents a unit of a Fortran program
	    (that is usally one subroutine)

        Attributes:
	    lines (list of fortran_lines)
	    continuation_pending (boolean) -- Did the last line end with a continuation sign ?
            ignore_errors (boolean) -- If true, ignore all lines where errors appeared

	    declared_variables  -- dictionary of all explicitly declared variables, indexed by name
	    last_decl_statement -- line numer of the last declaration statement

	"""

	def __init__ (self, ignore_errors=False) :
		self.lines = []
		self.continuation_pending = False
                self.ignore_errors = ignore_errors
                
	def unit_is_complete (self) :
		""" Call this when this program unit is complete, so that further setup can be done """
		if len(self.lines) > 0 :
                        try:
                                self.lines[-1].line_is_complete ()
                        except error.FortranException, ex:
                                if self.ignore_errors :
                                        print str(ex)
                                        print "Error is ignored"
                                        del self.lines[-1]
                                else :
                                        raise ex
                                
		## find all explicitly declared variables
		self.declared_variables = {}
		for l in self.lines :
			if l.statement.is_executable :
				break
			else :
				self.last_decl_statement = max(l.numbers)
			## we dont check for mutliple declarations here, the last one just overwrites all previous ones
			self.declared_variables.update( l.statement.variables_declared() )

	def is_declared(self, var, change_dimension=False) :
		""" Check if the variable 'var' is explicitly declared

		    if change_dimension == True, it add a dimension specification to the previous declaration 

		    raises:
                       FortranError -- when dimensions of previous declaration dont match

                """

		f_var = variable(tokenizer.tokenize(var))

		if not self.declared_variables.has_key (f_var.name) :
			return False
		else:
			if f_var.dimension == "" or f_var.dimension == self.declared_variables[f_var.name].dimension :
				return True
			else:
				if self.declared_variables[f_var.name].dimension == "" :
					self.declared_variables[f_var.name].dimension = f_var.dimension
					return True
				else :
                                	msg = "Variable is already declared, but with wrong dimensions\n " + \
                                      	      "  variable name: " + f_var.name + "\n" 

					raise error.FortranError(msg)
	

	## used for debugging only
        def output (self) :
                print "Program Unit : "
                for i in self.lines :
                        print i.line
                        print i.tokens
                print "Explicit Variables: "
                for i in self.explicit_declarations :
                        i.output ()
                print "\n"

	def append_line (self, line_number, new_line) :

		continuation_line = False
		new_line = new_line.rstrip(' \n').replace('\t', '')
		
		## is it a comment or empty line?
		if new_line == "" or new_line[0] in ["C", "c", "*"] or new_line.lstrip(' ')[0] == "!" :
			return

		# is this a continuation line ?
		# F90-style
		if self.continuation_pending :
			if new_line.lstrip (" ")[0] == "&" :
				new_line = new_line.lstrip (" ")[1:]
			continuation_line = True
			self.continuation_pending = False
		# F77 style:
		# we interpret it as continuation only if there is neither a space nor a letter
		# in column 6 to prevent trouble with free form source starting in column 6
		# If anybody used a letter in column 6 as continuation marker, the code
		# is pretty fucked up anyway 
		if new_line[:5].isspace () and not new_line[5].isalpha () and not new_line[5].isspace () :
			new_line = new_line[6:]
			continuation_line = True

		# continuation pending ?
		if len(new_line) > 0 and new_line[-1] == '&' :
			new_line = new_line[:-1]
			self.continuation_pending = True

		if not new_line == "" :
			if continuation_line :  
				self.lines[-1].continue_line (line_number, new_line)
			else:
				if not self.lines == [] :
                                        try:
                                                self.lines[-1].line_is_complete ()
                                        except error.FortranException, ex:
                                                if self.ignore_errors :
                                                        print str(ex)
                                                        print "Error is ignored"
                                                        del self.lines[-1]
                                                else :
                                                        raise ex

				self.lines.append ( line (line_number, new_line) )


	def add_explicit_declarations (self) :
		""" add a edxplicit declaration for all variables, that are only defined implicitly
		    (used in a stement, but not declared) """

		for l in self.lines :
			undeclared = []
			for var in l.statement.variables_used() :
                                try:
                                        if not self.is_declared(var, True) :
                                                undeclared.append(var)
                                except error.FortranError, ex :
                                        ex.set_line_number(l.numbers[-1])
                                        raise ex
					
			if len(undeclared) > 0 :
				statem = statement.new_declaration(undeclared)
				self.declared_variables.update (statem.variables_declared())

				if l.statement.is_declaration :
					new_line = inserted_line (min(l.numbers), statem, False)
				else :
					new_line = inserted_line (last_decl_statement, statem, True)

				self.lines.append (new_line)

			if l.statement.type == "common block" :
				l.statement.remove_dimensions = True
				
