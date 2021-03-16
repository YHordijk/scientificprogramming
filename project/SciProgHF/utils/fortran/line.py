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
## line.py
##------------------------------------------------------
##
## A class "line" to represent a line in a Fortran file
##
## The most important attributes are:
##    - str:       The line as a string
##    - tokens:    The tokenized line
##    - numbers:   The line numbers in the corresponsing file
##    - statement: The corresponsing statement
##
## Constructing a new line works as follows:
##    - call the constructor line(num, str)
##    - call line.continue_line(num, str) to append continuation lines
##    - call line.line_is_complete() :
##             * This tokenizes and parses the line
##
## Other important method:
##    - str: gives the line as a string
##
## There is also the class inserted_line, that is intended
## for lines that are newly inserted into the file.
## It only has a modified constructor.
##
## It is intended to be used only from file.py and program_unit.py, not by itself
##
## error handling:
##     possible exceptions from the parser or tokenizer are caught
##     and the line number is set, than it is passed on
##
##------------------------------------------------------

import tokenizer, statement, error

class line :
	""" A simple class representig a line (or somthing quite similar) in a fortran_file

	Attributes:
	    string (String)               -- The line itself
            tokens                     -- The tokenized line
	    numbers (list of Integers) -- The number(s) of the corresponding line(s) in the fortran_file

	    insert (bool)              -- is this a new line, that has not been in the fortran_file?
	    insert_at (Integer)        -- if True, insert after/before this line 
            insert_after (bool)        -- if True insert after line 'insert_at',
	                                  if False isert before line 'insert_at'  		

	    statement (derived from statement.unknown)
	                               -- The statement of this line
				          (availible after setup_statement was called)
	"""
        
	def __init__(self, num, str) :
		self.string = str
		self.numbers = [ num ]
		self.insert = False

        def continue_line (self, num, str) :
                self.string = self.string + str
                self.numbers.append (num)

	def line_is_complete (self) :
		""" Should be called we the line is complete, so we can tokenize and analyze it """
                
                self.tokens = tokenizer.tokenize(self.string)
                try:
                        self.setup_statement ()
                except error.FortranException, ex :
                        ex.set_line_number(self.numbers[0])
                        raise ex

	def setup_statement (self) :
		if isinstance(self.tokens[-1], str) and self.tokens[-1].startswith("!") :
			temp_tokens = self.tokens[:-1]
		else :
			temp_tokens = self.tokens

		self.statement = statement.new_statement (temp_tokens)

	def str (self) :
		## This is implemented for F77 fixed form at the moment!

                try:
                        lines = self.statement.str(66)
                except error.FortranException, ex :
                        ex.set_line_number(self.numbers[0])
                        raise ex

		result = ""
		for l in lines :
			result += (6*" " + l[0])
			for i in range(1, len(l)) :
				if not (i == len(l)-1 and l[i] == "") :
					## add continuation mark to previous line
					j = result.rfind('\n')+1
					result += (72 - len(result[j:]))*" " + "&\n"
					## the new line
					result += 5*" "+ "&" + l[i]
				else :
					result += "\n"
			result += "\n"
			
		return result

	## used for debugging only
	def output (self) :
		print "numbers:      ", self.numbers
		print "string:          ", self.string
		print "tokens:       ", self.tokens
		print "statement:"
		print "  type:       ", self.statement.type
		print "  str():     "
		##print "0        1         2         3         4         5         6         7   "
		##print "1234567890123456789012345678901234567890123456789012345678901234567890123"
		print self.statement.str(200)
		print "  variables_used: ", self.statement.variables_used()

		if self.insert :
			print "inserted line"
			print "insert_at", self.insert_at
			print "insert_after", self.insert_after


class inserted_line (line) :

	def __init__ (self, insert_at, statement, insert_after=False) :
		line.__init__(self, 0, "")
		self.numbers = []
		self.tokens = []
		self.statement = statement

		self.insert = True
		self.insert_at = insert_at
		self.insert_after = insert_after
