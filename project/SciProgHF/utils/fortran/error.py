##  -*- Mode: Python; py-indent-offset: 8 -*-
##
##------------------------------------------------------
##
## Python Package "fortran"
## A library for processing Fortran files 
##
## Written by
## Christoph Jacob <jacob@few.vu.nl>, May 2004
##
##------------------------------------------------------
## error.py
##------------------------------------------------------
##
## The Exceptions that are used in this package
##
## The base class is FortranException, all
## exceptions that are raised somewhere are derived
## from it
##
##------------------------------------------------------

class Bug (Exception) :
	""" Exception that is raised, when something obviously went wrong
	    It is NOT derived from FortranException, so that is it not catched
	    by normal error handling code
	"""
	pass

class FortranException (Exception) :
        """ Base class for exceptions used in this package
	Attributes:
	    message :      An appropriate error messgage
	    line_number :  the number of the line connected to the error or 0 if unknown

	"""

        def __init__ (self, message, num=0) :
		self.line_number = num
		self.msg = message

	def __str__ (self) :
		message = self.msg
		if not self.line_number == 0 :
			message = self.msg + "\n" + "The error occured on line: " + str(self.line_number) 
		return message

	def set_line_number (self, num) :
		self.line_number = num
        
class NotImplementedError (FortranException) :
	""" Exception that is raised, when a feature is requested, that is not implemented yet

       """

	def __init__ (self, message, num=0) :
		msg = "NotImplementedError:\nA feature was requested that is not implemented (yet) \n" + message
		FortranException.__init__(self, msg, num)

class TokenizerError (FortranException) :
	""" Exception used in the Tokenizer """

	def __init__ (self, message) :
		msg = "TokenizerError:\n" + message
		FortranException.__init__(self, msg, 0)

class ParsingError (FortranException) :
	""" Exception used in statement.py, when something cant be interpreted correctly
	    (usually it's something like a 'syntax error')
	"""

	def __init__ (self, message) :
		msg = "ParserError:\n" + message
		FortranException.__init__(self, msg, 0)

class FortranError (FortranException) :
	""" Exception used when something seems to be wrong with the original Fortran file
	    and its not just the tokenizer or the parse 
	"""

	def __init__ (self, message, num=0) :
		msg = "FortranError:\n" + message
		FortranException.__init__(self, msg, num)

class FortranIOError (FortranException) :
	""" Exception raised when IO-related things go wrong """

	def __init__ (self, message) :
		msg = "FortranIOError:\n" + message
		FortranException.__init__(self, msg, 0)
