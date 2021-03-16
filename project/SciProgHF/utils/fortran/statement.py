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
## statement.py
##------------------------------------------------------
##
## This is something like a parser:
## It takes a list of tokens and produces objects corresponding to these statements.
## The classes for the different statement types provide mainly:
##  - a constructor, that also does the parsing
##  - several attributes (and maybe methods) that provide information on the statement
##  - the line and string methods, that convert the statement to a string
## (for details, see comments in the implementation)
##
## new_statement(tokens) returns a new statement object
##
## availible statement classes are:
##   - unknown     -- for everything we dont recognize, base class for all other classes
##
##   - declaration -- for excplicit variable declarations
##   - parameter
##   - equivalence
##
##   - executable  -- base class for executable statements (not used yet)
##
## new_declaration (list of variable names) produces a new declaration of the given variables
##
##
## error handling :
##    - The constructors can raise "ParsingError"s
##    - The str and method uses tokenizer.join_tokens, which can raise "TokenizerError"s
##    
##------------------------------------------------------

import tokenizer, error
from variable import *

class unknown :
	""" Base class for fortran statements
	    Also used for statements, that are not understood (yet)

	Attributes:   !! Take care they exist in all derived classes !!
            type              -- A string identifying the type of the statement 
	    is_declaration    -- is this a (explicit) declaration statement ?
	                         (variable/parameter declarations, common blocks, etc.)
            is_executable     -- everything else

	    tokens            -- the corresponding list of tokens
	"""

	def __init__ (self, tokens) :
		self.is_declaration = False
		self.is_executable = False
		self.tokens = tokens
		self.type = "unknown"

	def str (self, max_len_first, max_len_following=0) :
		""" Gives this line as a of lists of strings. Each list is one 'program line'
		    (that means: actually belongs to one line, but is maybe split over
		    multiple line that have to be joined by appropriate continuation marks).
		    The first line with maximum length max_len_first,
                    the following lines of maximum length max_len_following
		    (The default values max_len_following=0 means 'same as first line')
		    For comments, the maximum lengths can be ignored! 
		"""

		## This version indents continuation lines by 1 space
		if not max_len_following == 0 :
			max_len_following -= 1
		
		result = tokenizer.join_tokens(self.tokens, max_len_first, max_len_following)
		result[1:] = [ (" " + r) for r in result[1:] ]

		return [result]
		
	def variables_used (self) :
		""" return a list of the names of the variables in this statement """
		return []

	def variables_declared (self) :
		""" return a dictionary of the variables declared in this statement """
		return {}
		
class declaration (unknown) :
	""" A statements explicitly declaring a variable

        Attributes (in addition to unknown):
	    variable type --
	    variables     --
	    
	"""

	def __init__ (self, tokens) :
		unknown.__init__ (self, tokens)
		self.is_declaration = True
		self.variables = []
		self.type = "declaration"

		## we want to buid a list like that: [['TYPE'], ['name1', 'name2', ...] ]

		## F90 syntax is easy:
		if "::" in self.tokens :
			declaration = tokenizer.split_token_list (self.tokens, ["::"])
			declaration[1] = tokenizer.split_token_list (declaration[1], [','])

		## F77 is harder:
		else :
			end = 1;
			if self.tokens[end] == '*' :
				if self.tokens[end+1] == '(' :
					end += 3
				else :
					end += 2 
			if self.tokens[end] == '(' :
				end += 2
			declaration = [[], []]
			declaration[0] = self.tokens[:end]
			declaration[1] = tokenizer.split_token_list (self.tokens[end:], [','])

		self.variable_type = tokenizer.join_tokens(declaration[0])

		## now process this list;
		for i in range(len(declaration[1])) :
			var = variable(declaration[1][i], declaration[0])
			self.variables.append (var)

	def output (self) :

		print "declaration :    "
		print "  variable_type: ", self.variable_type
		for i in self.variables :
			i.output()

	def str (self, max_len_first, max_len_following=0) :
		""" Gives this line as a list of strings. The first line with maximum length max_len_first,
                    the following lines of maximum length max_len_following
		    (The default values max_len_following=0 means 'same as first line')
		    For comments, the maximum lengths can be ignored! 
		"""

		## F77 style of variable declarations

		if self.variable_type in ["", "CHARACTER"] :
			types = []
			for var in self.variables :
				if var.type == "CHARACTER" :
					t = "CHARACTER*" + str(var.char_len) 
				else :
					t = var.type
				if not t in types :
					types.append (t)
		else :
			types = [self.variable_type]

		result = [ ]

		for t in types :
			result.append( [t + " "] )
			indent = len(result[-1][-1])
			mlen_first = max_len_first - indent                       ## we assume it fits on the line
			if not max_len_following == 0 :
				mlen_following = max_len_following - indent
			else :
				mlen_following = max_len_following

		## for the variable name build a token_list like that:
		## ["INT1", ",", "INT2(10)", ",", ..."
			variable_names = []
			for var in self.variables :
				if t.startswith("CHARACTER") and var.type == "CHARACTER" :
					if t[t.find('*')+1:] == var.char_len :
						variable_names += [var.name + var.dimension, ","] 
				elif var.type == t :
					variable_names += [var.name + var.dimension, ","]

			del variable_names[-1]

			vars = tokenizer.join_tokens(variable_names, mlen_first, mlen_following)

			result[-1][-1] += vars[0]
			for line in vars[1:] :
				result[-1].append (indent*" " + line) 

		return result

	def variables_declared (self) :
		""" return a dictionary of the variables declared in this statement """
		result = {}

		for var in self.variables :
			result[var.name.upper()] = var
		
		return result


class common_block (unknown) :
	""" A statements declaring a common

        Attributes (in addition to unknown):
	    common_name     --
	    variable_names  --

	    remove_dimensions --
	    
	"""

	def __init__ (self, tokens) :
		unknown.__init__ (self, tokens)
		self.is_declaration = True
		self.type = "common block"
		self.remove_dimensions = False

		if not tokens[0].upper() == "COMMON" and tokens[1] == "/" and tokens[3] == "/" :
			raise error.ParsingError("Invalid COMMON statement")
		if "/" in tokens[4:] :
			raise error.NotImplementedError("Two common blocks declared in one COMMON statement", False)

		self.common_name = tokens[2]

		name_list = tokenizer.split_token_list(tokens[4:], [","])
		self.variable_names = [ tokenizer.join_tokens(i) for i in name_list ]
		

	def str (self, max_len_first, max_len_following=0) :
		""" Gives this line as a list of strings. The first line with maximum length max_len_first,
                    the following lines of maximum length max_len_following
		    (The default values max_len_following=0 means 'same as first line')
		    For comments, the maximum lengths can be ignored! 
		"""

		result = ["COMMON /" + self.common_name + "/ "]
		indent = len(result[0])
		max_len_first -= indent                       ## we assume it fits on the line
		if not max_len_following == 0 :
			max_len_following -= indent

		## for the variable names build a token_list like that:
		## ["INT1", ",", "INT2", ",", ... ]
		var_names = []
		for var in self.variable_names :
			if not self.remove_dimensions :
				var_names += [var, ","]
			else :
				var_names += [tokenizer.tokenize(var)[0], ","]
		del var_names[-1]
		
		vars = tokenizer.join_tokens(var_names, max_len_first, max_len_following)

		result[0] += vars[0]
		for line in vars[1:] :
			result.append (indent*" " + line)

		return [result]

	def variables_used (self) :
		""" return a list of the names of the variables in this statement """
		## These names possibly contain dimension specification!
		return self.variable_names



class parameter (unknown) :
	""" A statement declaring a parameter

        Attributes (in addition to unknown):
	    parameters  -- list of (name, value)-pairs
	    
	"""

	def __init__ (self, tokens) :
		unknown.__init__ (self, tokens)
		self.is_declaration = True
		self.type = "parameter"

		if not tokens[0].upper() == "PARAMETER" and tokens[1] == "(" :
			raise error.ParsingError("Invalid PARAMETER statement")

		name_list = tokenizer.split_token_list(tokens[2], [","])
		self.parameters = [ (i[0], tokenizer.join_tokens(i[2:])) for i in name_list ]
		

	def str (self, max_len_first, max_len_following=0) :
		""" Gives this line as a list of strings. The first line with maximum length max_len_first,
                    the following lines of maximum length max_len_following
		    (The default values max_len_following=0 means 'same as first line')
		    For comments, the maximum lengths can be ignored! 
		"""

		result = ["PARAMETER ("]
		indent = len(result[0])
		max_len_first -= indent+1                       ## we assume it fits on the line
		if not max_len_following == 0 :
			max_len_following -= indent

		names = []
		for param in self.parameters :
			names += [param[0],"=", param[1], ","]
		del names[-1]

		## FIXME: maybe there is a really long right hand side in the parameter
		## statement. So catch the exeption and in nessacasy split the rh-sides
		
		params = tokenizer.join_tokens(names, max_len_first, max_len_following)

		result[0] += params[0]
		for line in params[1:] :
			result.append (indent*" " + line)

		result[-1] += ")"
		
		return [result]

	def variables_used (self) :
		""" return a list of the names of the variables in this statement """
		return [i[0] for i in self.parameters]

class equivalence (unknown) :
	""" A EQUIVALENCE statement declaring a parameter

        Attributes (in addition to unknown):
	    equ_lists  -- list of lists, each containing the variables in one EQUIVALENCE list
	    
	"""

	def __init__ (self, tokens) :
		unknown.__init__ (self, tokens)
		self.is_declaration = True
		self.type = "equivalence"

		if not tokens[0].upper() == "EQUIVALENCE" and tokens[1] == "(" :
			raise error.ParsingError("Invalid EQUIVALENCE statement")

		## FIXME no syntax checking is done here
		self.equ_lists = tokenizer.split_token_list(tokens[1:], [","])
		for i in range(len(self.equ_lists)) :
			self.equ_lists[i] = tokenizer.split_token_list(self.equ_lists[i][1], [","])

		
	def str (self, max_len_first, max_len_following=0) :
		""" Gives this line as a list of strings. The first line with maximum length max_len_first,
                    the following lines of maximum length max_len_following
		    (The default values max_len_following=0 means 'same as first line')
		    For comments, the maximum lengths can be ignored! 
		"""

		result = ["EQUIVALENCE "]
		indent = len(result[0])
		max_len_first -= indent+1                       ## we assume it fits on the line
		if not max_len_following == 0 :
			max_len_following -= indent

		names = []
		for l in self.equ_lists :
			names += ["(", []]
			for var_name in l :
				names[-1] += var_name + [","]
			del names[-1][-1]

		params = tokenizer.join_tokens(names, max_len_first, max_len_following)

		result[0] += params[0]
		for line in params[1:] :
			result.append (indent*" " + line)

		return [result]

	def variables_used (self) :
		""" return a list of the names of the variables in this statement """
		## These names do not contain dimension specification (everything in brackets
		## that comes after a name is am array index - either the arry was declared
		## correctly or it is wrong anyway, there is no implicit declaration of arrays) !

		result = []

		for l in self.equ_lists :
			for var_name in l :
				result.append(var_name[0])
		return result


class executable (unknown) :
	""" Base class for executable statements """

	def __init__ (self) :
		unknown.__init__ (self)
		self.is_executable = True
		self.type = "executable"
	


def new_statement (tokens) :
	""" generic interface to this module:
	    Take a list of tokens and returns a suitable statement-class
	"""

	## a declaration
	if tokens[0].upper() in ["INTEGER", "REAL", "LOGICAL", "DOUBLE PRECISION", "CHARACTER"] :
		return declaration (tokens)
	elif tokens[0].upper() == "COMMON" :
		return common_block (tokens)
	elif tokens[0].upper() == "PARAMETER" :
		return parameter (tokens)
	elif tokens[0].upper() == "EQUIVALENCE" :
		return equivalence (tokens)
	else :
		s = tokenizer.join_tokens(tokens)
		raise error.NotImplementedError("Unknown statement: "+s)

		return unknown (tokens)

def new_declaration (var_names) :
	""" Generates a new 'declaration', declaring all the variables in var_names (list of strings)
	    The implicit-rules are used to determine the type
        """

	tokens = ["::"]
	for n in var_names :
		tokens += tokenizer.tokenize(n) + [","]
	del tokens[-1]

	result = declaration (tokens)

	return result
