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
## variable.py
##------------------------------------------------------
##
## provides the class "variable" which represents one variable 
## 
##------------------------------------------------------

import tokenizer

class variable :
        """ Class to represent a variable in a Fortran program unit

        Attributes:
            name (String)      -- The variable name
            type (String)      -- The type of the variable
            dimension (String) -- The array dimensions or empty

	    char_len (String)  -- For CHARACTER variables: the length
        """

        
        def __init__ (self, decl_token_list, type_token_list=[]) :
                """ decl_token_list:    the list of tokens of the (implicit or explicit) declaration,
                                        beginnig with the variable name
                    type_token_list:    the token list declaring the type in an explicit declaration
                                        ( for example INTEGER, REAL*8 or 'INTEGER, DIMENSION(10)'
                                        if empty, implicit typing [IMPLICIT REAL*8 (A-H,O-Z)] is assumed
                """

                self.dimension = ""
		self.type = ""
                implicit = (type_token_list == [])

                ## first process the declaration part
                if not implicit :
                        type_decl = tokenizer.split_token_list (type_token_list, [','])

                        ## up to now we only consider the type name and a dimension statement
                        ## type name
                        self.type = tokenizer.join_tokens (type_decl[0])

			if type_decl[0][0].upper() == "CHARACTER" :
				self.char_len = '1'
				self.type = "CHARACTER"
				if "*" in type_decl[0] :
					i = type_decl[0].index("*")
					if type_decl[0][i+1] == '(' :
						self.char_len = tokenizer.join_tokens(type_decl[0][i+1:i+3])
					else :
						self.char_len = type_decl[0][i+1]

				elif len(type_decl[0]) > 1 and type_decl[0][1] == "(" :
					if type_decl[0][2][0:1] == ["LEN", "="] :
						self.char_len = int(type_decl[0][2][2])

                        for i in range(1, len(type_decl)) :
                                if type_decl[i][0].upper() == "DIMENSION" :
                                        if len(type_decl[i]) < 3 or not type_decl[i][1] == '(' :
                                                raise ParsingError("Invalid DIMENSION statement")
                                        else :
                                                self.dimension = tokenizer.join_tokens (type_decl[i][1:3])
				if type_decl[0][0].upper() == "CHARACTER" and \
				       type_decl[i][0].upper() == "LEN" and type_decl[i][1] == "=":
					char_len = typle_decl[i][2]

                ## now process the variable name
                self.name = decl_token_list[0]

                ## is there a dimesion specification following the variable name?
                ## (if so, it is overriding a possible previos DIMENSION statement
                if len(decl_token_list) >= 3 and decl_token_list[1] == '(' :
                        self.dimension = tokenizer.join_tokens (decl_token_list[1:3])

		## same for character length
		if self.type == 'CHARACTER' and "*" in decl_token_list:
			i = decl_token_list.index("*")
			if decl_token_list[i+1] == '(' :
				self.char_len = tokenizer.join_tokens(decl_token_list[i+1:i+3])
			else :
				self.char_len = decl_token_list[i+1]

                ## if needed, do implicit typing
                if implicit :
                        if self.name[0].upper() in "IJKLMN" :
                                self.type = "INTEGER"
                        else :
                                self.type = "REAL*8"

        ## for debugging only
        def output (self) :
                print "Type:      ", self.type
                print "Name:      ", self.name
                print "Dimension: ", self.dimension
		if self.type == "CHARACTER" :
			print "char_len:  ", self.char_len   
