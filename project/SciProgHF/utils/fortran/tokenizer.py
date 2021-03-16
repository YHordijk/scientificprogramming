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
## tokenizer.py
##------------------------------------------------------
##
## A simple tokenizer. It provides:
## - tokenize (line) :
##           string -> list of tokens
##
## - join_tokens (token_list) :
##           list of tokens -> string
##           (optional arguments allow to give maximum lengths for the first line and
##            for the following lines. If this is used, a list of strings is returned
##            instead of a single string)
##
## - split_token_list(token_list, separator_list) :
##           list of tokens -> list of lists of tokens, split at the given seaprators
##           The second argument is a list of separators!
##
## Error handling:
##
##------------------------------------------------------


##------------------------------------------------------
## KNOWN PROBLEMS:
##  * comments in lines, that are continued on the following line, are handled completely wrong!
##    [e.g.: INTEGER I, & ! one integer
##                   J    ! another one ]
##    Its really difficult to handle this, so I wont try.
##    FIXME: maybe add a check?
## 
##------------------------------------------------------

## FIXME: reorganize this whole file! It became far too complicated !

import error

def tokenize (line) :
        """ a simple tokenizer for Fortran lines.
            It takes a line and returns a (nested) list of tokens
            It assumes to get correct and 'resonable' Fortran code
            Example: 'INTEGER INTA, INTB(10), INTC(1:5,10)' gives
                     [ 'INTEGER', 'INTA', ',' , 'INTB', '(', ['10'], ',' , 'INTC', '(', ['1',':','5', ',' , '10'] ]

            This means: Separators are included as tokens (if they are non-whitespace)
                        Opening brackets are included and indicate that the next element in the list is
                        a sublist of the tokens enclosed in brackets
                        Closing brackets are not included

            The original line (or something quite similar to it) can be reconstructed from this list
            with join_tokens (list)

	    Error handling:
	                tokenize does not raise any exceptions (if it does, its a bug)

        """

        result = []

        sep = [' ', ',', ':', '\n', '+', '-', '*', '/', '\\', '=']
        start = 0
        in_token = False; in_brackets = False; in_string = False
        for i in range(len(line)) :
                if in_brackets :
                        if line[i] == '(' :
                                open_brackets += 1
                        elif line[i] == ')' :
                                if open_brackets > 0 :
                                        open_brackets -= 1
                                # we have found the matching closing bracket
                                else :
                                        result.append (tokenize (line[start:i]))
                                        start=i+1
                                        in_brackets = False
		elif in_string :
			if line[i] == string_start :
				result.append ( line[start:i+1] )
				start = i+1
				in_string = False

                elif in_token :
                        if line[i] in sep + ['(', '"', "'"] :
                                result.append ( line[start:i] )
                                start = i
                                in_token = False

                # we are looking for the start of a new token
                if i >= start and not in_token and not in_brackets and not in_string:
                        # start of section in brackets
                        if line[i] == '(' :
                                result.append ('(')
                                start = i+1
                                open_brackets = 0
                                in_brackets = True
			# start of new string
			elif line[i] in ["'", '"'] :
				start = i
				string_start = line[i]  # remember if string starts with ' or " 
				in_string = True
			# start of a ! comment
			elif line[i] == '!' :
				start = i;
				break
			# separators
			elif line[i] in sep :
				if not line[i] == ' ' :
					result.append(line[i])
				start = i+1
                        # start of a new token
                        elif not line[i] == ')':
                                start = i
                                in_token = True

        # append last token
        if not line[start:].strip ('\n ') == "":
                result.append ( line[start:] )

        ## now we do some post-processing:
        i = 0
        while i < len(result)-1 :
                ## "DOUBLE PRECISION" is one token, not two
                if result[i:i+2] == ["DOUBLE", "PRECISION"] :
                        result[i] = "DOUBLE PRECISION"
                        del result[i+1]
                ## "::" is one token
                elif result[i:i+2] == [":", ":"] :
                        result[i] == "::"
                        del result[i+1]
                i += 1

	return result;

def join_tokens (list, max_len_first = 0, max_len_following = 0) :
        """ Takes a list of tokens (same syntax as the ones produced by tokenize
            and reconstructs the corresponding string, nicely formatted

	    If max_len_first != 0, than a list of strings, each with maximum length max_len is returned
	    max_len_following is the maximum length of lines following the first one
	    (max_len_following defaults to max_len_first)
        """
        
        result = [""]
	if max_len_following == 0 :
		max_len_following = max_len_first

        def space (str) :
                # no space after these characters
                if not str == "" and not str[-1] in [' ', ':', '(', '+', '-', '*', '/', '\\'] :
                        res = ' '
		else :
			res = ''
                return res;

        in_brackets = False
	max_len = max_len_first

        for token in list :
                if in_brackets :
			if max_len == 0 :
				result[-1] += space(result[-1]) + '(' + join_tokens(token) + ')'
			else :
				## first we try if we can get the stuff in brackets (or a part of it)
				## on the current line
				try:
					str = space (result[-1]) + '('
					## the -2 is for the enclosing brackets

					token_str = join_tokens(token, max_len-len(result[-1]+str)-2, max_len_following)

					result[-1] += str + token_str[0]
					if len(token_str) > 1 :
						result += token_str[1:] 
					result[-1] += ')'

				except TokenizerError :
				## There is no way to put it on this line, use a new one
				## We dont catch exceptions here, if one happens, something is
				## wrong with the Fortran code

					## the -2 is for the enclosing brackets 
					token_str = join_tokens (token, max_len - 2, max_len_following)
					result.append( '(' + token_str[0] )
					if len(token_str) > 1 :
						result += token_str[1:] 
					result[-1] += ')'
					max_len = max_len_following
					
                        in_brackets = False
                # no space before these characters
                elif token in [',', ':', '+', '-', '*', '/', '\\'] :
                        result[-1] += token
                elif token == "(" :
                        in_brackets = True
                else :
		# a normal token
			token_str = space(result[-1]) + token
			## its a comment: we can ignore max_len
			if token[0] == '!' :
				result[-1] += token_str
			## the +1 is for a possibly following ','  
			elif max_len == 0 or len(result[-1] + token_str) + 1 <= max_len :
				result[-1] += token_str
			else :
				## handle strings seperatly, because we can break a string everywere!
				## (but we only want to do it if it is really nessacary
				## and if it is, we want to do it nicely)
				if token[0] in ['"', "'"] :
					## The easiest way is to put the complete string on the next line
					print token, len(token), max_len_following
					if len(token)+1 <= max_len_following :
						max_len = max_len_following
						result.append(token)
					else:
						space_pos =[0]; word_length = []
						res = 0
						while not res == -1 :
							res = token_str.find(' ', res+1)
							if not res == -1 : 
								word_length.append(res - space_pos[-1])
								space_pos.append(res)

						if max(word_length)+1 < max_len_following :
							## FIXME: handle string breaking over mutliple lines here!

							max_len_following
							result.append(token)
							
							raise error.NotImplementedError("tokenizer.join_tokens: " \
										  "breaking of strings over mutiple lines")

						else :
					                ## if nothing helps, we use brute force
							i = max_len - len(result[-1])
							result[-1] += token_str[0:i]
							token_str = token_str[i:] 
							while not token_str == "" :
								max_len = max_len_following
								result.append (token_str[0:max_len_following])
								token_str = token_str[max_len_following:] 	
				else :
					max_len = max_len_following
					if len(token)+1 > max_len :
						raise error.TokenizerError("max_len too small in join_tokens")
					result.append (token)

	if max_len == 0 :
		return result[0]
	else:
		return result

def split_token_list (list, sep_list) :

        result = []
        start = 0

        try:
                for i in range(len(list)) :
                        if list[i] in sep_list :
                                result.append (list[start:i])
                                start = i+1
        except TypeError:
                raise error.Bug, "split_token_list has to be called with a LIST of separators"
                        
        result.append (list[start:])

        return result
