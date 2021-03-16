#!/usr/bin/env python

# script to generate a C header and Fortran module file
# with the current settings of XCFun

import re

f = open('src/functionals/list_of_parameters.h', 'r')
s = f.readlines()
f.close()

regex = re.compile('PARAM\((.*)\),')

para_list = []
for line in s:
    line_without_comment = line.split(r'//')[0]
    if regex.search(line_without_comment):
        para_list.append(regex.search(line_without_comment).group(1))

s  = 'module xcfun_autogen\n'
s += '!  this file is generated by generate_interfaces.py - do not edit\n'
s += '   implicit none\n'
s += '   integer, parameter :: XC_NR_PARAMS = %i\n' % len(para_list)
for i in range(len(para_list)):
    s += '   integer, parameter :: %s = %i\n' % (para_list[i], i+1)
s += 'end module\n'

f = open('fortran/xcfun_autogen.F90', 'w')
f.write(s)
f.close()

s  = '// this file is generated by generate_interfaces.py - do not edit\n'
s += 'enum xcfun_parameters {\n'
for para in para_list:
    s += '%s,\n' % para
s += 'XC_NR_PARAMS\n};\n'

f = open('include/xcfun_autogen.h', 'w')
f.write(s)
f.close()