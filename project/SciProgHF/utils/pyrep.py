#! /usr/bin/env python        
############################################################################################
#
# Testing script on analyzing xml output files from Dirac
#
# Writen by Miro ILIAS, Banska Bystrica Nov 2011
# Development started upon Ulf's suggestion
#
# http://www.learningpython.com/2008/05/07/elegant-xml-parsing-using-the-elementtree-module/
#
# Further improvement, suggestion etc are warmly welcomed
#
############################################################################################
import  sys
import  os
import  string
import  optparse

import xml.etree.ElementTree as et
#from numpy import *

#########################################################################
#
#                        MAIN FUNCTION
#
#########################################################################
def main(*args):

    dirac_xml="dirac.xml"
    f = open(dirac_xml,'r')

    #f_read = f.read()
    #print f_read 

    xml_data='<?xml version="1.0"?>\
<parent id="top">\
<child1 name="paul">Text goes here</child1>\
<child2 name="fred">More text</child2>\
</parent>'

    xml = et.XML(f.read())
    #element = et.XML(f_read)
    #xml = et.XML(xml_data)

    #dirac_tag='.//mytag'
    #dirac_tag="name"
    dirac_tag='child1'
    #dirac_tag='child2'

    #print xml.findall(dirac_tag)

    for t in xml.findall(dirac_tag):
        print t.text

###########################################################################
###                   process the main function
###########################################################################
if __name__ == '__main__':
    sys.exit(main(*sys.argv))

