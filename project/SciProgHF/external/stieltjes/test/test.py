#!/usr/bin/env python

import sys
import subprocess
import numpy
import os

#prog_dir = '/home/elke/programme/stieltjes/'

binary  = '../build/stieltjes.x'
input   = '../test/sip.dat'
#input   = '../FPSPEC'
output  = 'stieltjes_out.own'

allowed = 5.0e-8


os.chdir('../test')
#subprocess.call('pwd',shell=True)

subprocess.call([binary, input])
subprocess.call('rm -rf fort.*',shell=True)

ref    = numpy.loadtxt('reference')
own    = numpy.loadtxt('inp.dat')

for i in range (0,100):
   diff    = abs(ref[i] - own[i])
   reldiff = diff / ref[i]
   #print reldiff
   if reldiff > allowed:
      print 'exit'
      sys.exit(1)

sys.exit(0)
