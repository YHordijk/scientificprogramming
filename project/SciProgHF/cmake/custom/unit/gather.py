#!/usr/bin/env python

import os
import sys

# get directories in the "folder" in alphabetical order
def listdirs(folder):
  result = []
  for d in os.listdir(folder):
    if os.path.isdir(os.path.join(folder, d)): result.append(os.path.join(folder, d).replace("\\","/"))
    result.sort()
  return result

# set location of unit tests
unit_test_dir = sys.argv[1]

# open file for gathering unit tests
f = open(sys.argv[2], 'w')

# gather unit tests
for subd in listdirs(unit_test_dir):
  if not '.git' in subd:
    for path in listdirs(subd):
      if os.path.isdir(path):
        if os.path.isfile(path + '/test.F90'):
          name = path.split('/')[-2] + '__' + path.split('/')[-1]
          f.write('add_unit_test(%s %s)\n' % (path, name))

f.close()
