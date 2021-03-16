#! /usr/bin/env python
#
# Inform user who compiles DIRAC about (non-)existing dirarc file(s)
#
#  install_dir comes in as parameter
#     #install_dir = '${CMAKE_SOURCE_DIR}'
#
import sys
import os
import string
import glob
import optparse
import shlex
import shutil
#
global install_dir
#
def handle_diracrc():
    global install_dir
    fullargs = []
    rc_file_found=None
    ret_var=None
# hierarchy - from closest to the farthest
    paths = ['.diracrc','diracrc',os.path.join(install_dir,'.diracrc'),
          os.path.join(install_dir,'diracrc'),
          os.path.join(os.path.expanduser('~'),'.diracrc')]
    for path in paths:
        try:
            f = open(path,'r')
            read_args=False
            rc_file_found=path
            for l in f.readlines():
                l = l.strip()
                if len(l) > 0:
                    if l[0] != '#':
                        for a in shlex.split(l):
                            fullargs.append(os.path.expandvars(a))
                            read_args=True
            if read_args:
                break
        except Exception:
            pass
    if rc_file_found:
        print '-- cmake python script found DIRAC configuration file:',rc_file_found
    if fullargs: # save info about the diracrc file
        print '-- this DIRAC config file contains these flags for pam:',fullargs
    else:
        if rc_file_found:
            print '-- but this file contains no active flags for pam - fix it !'
        else:
            dest=os.path.join(install_dir,"cmake","diracrc.in")
            src=os.path.join(install_dir,"diracrc")
            print '-- no dirac cofiguration file found - created default template:',src
            shutil.copy(dest,src)

def read_flags():
    global install_dir
    from optparse import OptionParser, OptionGroup
    usage = "%prog [options]"
    parser = OptionParser(usage)
    group = OptionGroup(parser, 'mandatory specifications:')
    group.add_option('--install', type='string', action='store',
           dest='dest_install_dir', help='installation directory',metavar='DIRECTORY_STRING')
    parser.add_option_group(group)
    (options, args) = parser.parse_args()
    if options.dest_install_dir:
        install_dir=options.dest_install_dir
        if not os.path.isdir(install_dir):
            print 'given install dir:',install_dir
            print 'good install dir ? ',os.path.isdir(install_dir)
            sys.exit("install dir is wrong!")
    else:
        sys.exit("install dir not specified by flag !")

# "main" part
read_flags()
handle_diracrc()
