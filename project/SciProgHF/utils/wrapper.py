#!/usr/bin/env python

"""
Universal wrapper

Copyright (c) 2011-2016 Ulf Ekstrom

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os, subprocess, shutil, glob, signal, optparse, sys, shlex, tarfile, tempfile, zipfile, socket

def wrap_dirac():
    usage = """%prog [options] MOLECULE.{mol,xyz} INPUT.inp [ARCHIVE.tgz] [FILE ..]
Run Dirac starting from the given mol and inp files. Additional files will be copied
to the temporary directory, archives will be extracted.

If a list of remote hosts is given, either with --host or --hostfile, a temporary
directory will be created on all unique hosts, all files will be copied there,
and the directories will be deleted together with the main tmpdir. No files are
saved from the remote hosts.

Use --launcher="mpirun OPTIONS" with the appropriate OPTIONS to run an mpi job.

Use --launcher="valgrind OPTIONS" to debug with valgrind.

Create a file named 'keep' in the temporary directory if you don't want it to
be deleted at the end of a run.
"""
    W = Wrapper()
    parser = optparse.OptionParser(usage)

    parser.add_option('-v','--verbose', action='store_true', dest="verbose",
                     help='Print more information')

    group = optparse.OptionGroup(parser, 'General options')
#    group.add_option('--mb', type='string', default=None, dest="mb",
#                     help='Amount of memory for Dirac, in MB',metavar="MEGABYTES")
    group.add_option('--executable', type='string', default=None, dest="executable",
                     help='Path to main executable',metavar="FILE")
    group.add_option('-D','--define', type='string', action='callback',callback=W.define_callback,
                     help='Define the environment variable NAME as VALUE',
                     metavar='NAME=VALUE')
    group.add_option('--no-archive', action='store_true', dest="noarchive",default=False,
                     help="Do not create an archive file.")
    group.add_option('--backup', type='int', default=1, dest='backup_level',
                     help="Set backup level for output files [default: %default]. Set to 0 for no backup.")
    group.add_option('-i','--infile', type='string', action='callback',callback=W.infile_callback,
                     help='Copy FILE to the temporary directory, optionally renaming it to NAME',
                     metavar='[NAME=]FILE]')
    group.add_option('-r','--replace', type='string', action='append',default=[],dest="replace",
                     help='Replace NAME by VALUE in the mol and inp files',metavar="NAME=VALUE")
    group.add_option('-g','--getfile', type='string', action='callback',callback=W.getfile_callback,
                     help='Copy FILE from temporary directory after the calculation',
                     metavar='NAME[=FILE]')
    group.add_option('-p', '--pcm', type='string', action='callback', callback=W.pcm_callback,
                     help='PCM calculation: parse PCM input and copy it to scratch')
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Temporary directory options')
    group.add_option('-t','--tmpdir', type='string', default=None,dest="tmpbase",
                     help='Base directory for temporary directory',metavar="PATH")
    group.add_option('-T','--tmpfull', type='string', default=None,dest="tmppath",
                     help='Full path to temporary directory. Will be created if not there.',metavar="PATH")
    group.add_option('-k','--keep-tmp', action='store_true', dest="keeptmp",default=False,
                     help="Keep temporary directory, default if it's preexisting.")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Parallel options')
    group.add_option('--launcher', type='string', default=None,dest="launcher",
                     help='Command used to launch the main executable, i.e. "mpirun -np 2"',metavar="COMMAND")
    group.add_option('--host', type='string', default=None,dest="host",
                     help='Comma separated list of hosts with separate file systems',metavar='"NODE[,NODE, ..]"')
    group.add_option('--hostfile', type='string', default=None,dest="hostfile",
                     help='File to read host list from, one host per line',metavar='FILE')
    parser.add_option_group(group)

    # Read .rc files
    rcargs = []
    rcpaths = ['.wrapperrc','wrapperrc',
               os.path.join(os.path.expanduser('~'),'.wrapperrc')]
    for path in rcpaths:
        try:
            f = open(path,'r')
            for l in f.readlines():
                l = l.strip()
                if len(l) > 0:
                    if l[0] != '#':
                        for a in shlex.split(l):
                            rcargs.append(os.path.expandvars(a))
            break
        except Exception:
            pass

    (options, args) = parser.parse_args(rcargs+sys.argv[1:])
    W.verbose = options.verbose
    W.tmppath = options.tmppath
    W.tmpbase = options.tmpbase
    W.keeptmp = options.keeptmp
    W.archpatterns = ['*.xml','DFCOEF','MOLECULE.*','DIRAC.INP','cavity.off','PEDRA.OUT', 'PCM_mep_asc', 'molec_dyadic.dat', 'cavity.inp', 'dumped1.dat', 'dumped2.dat', 'create_cavity.out']

    # Detect file types in args from the file names, set up job name
    molpath = None
    inppath = None
    jobname = None
    def stripsuffix(a,ss):
        "If a ends in one of the strings in ss, return a with the suffix removed, else return None"
        for s in ss:
            if a[-len(s):] == s:
                return a[:-len(s)]
        return None
    for f in args:
        if stripsuffix(f,['.mol','.xyz']):
            if molpath is not None:
                print >> sys.stderr, "Wrapper: Multiple molecule files given, aborting."
                sys.exit(-1)
            molpath = f
        elif stripsuffix(f,['.inp']):
            if inppath is not None:
                print >> sys.stderr, "Wrapper: Multiple .inp files given, aborting."
                sys.exit(-1)
            inppath = f
        elif stripsuffix(f,['.tgz','.tar.gz','.zip']):
            W.extract.append(f)
            jobname = stripsuffix(os.path.basename(f),['.tgz','.tar.gz','.zip'])
        else: # Just copy in unknown files
            W.copyin.append([f])
    if molpath is not None and inppath is not None:
        jobname = (stripsuffix(os.path.basename(inppath),['.inp']) + '_' +
                   stripsuffix(os.path.basename(molpath),['.mol','.xyz']))
    elif jobname is None:
        print >> sys.stderr, "Wrapper: Need either molecule and input files or archive, aborting."
        sys.exit(-1)
    jobname = "_".join([jobname]+options.replace)

    replacements = reduce(list.__add__,[x.split("=") for x in options.replace],[])
    if molpath is not None:
        if stripsuffix(molpath,['.mol']):
            W.copyin.append([molpath,'MOLECULE.MOL']+replacements)
        else:
            W.copyin.append([molpath,'MOLECULE.XYZ']+replacements)
    if inppath is not None:
        W.copyin.append([inppath,'DIRAC.INP']+replacements)

    stdout_file_name = jobname + '.out'
    if options.backup_level > 0:
        # backup output file before writing to it
        if os.path.isfile(stdout_file_name):
            if options.backup_level > 1:
                for i in range(options.backup_level-2, -1, -1):
                    src = stdout_file_name + '_' + str(i)
                    dst = stdout_file_name + '_' + str(i+1)
                    if os.path.isfile(src):
                        shutil.copy2(src, dst)
            src = stdout_file_name
            dst = stdout_file_name + '_0'
            shutil.copy2(src, dst)

    W.stdout = open(stdout_file_name,'w')
    W.stderr = open(jobname+'.err','w')

    W.stdin = None
    if not options.noarchive:
        W.archive = jobname+'.tgz'

    # Append executable itself to the commands to be run
    if options.executable is None:
        print >> sys.stderr,"Wrapper: No executable given (with --executable=), aborting."
        sys.exit(-1)
    else:
        exename = os.path.basename(options.executable)
        W.copyin.append([options.executable,exename])
        if options.launcher is None:
            W.cmdlines.append(['./'+exename])
        else:
            W.cmdlines.append(shlex.split(options.launcher)+['./'+exename])
    # Parse remote host specifications
    if options.host is not None:
        W.remote_hosts = set(options.host.split(','))
    if options.hostfile is not None:
        for items in [l.split() for l in open(options.hostfile,'r').readlines()]:
            if items:
                W.remote_hosts.add(items[0])
    if socket.gethostname() in W.remote_hosts: # Don't copy to this machine
        W.remote_hosts.remove(socket.gethostname())

    # Now start the real work
    try:
        W.setup()
    except:
        raise # TODO: when bug free just exit with an error code
    try:
        W.stagein()
        try:
            W.run()
        finally:
            W.stageout()
    except:
        raise # TODO: when bug free just exit with an error code
    finally:
        W.cleanup()
        W.stderr.close()
        if os.path.getsize(jobname+'.err') == 0:
            os.remove(jobname+'.err')
        W.stdout.close()

class Wrapper:
    def __init__(self):
        # User modifiable variables
        self.tmpbase = None
        self.tmppath = None
        self.keeptmp = False
        self.stdin = None
        self.stdout = None
        self.stderr = None
        self.cmdlines = []
        self.copyin = []
        self.copyout = []
        self.extract = []
        self.archive = None
        self.archpatterns = []
        self.verbose = False
        self.workdir = None
        self.env = os.environ.copy()
        self.remote_hosts = set() # Duplicate tmpdir and stage in on these hosts
        self.rsh = 'ssh'
        self.rcp = 'scp'
        self.do_pcm = False

    def define_callback(self, option, opt, value, parser):
        if '=' in value:
            ff = value.split('=')
            self.env[ff[0]] = '='.join(ff[1:])
        else:
            print >> sys.stderr, "Wrapper: Illegal argument to",opt,"expected NAME=VALUE"
            sys.exit(1)

    def infile_callback(self, option, opt, value, parser):
        if '=' in value:
            ff = value.split('=')
            if len(ff) != 2:
                parser.error("Illegal argument to "+opt)
                print >> sys.stderr, "Wrapper: Illegal argument to",opt
                sys.exit(1)
            self.copyin.append(ff)
        else:
            self.copyin.append([value])

    def getfile_callback(self, option, opt, value, parser):
        if '=' in value:
            ff = value.split('=')
            if len(ff) != 2:
                print >> sys.stderr, "Wrapper: Illegal argument to",opt
                sys.exit(1)
            self.copyout.append(ff)
        else:
            self.copyout.append([value])

    def command_callback(self, option, opt, value, parser):
        self.cmdlines.append(shlex.split(value))

    def pcm_callback(self, option, opt, value, parser):
        if value is None:
            sys.exit('No .pcm file specified')
        else:
            self.copyin.append([value,'pcmsolver.inp'])
        if 'PCMSOLVER_PARSE' in os.environ.keys():
            self.copyin.append([os.environ['PCMSOLVER_PARSE'], 'pcmsolver'])
        else:
            self.copyin.append(['pcmsolver','pcmsolver'])
#        self.archpatterns.append('PEDRA.OUT')
#        self.archpatterns.append('cavity.off')
        self.do_pcm = True

    def setup(self):
        try:
            if self.tmppath is not None:
                self.tmppath = os.path.abspath(self.tmppath)
                if os.path.isdir(self.tmppath):
                    self.keeptmp = True
                else:
                    os.makedirs(self.tmppath)
            else:
                if self.tmpbase is None:
                    self.tmppath = tempfile.mkdtemp()
                else:
                    self.tmppath = tempfile.mkdtemp(dir=self.tmpbase)
            if self.verbose: print "Wrapper: Using temporary directory '%s'" % self.tmppath
            self.workdir = os.getcwd()
        except:
            print >> sys.stderr, "Wrapper: Cannot create temporary directory '%s', quitting." % self.tmppath
            raise
        try:
            for h in self.remote_hosts:
                self.run_cmd([self.rsh,h,'mkdir -p "'+self.tmppath+'"'])
        except subprocess.CalledProcessError:
            print >> sys.stderr, "Wrapper: Cannot create remote temporary directory on %s, quitting." % h
            for hh in self.remote_hosts:
                try:
                    self.run_cmd([self.rsh,hh,'rm -rf "'+self.tmppath+'"'],silent=True)
                except subprocess.CalledProcessError:
                    pass
            raise

    def stagein(self):
        "Copy in files and extract archives. Stay in cwd. Extract first to allow copyin to override archive files."
        # Extract archives, stop if problem
        def hassuffixes(a,ss):
            for s in ss:
                if a[-len(s):] == s:
                    return True
            return False
        for a in self.extract:
            if self.verbose: print "Wrapper: Extracting archive",a
            try:
                if hassuffixes(a,['.tar']):
                    t = tarfile.open(a,'r')
                elif hassuffixes(a,['.tgz','.tar.gz']):
                    t = tarfile.open(a,'r:gz')
                elif hassuffixes(a,['.tbz','.tar.bz2']):
                    t = tarfile.open(a,'r:bz2')
                elif hassuffixes(a,['.zip']):
                    t = zipfile.ZipFile(a)
                else:
                    print >> sys.stderr, "Wrapper: Cannot determine type of archive file",a,"quitting."
                    sys.exit(1)
                t.extractall(path=self.tmppath)
            except:
                print >> sys.stderr, "Wrapper: Error extracting files from '"+a+"', quitting."
                raise
        # Copy in files, stop if problem
        for f in self.copyin:
            src = f[0]
            if os.path.isfile(src):
                if len(f) % 2 == 0: # dst name given
                    dst = os.path.join(self.tmppath,f[1])
                    rep = f[2:]
                else:
                    dst = os.path.join(self.tmppath,os.path.basename(src))
                    rep = f[1:]
                rep = zip(rep[0::2],rep[1::2])
                try:
                    self.copy_with_replacement(src,dst,rep)
                except:
                    print >> sys.stderr, "Wrapper: File copy failed, quitting."
                    raise
            else:
                print >> sys.stderr, "Wrapper: File copy of %s failed, source file nonexistent." % src
        for h in self.remote_hosts:
            try:
                self.run_cmd([self.rcp]+glob.glob(os.path.join(self.tmppath,'*'))+[h+':'+self.tmppath])
            except subprocess.CalledProcessError:
                print >> sys.stderr, "Wrapper: Cannot copy files to host",h,"quitting."
                raise

    def copy_with_replacement(self,src,dst,rep=[]):
        try:
            if self.verbose: print "Wrapper: Copying '"+src+"' to '"+dst+"'"
            if len(rep) == 0:
                if src.endswith('.gz') and not dst.endswith('.gz'):
                    # -ifile.gz=file
                    # guzip file after copyin
                    shutil.copy2(src, dst+'.gz')
                    os.system('gunzip %s.gz' % dst)
                else:
                    shutil.copy2(src,dst)
            else:
                if self.verbose:
                    for r in rep:
                        print "Wrapper: While replacing '"+r[0]+"' by '"+r[1]+"'"
                of = open(dst,'w')
                for l in open(src,'r').readlines():
                    for r in rep:
                        l = l.replace(r[0],r[1])
                    of.write(l)
                of.close()
        except IOError,e:
            print >> sys.stderr, "Wrapper: Error copying file '"+src+"'"
            raise

    def stageout(self):
        # Create archive, easiest from tmppath
        os.chdir(self.tmppath)
        if self.archive is not None:
            try:
                self.archive = os.path.join(self.workdir,self.archive)
                if self.verbose: print "Wrapper: Creating archive",self.archive
                t = tarfile.open(self.archive,'w:gz')
                for p in self.archpatterns:
                    for f in glob.glob(p):
                        t.add(f)
                t.close()
            except:
               print >> sys.stderr, "Wrapper: Error creating tar archive"
               raise # TODO: Should try to copy out files as well
        # Copy back files
        def cpout(src,dst=None):
            if self.verbose: print "Wrapper: Copying back",src
            if os.path.isfile(src):
                if dst is None: dst='.'
                shutil.copy2(src,dst)
            elif os.path.isdir(src):
                if dst is None: dst=os.path.basename(src)
                if os.path.exists(dst):
                    shutil.rmtree(dst)
                shutil.copytree(src,dst)
            elif self.verbose: # Just give a warning, don't panic
                print "Wrapper: When copying out, file",src,"not found."
        os.chdir(self.workdir)
        for f in self.copyout:
            try:
                if len(f) == 2:
                    if f[1].endswith('.gz') and not f[0].endswith('.gz'):
                        # -gfile=file.gz
                        # gzip file and add file.gz=file.gz to cpout
                        os.system('gzip %s' % os.path.join(self.tmppath, f[0]))
                        cpout(os.path.join(self.tmppath,f[0]+'.gz'),f[1])
                    else:
                        cpout(os.path.join(self.tmppath,f[0]),f[1])
                else: # allow patterns in copyout if there is no =
                    for ff in glob.glob(os.path.join(self.tmppath,f[0])):
                        cpout(ff)
            except:
                print >> sys.stderr,"Wrapper: Error copying out file",f
                raise # TODO: Must try to copy out all files

    def cleanup(self):
        os.chdir(self.workdir)
        if not self.keeptmp and not os.path.isfile(os.path.join(self.tmppath,'keep')):
            shutil.rmtree(self.tmppath)
            if self.verbose: print "Wrapper: Temporary directory deleted"
            for h in self.remote_hosts:
                try:
                    self.run_cmd([self.rsh,h,'rm -rf "'+self.tmppath+'"'])
                except subprocess.CalledProcessError:
                    pass
        else:
            if self.verbose: print "Wrapper: Temporary directory kept:",self.tmppath

    def run_cmd(self,cmd,silent=False):
        "Run command line cmd (a list of arguments), return exit code"
        pgid = None
        oldhand = []
        try:
            cmd = "'"+"' '".join(cmd)+"'" # Need this because we run with shell=True
            def handler(signum,frame):
                if pgid is not None:
                    if signum == signal.SIGALRM: signum = signal.SIGTERM
                    os.killpg(pgid,signum)
            catch = [signal.SIGTERM, signal.SIGHUP, signal.SIGINT, signal.SIGALRM]
            for s in catch:
                oldhand.append(signal.signal(s,handler))
            if self.verbose:
                print "Wrapper: Executing",cmd
            p = subprocess.Popen(cmd,stdin=self.stdin,stdout=self.stdout,stderr=self.stderr,shell=True,
                                 preexec_fn=os.setsid,env=self.env)
            pgid = p.pid
            p.wait() # OSError will be raised here if killed
            if p.returncode != 0:
                raise subprocess.CalledProcessError(p.returncode,cmd)
        except subprocess.CalledProcessError, e:
            if not silent:
                print >> sys.stderr, "Wrapper: Command "+str(cmd)+" failed with returncode",e.returncode
            raise
        except OSError, e: # Typically signal
            if not silent:
                print >> sys.stderr, "Wrapper: Got OSError %i while executing" % e.errno,cmd
            raise subprocess.CalledProcessError(e.errno,cmd)
        finally: # restore signal handlers
            if len(oldhand) > 0:
                for s,h in zip(catch,oldhand):
                    signal.signal(s,h)

    def run(self):
        try:
            os.chdir(self.tmppath)
            if self.do_pcm:
                os.system('./pcmsolver')
            for cmd in self.cmdlines:
                self.run_cmd(cmd)
        finally:
            os.chdir(self.workdir)

wrap_dirac()
