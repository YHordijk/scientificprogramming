#! /usr/bin/env python        
############################################################################################
#
# Script for preparing and submitting DIRAC pam job into a specified  queue system.
#
# Writen by Miro ILIAS, Tel Aviv, 2008 after Luuk's pamq script
#
# Further improvement, suggestion etc are warmly welcomed
#
############################################################################################
import  sys
import  os
import  string
import  optparse

global install_dir
install_dir_default = os.path.realpath(__file__)[:-8]
install_dir = install_dir_default

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class PamQVariables:
  # ... class containing all important variables, methods for submitting 
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #              ... initialize set of variables ...
  #
    global install_dir # comes from configure

    if ( not os.path.isdir(install_dir)) :
      print 'Error exit ! Installation directory install_dir does not exit !'
      print '....variable install_dir=',install_dir
      sys.exit()

    ### contains installation directory ...
    self.dirac_home=install_dir
    
    self.dirac_home=os.path.dirname(self.dirac_home)

    # print level
    self.iprint=1

    self.pam_arguments=''
    self.submit_job=True

    # ... specify the default job name ...
    self.job_name='queue_job.'+str(os.getpid())

    # ... specify the default name of the batch script ...
    self.queue_script="queue_job."+str(os.getpid())+".batch_script"

    # ... specify the default batch script output ...
    self.job_out='queue_job.'+str(os.getpid())+'.log_out'

    self.queue_system='' # must be configured

    self.queue_name='' # must be configured
    self.queue_name_specified=False

    # walltime, deafult None 
    self.queue_time=None

    # ... is pure serial job (no mpi added to pam)
    self.is_serial=False

    self.nnodes=0  # must be configured 
    self.proc_per_node=None # can be None
    self.mpi_procs=None # calculated later as nnodes*proc_per_node

    # ... dirac tests 
    self.dirac_tests_run=False
    self.testrun_arguments=None

    # ... what is going to be in script
    self.batch_script_command=None

    # ... set home directory 
    self.home_dir=os.getenv("PWD")

    # ... default DIRAC pam from the current (as best, from the build) directory
    self.pam_script=os.path.join(self.home_dir,'pam') 

    # ... default: script uses own -oe stdout+stderr log file
    self.is_own_queue_stdout=True

    # ... configuration file with selected defaults ...
    self.prepare_config_file=False

    #-------------------------------------------------------------------
    #   ... list of implemented queue system for the script:
    #-------------------------------------------------------------------
    self.implemented_queue_systems=['pbs','bsub','loadleveler','sge']

    # ... set up file with default variables
    self.pamq_defaults_file=self.dirac_home+'/PAMQ_DEFAULTS'

    #  ... read defaults from configure file ...
    self.read_from_config_file=False
    if os.path.isfile(self.pamq_defaults_file):
      self.read_configure_file()
      self.read_from_config_file=True

    if self.iprint >8 :
      print ' setting home directory:',self.home_dir
      print 'INITIALIZATION of class object done'

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def set_pam_mpi(self):
  #
  # class method: set up 'mpi' parameter for the pam script
  #
    if self.iprint > 7:
      print '\nsetting mpi for pam' 
      print 'set_pam_mpi: extracted pam arguments:', self.pam_arguments

    # extend pam flags with the necessary -mpi arguments
    if self.proc_per_node:
      self.mpi_procs=self.nnodes*self.proc_per_node
    else: # self.proc_per_node not defined, so take all CPU's at hand:)
      self.mpi_procs=str("$NUM_ALL_CORES")

    if self.iprint >= 5:
      print 'set_pam_mpi: calculated # mpi=',self.mpi_procs
      print 'using self.nnodes=',self.nnodes
      print 'and self.proc_per_node=',self.proc_per_node

    # ... add mpi stuff, always ?
    self.pam_arguments=' --mpi='+str(self.mpi_procs)+self.pam_arguments

    if self.iprint >= 7:
      print 'set_pam_mpi: mpi-extended pam arguments:', self.pam_arguments

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def create_configure_file(self):
  #-----------------------------------------------------------------------------------------------
  #
  #      class method: create configure file with important parameters as DEFAULTS
  #
  # since we are inside class method definition, so we are to use 'self.' specification
  #
  #-----------------------------------------------------------------------------------------------
    if self.iprint >= 12:
      print 'preparing PAMQ_DEFAULTS file'
    
    error_exit=False
 
    f = open(self.pamq_defaults_file, 'w') 

    print '\nOpening file for saving defaults in:',self.pamq_defaults_file

    f.write('###  \n')
    f.write('###  SOME IMPORTANT DEFAULT VARIABLES FOR THE "pamq.py" SCRIPT ###\n')
    f.write('###  \n')

    if ( (len(self.queue_system)>0) and  (self.queue_system in self.implemented_queue_systems)) :
      f.write('# ... default queue system \n')
      f.write('QUEUE_SYSTEM='+self.queue_system+'\n')
      print 'saved QUEUE_SYSTEM='+self.queue_system
    else :
      print '\ncreate_configure_file: YOU MUST SPECIFY QUEUE SYSTEM!\n'
      print 'Implemented queue systems=',self.implemented_queue_systems
      error_exit=True 

    if self.queue_name:
      f.write('# ... default queue name \n')
      f.write('QUEUE_NAME='+self.queue_name+'\n')
      print 'saved QUEUE_NAME=',self.queue_name

    if self.nnodes > 0:
      f.write('# ... default number of nodes for (parallel) run\n')
      f.write('NODES='+str(self.nnodes)+'\n')
      print 'saved NODES=',self.nnodes
    else :
      print 'FYI: you did not set default ppn value, never mind... '

    if self.proc_per_node:
      f.write('# ... default number of processors (CPU) per 1 node\n')
      f.write('PROCESSORS_PER_NODE='+str(self.proc_per_node)+'\n')
      print 'saved PROCESSORS_PER_NODE=',self.proc_per_node
    else :
      print 'FYI:default number of processors (CPU) per node - ppn - not set, never mind...'

    if (len(self.pam_script)>0) :
      f.write('# ... default location of the DIRAC pam script \n')
      f.write('PAM_SCRIPT='+self.pam_script+'\n')
      print 'saved PAM_SCRIPT=',self.pam_script

    if (self.is_serial) :
      f.write('# .... default for serial jobs !\n')
      f.write('SERIAL_JOBS_ONLY\n')
      print 'saved SERIAL_JOBS_ONLY (ie set only serial runs)'

    f.close()

    if (error_exit) :
      print('script error exit !')
      sys.exit()
    else :
      print(self.pamq_defaults_file+' file with default ready.')


  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def read_configure_file(self):
  #
  # class method: read parameters (DEFAULTS) from the configure file
  #
  #  we are inside class method definition, so we are to use 'self.' specification
  #
    if self.iprint >= 12 :
      print 'in ... read_configure_file'

    file = open(self.pamq_defaults_file,'r')

    # read line after line
    for line in file:
      #print line
      if (line.find('QUEUE_SYSTEM') >= 0 and not line.find('#') >= 0):
         queue_system=line[line.find('=')+1:len(line)]
         #print 'read QUEUE_SYSTEM=',queue_system
         self.save_queue_system(queue_system) 
    
      if (line.find('QUEUE_NAME') >= 0 and not line.find('#') >= 0):
         queue_name=line[ line.find('=')+1:len(line)]
        # print 'read QUEUE_NAME=',queue_name
         self.queue_name=queue_name
         self.queue_name_specified=True

      if (line.find('NODES') >= 0 and not line.find('#') >= 0):
         nodes=line[ line.find('=')+1:len(line)]
         #print 'read NODES=',int(nodes)
         self.nnodes=nodes

      if (line.find('PROCESSORS_PER_NODE') >= 0 and not line.find('#') >= 0):
         ppn=line[line.find('=')+1:len(line)]
         #print 'read PPN=',int(ppn)
         self.proc_per_node=ppn

      if (line.find('PAM_SCRIPT') >= 0 and not line.find('#') >= 0):
         pam_script=line[line.find('=')+1:len(line)]
        # print 'read PAM_SCRIPT dir=',pam_script
         self.save_pam_script(pam_script,True)

      if (line.find('SERIAL_JOBS_ONLY') >= 0 and not line.find('#') >= 0):
        # print 'serial runs as default !'
         self.is_serial=True

    file.close()

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def save_queue_system(self,qsys_in):
  #
  #  class method: check and save queue_system into proper variable
  #
  #  we are inside class method definition, so we are to use 'self.' specification
  #

    # check the ending character !
    if qsys_in.endswith('\n'):
     # print 'ending with EOL !' 
      qsys_in=qsys_in[:-1]

    # check if the queue system is in the list ...
    if (qsys_in not in self.implemented_queue_systems) :
      print '\nsave_queue_system: entering queue system, qsys_in=',qsys_in
      print 'self.implemented_queue_systems=',self.implemented_queue_systems
      print 'This queue system is not implemented ! error exit !'
      sys.exit()
    else :
      #print 'save_queue_system: entering queue system VERYFIED!'
      self.queue_system=qsys_in

    if self.iprint >= 12:
      print 'save_queue_system: saved class member queue_system=',self.queue_system

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def save_pam_script(self,pam_script_in,from_file=False):
  #
  #  class method: check and save queue_system into proper variable
  #
  #  we are inside class method definition, so we are to use 'self.' specification
  #

    # check the ending character (if reading from file)!
    if pam_script_in.endswith('\n'):
      #print 'pam_script ending with EOL !' 
      pam_script_in=pam_script_in[:-1]

    # ... verify the existence of the pam script !
    if os.path.isfile(pam_script_in) :
      #print 'pam script is present!'
      self.pam_script=pam_script_in
    else:
      print 'save_pam_script: entering pam_script_in=',pam_script_in
      print 'save_pam_script: this pam script (with full path) does not exist! Error exit !'
      if (from_file) :
        print 'Bad pam_script value in the file PAMQ_DEFAULTS !' 
      sys.exit()

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def prepare_batch_script_command(self):
  #
  #   class method: prepare final script command
  #
  #  we are inside class method definition, so we are to use 'self.' specification
  #
  # ... either run test or standard dirac run
  #

  # unite into one local string
    self.pam_script_args=self.pam_script+self.pam_arguments
 
    if self.dirac_tests_run :
      if os.path.isfile(self.testlast_path): 
        if self.is_serial :
          self.batch_script_command=self.testlast_path+self.testrun_arguments
        else :
          # ... parallel run of tests ...
          self.batch_script_command=self.testlast_path+' --mpi='+str(self.mpi_procs)+' '+self.testrun_arguments
          self.submit_job=True
      else :
        print 'prepare_batch_script_command: DIRAC testlast script NOT found in:\n',self.testlast_path  
        sys.exit()
    else :
      self.batch_script_command=self.pam_script_args # complete pam command with arguments

    if self.iprint > 7:
      print 'prepare_batch_script_command:',self.batch_script_command
  

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pbs_script():
#==========================================================================================
#
#
#   prepare the script file for PBS queueing system script and submit it (if specified)
#
#
#==========================================================================================
  global GlobalVars

  if GlobalVars.iprint >8 :
    print '\nin prepare_pbs_script...'
    print ' GlobalVars.pam_script_args=',GlobalVars.batch_script_command

  f=open(GlobalVars.queue_script,'w') # open file with many queue commands

  s = ''
  s += '''\
#!/bin/bash
#
#         Queue controls are directly readable from the submitting batch script
#----------------------------------------------------------------------------------------
#
#
#     job name
#
#PBS -N %s''' % GlobalVars.job_name+str('\n')
# continue...
  if GlobalVars.proc_per_node:
    s += '''\
#
# number nodes (NOTE: different from number of processors specified by -mpi pam flag !)
# plus number of CPUs per node
#
#PBS -l nodes='''+str(GlobalVars.nnodes)+''':ppn='''+str(GlobalVars.proc_per_node)+str('\n')
  else:
    s += '''\
#
# number nodes only (NOTE: different from number of processors specified by -mpi pam flag !)
#
#PBS -l nodes='''+str(GlobalVars.nnodes)+str('\n')
#
#   if the queue_name is specified in the input, add it
#
  if GlobalVars.queue_name:
    s += '''\
#
# specify the queue name
#
#PBS -q %s ''' % GlobalVars.queue_name
  else:
    s += '''\
#
#
# queue name not specified by the pamq.py user (command incative), using queue system default 
#
##PBS -q %s ''' % GlobalVars.queue_name+str('\n')
#
#  if time is specified in the input, add it
#
  if GlobalVars.queue_time:
    s += '''\
#
#
# specify the walltime 
#
#PBS -l walltime=%s'''%GlobalVars.queue_time+str('\n')
  else:
    s+='''\
#
# queue time not specified by the user, takes system default values in 
#
##PBS -l walltime=%s'''%GlobalVars.queue_time+str('\n')
# ...continue ...
  s += '''\
#
# unite log outputs stdout with stderr
#
#PBS -j oe
#
'''
  if GlobalVars.is_own_queue_stdout:
    s += '''\
#
#   specify queue job standard output - stdout
#
#PBS -o %s''' % GlobalVars.job_out+str('\n')

  s += '''\
#
# employ already set enviroment variables
#
#PBS -V
#
# ... print out important PBS variables into standard output
#
echo '------------------------------------'
echo '  *** PBS control variables ***  ' 
echo '------------------------------------'
echo '   $PBS_O_HOME='$PBS_O_HOME
echo '   $PBS_O_HOST='$PBS_O_HOST
echo '   $PBS_SERVER='$PBS_SERVER
echo '  $PBS_O_QUEUE='$PBS_O_QUEUE
echo '$PBS_O_WORKDIR='$PBS_O_WORKDIR
echo '  $PBS_JOBNAME='$PBS_JOBNAME
echo '    $PBS_JOBID='$PBS_JOBID
echo ' $PBS_NODEFILE='$PBS_NODEFILE
echo '------------------------------------'
#
#
echo ''
echo '****  assigned nodes (content of the $PBS_NODEFILE): ****'
cat $PBS_NODEFILE
#
# ... for MPICH2 library boot hosts, DIRAC must use -machinefile <file> before -np <n_procs>
#
#mpdboot -f $PBS_NODEFILE
#
NUM_OF_CPUS=$(cat /proc/cpuinfo | grep processor | wc -l)
echo -e  'Number of all processors available per core, NUM_OF_CPUS='$NUM_OF_CPUS
NUM_OF_PROCESSES=$(cat ${PBS_NODEFILE} | wc -l)
echo -e 'Number of cores requested, NUM_OF_PROCESSES='$NUM_OF_PROCESSES
# parameter - number of all cores
let NUM_ALL_CORES="$NUM_OF_CPUS*$NUM_OF_PROCESSES"
echo -e 'Number of all processors at disposal, NUM_ALL_CORES='$NUM_ALL_CORES

# ... always change to the home directory 
#  from where you execute DIRAC pam - this is PBS variable 
#
cd $PBS_O_WORKDIR
#
#   ... execute own DIRAC pam script
#
echo 'command to be executed'
echo " %s "
#
''' % GlobalVars.batch_script_command
#
  s += '''\
%s 
''' % GlobalVars.batch_script_command
  s += '''\
#
#
# ... that's whole PBS script, exit
#
exit
'''
# ... write all string to the file
  f.write(s)
  f.close()

 # submit job ...
  if GlobalVars.submit_job:
    if GlobalVars.iprint >= 0 :
      print '\ngoing to submit job with prepared script:',GlobalVars.queue_script
      print '            standard output file will be:',GlobalVars.job_out
      print '                  queue job name will be:',GlobalVars.job_name,'\n'

    os.system("qsub "+GlobalVars.queue_script) # submitting !
    # ... some info ...
    os.system("echo '  -  Current user-status on queueing-system : '")
    os.system("echo 'Job id           Name             User             Time Use S Queue'")
    os.system("echo '---------------- ---------------- ---------------- -------- - -----'")
    os.system("qstat | grep $USER")
  else:
    print '\n not submitting,prepared only the script for the batch system:',GlobalVars.queue_script
    print ''

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_grid_engine_script():
#
#
#   prepare the file of Grid Engine queueing system script and submit it (if specified)
#
#
  global GlobalVars

  if GlobalVars.iprint >8 :
     print '\nin prepare_grid_engine_script...'

  if GlobalVars.iprint >8 :
    print ' home directory:',GlobalVars.home_dir

  f=open(GlobalVars.queue_script,'w') # open file with queue commands

  s = ''
  s += '''\
#!/bin/bash
#
#         Queue controls are directly readable from the submitting batch script
#----------------------------------------------------------------------------------------
#
#     The job is located in the current directory
#
# NOTE: check if woking nodes share this home directory !
# (on power.tau.ac.il not)
#
#$ -cwd
#
#  Shell script  
#NOTE: Check if the script keeps all envormental variables for nodes-spce
# (on power.tau.ac.il not)
#
#
#  on some settings (like power.tau.ac.il not) consider 
# that default batch memory limits might not sufficient for DIRAC  !
#
##$ -l h_vmem=3000M
##$ -l h_vmem=15G  
#
##$ -l h_data=3000M
##$ -l h_data=15G      
#
#
##$ -S /bin/bash
#
#     Job name
#
#$ -N %s
  ''' % GlobalVars.job_name
  s += '''\
#
# number nodes (NOTE: different from number of processors specified by -mpi pam flag !)
#
#$ -pe mpich %i
''' % GlobalVars.nnodes
#
  s += '''\
#
#
# Join both output and standard error streams
#$ -j y
#
# ... queue job standard output - stdout
#$ -o %s
''' % GlobalVars.job_out
#
  s += '''\
#
# Export these environmental variables (for quantanamera.u-strabg.fr)
#
#$ -v P4_GLOBMEMSIZE=171966464
#
#
# ... print out important qeueu variables into standard output
#
echo '-----------------------------------------'
echo '  *** Grid Engine control variables ***  ' 
echo '-----------------------------------------'
echo '       $JOB_NAME='$JOB_NAME
echo '         $JOB_ID='$JOB_ID
echo '       $HOSTNAME='$HOSTNAME
echo '    $SGE_TASK_ID='$SGE_TASK_ID
echo '  $SGE_O_WORKDIR='$SGE_O_WORKDIR
echo '     $SGE_O_HOME='$SGE_O_HOME
echo '  $SGE_O_LOGNAME='$SGE_O_LOGNAME
echo '     $SGE_O_PATH='$SGE_O_PATH
echo '    $SGE_O_SHELL='$SGE_O_SHELL
echo '$SGE_STDOUT_PATH='$SGE_STDOUT_PATH
echo '         $NSLOTS='$NSLOTS
echo '         $NHOSTS='$NHOSTS
echo '        $NQUEUES='$NQUEUES
echo '          $QUEUE='$QUEUE
echo '------------------------------------'
#
#
'''
  s += '''\
#
#   ... execute own DIRAC pam script
#
%s 
''' % GlobalVars.batch_script_command
  s += '''\
#
#
# ... that's whole script, exit
exit
'''
# ... write all string to the file
  f.write(s)
  f.close()

# submit job ... TODO: perhaps it is possible to make it universal ?
  if GlobalVars.submit_job:
    if GlobalVars.iprint >= 0 :
      print '\ngoing to submit job with prepared script:',GlobalVars.queue_script
      print '            standard output file will be:',GlobalVars.job_out
      print '                  queue job name will be:',GlobalVars.job_name,'\n'

    os.system("qsub "+GlobalVars.queue_script) # submitting !
    # ... some info ...
    os.system("echo '  -  Current user-status on queueing-system : '")
    os.system("echo 'Job id           Name             User             Time Use S Queue'")
    os.system("echo '---------------- ---------------- ---------------- -------- - -----'")
    os.system("qstat | grep $USER")
  else:
    print '\n not submitting,prepared only the script for the batch system:',GlobalVars.queue_script
    print ''


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_loadleveler_script():
#
# TODO : add !
#
  global GlobalVars

  if GlobalVars.iprint > 7:
    print 'prepare_loadleveler_script...'

  if GlobalVars.iprint >8 :
    print ' home directory:',GlobalVars.home_dir
   # print ' GlobalVars.pam_script_args=',GlobalVars.pam_script_args

  f=open(GlobalVars.queue_script,'w') # open file with queue commands

  s = ''
  s += '''\
#
# loadl stuff
#
# @ shell            = /bin/sh
# @ network.MPI      = sn_all,not_shared,USnetwork.MPI   
# @ notification     = never
# @ queue      
# 
# 
# @ initialdir       = %s
''' % GlobalVars.home_dir
  s += '''\
#
# @ output           = %s
''' % GlobalVars.job_out
  s += '''\
#
# @ error            = %s
''' % GlobalVars.job_out
  s += '''\
#
# @ node             = %s
''' % str(GlobalVars.nnodes)
  s += '''\
# @ tasks_per_node   = %s
''' % str(GlobalVars.proc_per_node)
  s += '''\
#
# @ wall_clock_limit = %s
''' % GlobalVars.queue_time
  s += '''\
#
# @ environment      = LLSHELL=%s ;ENVIRONMENT=BATCH;
''' % os.getenv("SHELL")
  s += '''\
#
#  change to original directory
#
cd %s
''' % GlobalVars.home_dir
  s += '''\
#
export XLFRTEOPTS="namelist=old"
#
#
# ... execute own DIRAC pam script
#
%s 
''' % GlobalVars.batch_script_command
  s += '''\
#
#
# ... that's whole script, exit
exit
''' 

# ... write all string to the file
  f.write(s)
  f.close()

 # submit job ...
  if GlobalVars.submit_job:
    if GlobalVars.iprint >= 0 :
      print '\ngoing to submit job with prepared script:',GlobalVars.queue_script
      print '            standard output file will be:',GlobalVars.job_out
      print '                  queue job name will be:',GlobalVars.job_name,'\n'

    os.system("llsubmit "+GlobalVars.queue_script) # submitting !
    # ... some info ...
 #   os.system("echo '  -  Current user-status on queueing-system : '")
 #   os.system("echo 'Job id           Name             User             Time Use S Queue'")
 #   os.system("echo '---------------- ---------------- ---------------- -------- - -----'")
 #   os.system("qstat | grep $USER")
  else:
    print '\n not submitting,prepared only the script for the batch system:',GlobalVars.queue_script
    print ''


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_bsub_script():
#
# TODO : add !
#
  if GlobalVars.iprint > 7:
    print 'preparing bsub script...'

  sys.exit()

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_queue_script():
#
# prepare the script according to queue syustem and submit it
#
  global GlobalVars

  if GlobalVars.iprint > 7:
    print 'Going to prepare script for batch system...'

  if GlobalVars.queue_system=='pbs':
    process_pbs_script()	
  elif GlobalVars.queue_system=='loadleveler': 
    process_loadleveler_script()
  elif GlobalVars.queue_system=='sge': 
    process_grid_engine_script()
  elif GlobalVars.queue_system=='bsub': 
    process_bsub_script()
  elif GlobalVars.queue_system=='xgrid': 
    print
#	process_xgrid_script()	
  else:
   print '\nprocess_queue_script: NOT defined queue_system, GlobalVars.queue_system=',GlobalVars.queue_system
   print 'process_queue_script: list of queue systems:\n',GlobalVars.implemented_queue_systems
   sys.exit()


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pamq_arguments():
#----------------------------------------------------------------------------------------
#
#   process script arguments 
#
#----------------------------------------------------------------------------------------
  global GlobalVars
  from optparse import OptionParser

  if GlobalVars.read_from_config_file :
    usage = "usage: %prog arg1 arg2...; defaults read from file:"+GlobalVars.pamq_defaults_file
  else :
    usage = "usage: %prog arg1 arg2...; no defaults get from config file."

  usage=usage+'\n   script dwells in :'+GlobalVars.dirac_home   
  parser = OptionParser(usage)

  #  take care  GlobalVars.pam_script
  if (len(GlobalVars.pam_script)==0) :
    pam_script_add="no default location of DIRAC pam script - YOU MUST SPECIFY IT !"
  else:
    pam_script_add="location of the DIRAC pam script, default:"+GlobalVars.pam_script
  parser.add_option("-p", "--pam", type="string", nargs=1,
                  dest="dest_pam",default=GlobalVars.pam_script,
                  help=pam_script_add, metavar="STRING")

  #  fill GlobalVars.pam_arguments
  dest_pam_arguments=" pam arguments (like '"' --noarch --mol=[*.mol] --inp=[*.inp]'"')"
  parser.add_option("--pam_args", type="string", dest="dest_pam_args",
                  help=dest_pam_arguments, metavar="'"'arg1 arg2 ... '"' ")

  #  fill GlobalVars.testrun_arguments
  dest_test_arguments="'"'testlast'"' arguments (like '"' --tests=short'"')"
  parser.add_option("--test_args", type="string", dest="dest_test_args",
                  help=dest_test_arguments, metavar="'"'args...'"' ")

  # ... get in the print level 
  iprint_label='script verbose print level, default:'+str(GlobalVars.iprint)
  parser.add_option("-i", "--iprint", type="int", nargs=1,  
                  action="store", dest="dest_iprint", default=GlobalVars.iprint, 
                  help=iprint_label, metavar="INT_NUM")

  # ... do not submit flag
  bool_submit=str(GlobalVars.submit_job)
  nosub_label='do not submit job, prepare only the queue script for batch, default-submit:'+bool_submit
  parser.add_option("--nosub", action="store_false", default=GlobalVars.submit_job,
                     dest="dest_nosub", help=nosub_label) 

  # ... specify as serial (no mpi)
  bool_serial=str(GlobalVars.is_serial)
  serial_label='submit job as serial (no mpi), default:'+bool_serial
  parser.add_option("-s","--serial", action="store_true",default=GlobalVars.is_serial,
                     dest="dest_is_serial", help=serial_label) 

  # ... no concrete stdout, use the stdout+stderr system default ... #
  queue_stdout_label='Use a queue system default stdout+stderr log file,not the own file '+GlobalVars.job_out
  parser.add_option("--nso",action="store_false",default=GlobalVars.is_own_queue_stdout,
                      dest="dest_is_own_queue_stdout",help=queue_stdout_label)

  # ... get number of nodes
  if (int(GlobalVars.nnodes) <=0) :
    nodes_label='number of reserved nodes, default:'+str(GlobalVars.nnodes)+'=> YOU MUST SET THE NUMBER !'
  else :
    nodes_label='number of reserved nodes, default:'+str(GlobalVars.nnodes)
  parser.add_option("-n", "--nodes", type="int",  nargs=1,
                  action="store",dest="dest_nnodes", default=GlobalVars.nnodes, 
                  help=nodes_label, metavar="INT_NUM")

  # ... GlobalVars.queue_system
  if (len(GlobalVars.queue_system)>0) :
    queue_system_label='specify the queuing system, default :'+"'"+GlobalVars.queue_system+"'"
    parser.add_option("--qsys", type="string", nargs=1,
                  action="store",dest="dest_queue_system", default=GlobalVars.queue_system, 
                  help=queue_system_label, metavar="STRING")
  else :
    queue_system_label='no default queuing system, YOU MUST SPECIFY IT !'
    parser.add_option("--qsys",type="string",nargs=1,action="store",dest="dest_queue_system",help=queue_system_label,metavar="STRING")

  # ... GlobalVars.job_name
  queue_job_name_label='specify the name for the queue run, default:'+"'"+GlobalVars.job_name+"'"
  parser.add_option("--name", type="string", nargs=1,
                  action="store",dest="dest_queue_job_name",  
                  help=queue_job_name_label,metavar="STRING")

  # ... GlobalVars.queue_name
  if (GlobalVars.queue_name_specified) :
    queue_label='specify the queue name, default :'+GlobalVars.queue_name
    parser.add_option("-q", "--queue", type="string", nargs=1,
                  action="store", dest="dest_queue_name",default=GlobalVars.queue_name,
                  help=queue_label,metavar="STRING")
  else :
    queue_label='specify the queue (otherwise uses system default)'
    parser.add_option("-q", "--queue", type="string", nargs=1,
                  action="store", dest="dest_queue_name", 
                  help=queue_label,metavar="STRING")

  # ... GlobalVars.proc_per_node ... "--ppn" not accepting !
  if GlobalVars.proc_per_node==None:
    queue_ppn_label='specify number of processor per node, default:'+str(GlobalVars.proc_per_node)+'<= takes all CPUs per core'
  else:
    queue_ppn_label='specify number of processor per node, default:'+str(GlobalVars.proc_per_node)

  parser.add_option("--ppn", type="int", nargs=1, 
                   action="store", dest="dest_ppn", default=GlobalVars.proc_per_node,
                  help=queue_ppn_label, metavar="INT")

  # ... GlobalVars.queue_time
  queue_time_label='specify the time in the queue (default set by the system)'
  parser.add_option("-t", "--qtime", type="string", nargs=1,
                  action="store", dest="dest_queue_time", 
                  help=queue_time_label,metavar="STRING")

  # ... prepare file with selected default variables (configuration file) ...
  config_file_label='save read parameters (like qsys, queue... ) as defaults into file ' \
               +GlobalVars.dirac_home+'/PAMQ_DEFAULTS; default - do not save'
  parser.add_option("--config", action="store_true", default=GlobalVars.prepare_config_file,
                     dest="dest_prep_config", help=config_file_label) 

  ###################################################################################################
  #
  #                                    ... finally assign
  #
  ###################################################################################################
  (options, args) = parser.parse_args()

# .. now assign individual arguments 

  if GlobalVars.iprint > 5:
    print '\n processing arguments:  len(sys.argv)=', len(sys.argv) 

  # ... handle queue job name
  if (options.dest_queue_job_name is not None) :
    GlobalVars.job_name=options.dest_queue_job_name
    #print '* assigned GlobalVars.job_name=',GlobalVars.job_name
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.job_name=',GlobalVars.job_name

  # ... handle iprint
  if (options.dest_iprint is not None) :
    GlobalVars.iprint=int(options.dest_iprint)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.iprint=',GlobalVars.iprint

  # number of nodes
  if (options.dest_nnodes is not None) :
    GlobalVars.nnodes=int(options.dest_nnodes)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.nnodes=',GlobalVars.nnodes

  # about submitting job
  if (options.dest_nosub is not None) :
    GlobalVars.submit_job=options.dest_nosub
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.submit_job=',GlobalVars.submit_job

  #  .... no mpi stuff for pam 
  if (options.dest_is_serial is not None) :
    GlobalVars.is_serial=options.dest_is_serial
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.is_serial=',GlobalVars.is_serial

  #  ... take care of your own stdout+stderr log file ...
  if (options.dest_is_own_queue_stdout is not None) :
    GlobalVars.is_own_queue_stdout=options.dest_is_own_queue_stdout
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.is_own_queue_stdout=',GlobalVars.is_own_queue_stdout
  else :
    GlobalVars.is_own_queue_stdout=True
    if GlobalVars.iprint >= 7:
      print 'assigned boolean GlobalVars.is_own_queue_stdout=',GlobalVars.is_own_queue_stdout

  # ... save the pam directory
  if (options.dest_pam is not None) :
   # GlobalVars.pam_script=options.dest_pam
    GlobalVars.save_pam_script(options.dest_pam)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.pam_script=',GlobalVars.pam_script

  # ... deal with pam arguments
  if (options.dest_pam_args is not None) :
    GlobalVars.pam_arguments=options.dest_pam_args
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.pam_arguments=',GlobalVars.pam_arguments
  else :
    GlobalVars.submit_job=False

  # ... deal with testlast arguments
  if (options.dest_test_args is not None) :
    GlobalVars.testrun_arguments=options.dest_test_args
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.testrun_arguments=',GlobalVars.testrun_arguments
    GlobalVars.dirac_tests_run=True
  else :
    GlobalVars.dirac_tests_run=False

  # ... deal with queue_name
  if (options.dest_queue_name is not None) :
    GlobalVars.queue_name=options.dest_queue_name
    GlobalVars.queue_name_specified=True
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_name=',GlobalVars.queue_name
      print 'assigned GlobalVars.queue_name_specified=',GlobalVars.queue_name_specified

  # ... deal with the queue_system
  if (options.dest_queue_system is not None) :
    GlobalVars.queue_system=options.dest_queue_system
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_system=',GlobalVars.queue_system

  # ... deal with  number of processors per 1 cluster node
 # print '* options.dest_ppn=',options.dest_ppn
  if (options.dest_ppn is not None) :
  #  print 'options.dest_ppn=',options.dest_ppn
    GlobalVars.proc_per_node=int(options.dest_ppn)
    if GlobalVars.iprint >= 7:
      print 'assigned GlobalVars.proc_per_node=',GlobalVars.proc_per_node

  # ... deal with queue_time
  if (options.dest_queue_time is not None) :
    GlobalVars.queue_time=options.dest_queue_time
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_time=',GlobalVars.queue_name

  # ... deal with creating variables_defaults file
  if (options.dest_prep_config is not None) :
    GlobalVars.prepare_config_file=options.dest_prep_config
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.prepare_config_file=',GlobalVars.prepare_config_file

  # ... must be at leat one argument - pam_args !
  #if len(args) < 1:
  #  parser.error("incorrect number of arguments")

  if (len(parser.largs)>0) :
    print 'error: leftover arguments found :',parser.largs 
    parser.error("incorrect passing arguments, all must be under flags !")

  if GlobalVars.iprint > 8 :
    print '\n processing arguments' 
    print 'len(sys.argv)=', len(sys.argv)
    for i in range(len(sys.argv)) :
      print i,sys.argv[i]
		
  #------------------------------------------------------------------------
  #      collect all left pamq.py argumnents into  pam_arguments
  #------------------------------------------------------------------------
  for item in parser.largs:
    GlobalVars.pam_arguments= GlobalVars.pam_arguments+' '+item
  GlobalVars.pam_arguments=' '+GlobalVars.pam_arguments

  if GlobalVars.iprint > 7:
    print 'extracted pam arguments:',GlobalVars.pam_arguments

  # set mpi for pam 
  if not GlobalVars.is_serial :
    GlobalVars.set_pam_mpi()

  # if prepare file, do not submit
  if GlobalVars.prepare_config_file:
    GlobalVars.submit_job=False
    GlobalVars.create_configure_file()
    if GlobalVars.iprint > 7:
      print 'preparing configure file-not submitting job '

  # ... prepare final script command
  GlobalVars.prepare_batch_script_command()
  
#########################################################################
#
#                        MAIN FUNCTION
#
# MI: I had to remove error catching to be able to localize bugs in the script
#
#########################################################################
def main(*args):
#  try:
    global GlobalVars
    GlobalVars = PamQVariables() 

  #  if len(sys.argv) < 2 :
  #    pamq_usage()
  #    print 'error exit !'
  #    sys.exit(2)

    process_pamq_arguments()

    if (not GlobalVars.prepare_config_file) :
      process_queue_script()

#  except:
    # handle some exceptions
#    print 'caught exception...'
#  else:
#    return 0 # exit errorlessly
 
###########################################################################
###                   process the main function
###########################################################################
if __name__ == '__main__':
    sys.exit(main(*sys.argv))
