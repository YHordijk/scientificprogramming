:orphan:

.. _windows_developers:


DIRAC on Windows - developers section
=====================================

After completing the explanation of DIRAC installation and first steps on the Microsoft (MS) Windows (shortly Windows) platform 
for users (see :ref:`windows_users`)  we continue this manual for the developers.

Motivation
----------

DIRAC is being developed and works mainly on machines with Linux/Unix/Mac operating systems.
Since many desktop computers are sold with preinstalled Windows operating system, it is worth, we believe, to adapt the DIRAC software   
for the Windows environment. Also, the Windows operating system has great potential for the
`high-performance computing <http://www.microsoft.com/education/en-au/solutions/Pages/high-performance-computing.aspx>`_.
This is the motivation for investing efforts in the code's adaptation for the Windows platform.

Thanks to the modern, platform-universal *CMake* buildup system and related Windows supported software 
you can fully set up the DIRAC software under the Windows operating system.

We have tested the DIRAC downloading from its git-repositories, configuration, compilation and tests running, 
all in the framework of the Windows 7/8 operating systems.


Dependencies
------------

In addition to the necessities listed in the section for users, DIRAC developers have to install:

5. The Git version control system (http://git-scm.com). When ask you can choose options "Use Git from Git Bash only" or "Use Git from the Windows Command Prompt" (recommended). If you choose the first option then you have to manually add into your %PATH% *"...\\Git\\cmd"* (depends on your installation). It is not possible to successfully run CMake with optional Unix tools in %PATH%.

6. The Sphinx documentation generator and related packages, see below.

Documentation tools
-------------------

To be able to generate the html documentation (*mingw32-make html*), the user has to download and install these Python packages:

*Sphinx* (http://sphinx-doc.org/). Python Wheel can be found at https://pypi.python.org/pypi/Sphinx/ and then run (for example): ::

  pip install Sphinx-1.3b1-py2.py3-none-any.whl

You might need also the *python-dateutil* and *pyparsing* packages which you can obtain from https://pypi.python.org/pypi/python-dateutil/ and 
https://pypi.python.org/pypi/pyparsing/ respectively. You can install their Python Wheels using *pip install* command.

Next download Python Wheel file of *sphinxcontrib-bibtex* from https://pypi.python.org/pypi/sphinxcontrib-bibtex and Python Wheel file of *sphinx_rtd_theme* from https://pypi.python.org/pypi/sphinx_rtd_theme. Install both using *pip install* command.

*numpy* (http://www.numpy.org/), which provides own (self-)installation file. For 32-bit Python you can download installer from
http://sourceforge.net/projects/numpy/files/NumPy/ and for 64-bit Python there are only unofficial builds for example at http://www.lfd.uci.edu/~gohlke/pythonlibs/.

*matplotlib* (http://matplotlib.org/). Windows installer or Python Wheel can be found at http://matplotlib.org/downloads.html.
Please download version suitable for your version (2.x or 3.x) and type (32/64 bit) of Python.

Install *doxygen* if you want to generate documentation using *mingw32-make doxygen*. You can find self-installer (*.exe*) at http://www.stack.nl/~dimitri/doxygen/download.html. We recommend to use self-installer because it adds doxygen binary directory to system path variable automatically.

To create html slides (*mingw32-make slides*). Install *hieroglyph* (https://pypi.python.org/pypi/hieroglyph) by command: ::

  pip install hieroglyph

Notes:

Check that the Windows *%Path%* environment variable points to both Python and its extensions, and has the the form of *"...C:\\Python27;C:\\Python27\\Scripts;..."*. Doxygen binary directory *""...C:\\DOXYGEN\\BIN;...* should also be present in *%Path%*. Do not use 32 and 64 bit packages together. You can obtain *pip* from https://pip.pypa.io/en/latest/.


Git & PuTTy
-----------

To enable server connectivity please install full PuTTy package for handling your SSH key pair.
The best way is through the PuTTy self installer, see
http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html.

Prepare your own pair of SSH keys - private and public. You can generate them
with the PuTTyGen generator which writes private key in wanted ppk-format.

For downloading/uploading from/to the DIRAC repository you have to have the *PuTTy "Pageant"* authentication agent 
active with loaded private key, while the public key has to be placed in the DIRAC git repository user space.
Pageant can be activated during the Windows startup, for example, launched at the start of user's login on Windows machine: ::

  "C:\Program Files\PuTTY\pageant.exe"  C:\Users\<user_name>\.ssh\id_rsa_Windows7.ppk

The *"pageant"* agent, running in the background, also enables passwordless terminal connections (or PuTTy sessions)
onto server which are provided with users's public key. 
  
For communication with the DIRAC repository activate the PuTTy's *"plink.exe"* executable through the "GIT_SSH" 
environment variable (please do not add this variable into system environment variables because people who do not want to 
use *"pageant"* can have problems to normally use *"git"* - create user environment variable instead). Its typical value is ::

 GIT_SSH=C:\Program Files (x86)\PuTTY\plink.exe

assuming that the PuTTy package is placed in the default installation space: ::

  C:\Program Files (x86).

Also, before using *"git"* refresh your private key in your computer's cache, use ::

  "C:\Program Files (x86)\PuTTY\plink.exe" -agent git@gitlab.com:dirac/dirac.git
  plink -agent git@bitbucket.org

For comfortable commits with the git version control system you should set up the editor.
For the popular *Notepad++* editor setting for git `follow  <http://stackoverflow.com/questions/10564/how-can-i-set-up-an-editor-to-work-with-git-on-windows>`_ :  ::

 git config --global core.editor "'C:/Program Files/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"

To update your local DIRAC source (must have the ssh connection enabled, see above):  ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk>git pull

Do not configure DIRAC in the Linux-like "git-bash" window, which is part of common Git installation.
Executing cmake configuration commands in this windows causes unhomogenous treatment of path-strings: mixing
of Linux and Windows path-delimiters.

Some caveats
------------

Windows operating system does not accept certain reserved words for file and folder names, like aux, con ...
(google after forbidden file and folder names on Windows).

Likewise Windows does not like "*"-characters in file names and due to this is not able to download
such files from the repository.
For that reason several traditional basis set files have been renamed accordingly, see the following table:

==================       ==================
Windows new name         Previous file name
==================       ==================
3-21++G-star             3-21++G*
3-21G-star               3-21G*
6-311+G-star             6-311+G*
6-311G-star              6-311G*
6-311++G-star-star       6-311++G**
6-311G-star-star         6-311G**
6-31++G-star             6-31++G*
6-31+G-star              6-31+G*
6-31G-star               6-31G*
6-31++G-star-star        6-31++G**
6-31G-star-star          6-31G** 
==================       ==================

The solution is that the basis set reader automatically translates "*" basis set names
to the "star" filenames.  

Testing
-------

Run specific tests: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -VV -R fscc
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -V -D ExperimentalTest -R fscc -D ExperimentalSubmit

Perform the complete CDash buildup in one step: ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make Experimental

or ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>ctest -D Experimental

You can perform the CDash buildup in separate steps:  ::

 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalUpdate
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalConfigure
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalBuild
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalTest
 C:\Users\<username>\Documents\Work\Dirac_git\trunk\build>mingw32-make ExperimentalSubmit


Automatic test runs 
-------------------

You have to set up the *Task Scheduler* accordingly for DIRAC automatic buildups. 

Don't forget that in order to perform the git-synchronization you must have your private key loaded 
in the *pageant* program (see above). 
This can be achieved in this way (bar Actions of the Task Scheduler): :: 

  Action: Start a program

  Program/script:
  "C:\Program Files (x86)\PuTTY\pageant.exe"

  Add arguments (optional):
  "C:\Users\milias\Documents\milias_private_key.ppk" -c  "C:\Users\milias\Documents\Dirac\software\devel_trunk\maintenance\cdash\cdash.Win8.0.FPV-Miro.bat"

In the buildup script (here cdash.Win8.0.FPV-Miro.bat) - in order to kill the zombie pageant process - you may to choose either to kill 
the thread directly or to restart your Windows PC. In the latter case please check the sleeping mode is working not have the office 
PC running in the night. The buildup script wakes the PC up, which after the whole build switches back to sleeping.


Further development plans
-------------------------

So far we have tested DIRAC in sequential buildups (only) with the MinGW64 (free of charge) suite of compilers.
To profit from multithreading under Windows, enthusiastic developers (where are you?) should investigate Open MPI parallelization, www.open-mpi.org. 
Also, one should try some commercial compilers, like Intel, PGI...

Concerning mathematical libraries, we have to resort to DIRAC own BLAS & LAPACK libraries. 
The ATLAS library preparation for Windows is only in the experimental phase (http://math-atlas.sourceforge.net/errata.html#win64).
There is the Intel-MKL library at hand, but costs some money ;). For sure we should investigate the Windows version of the ACML library.
OpenBlas http://www.openblas.net/.

Another topic is the flavour of various CMake generators. On Linux systems we are using the default choice, "Unix Makefiles". 
On the MS Windows operating system we are employing the "MinGW Makefiles" generator - part of the GNU MinGW suite. However, there are many more generators worth to try, 
like NMake, MSYS Makefile, Visual Studio, Codeblocks...
