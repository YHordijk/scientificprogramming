:orphan:
 

Expert options
==============


Compiling in verbose mode
-------------------------

Sometimes you want to see the actual compiler flags and definitions when compiling the code::

  $ make VERBOSE=1


How can I change optimization flags?
------------------------------------

You can turn optimization off (debug mode) like this::

  $ ./setup --type=debug [other flags]
  $ cd build
  $ make

You can edit compiler flags under the path ``cmake/custom/compiler_flags``-

Alternatively you can edit compiler flags through ccmake::

  $ cd build
  $ ccmake ..
