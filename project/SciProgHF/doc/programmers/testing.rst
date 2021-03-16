:orphan:
 

Never commit functionality to the main development line without tests
---------------------------------------------------------------------

If you commit functionality to the main development line without tests then
this functionality will break sooner or later and we have no automatic
mechanism to detect it. Committing new code without tests is bad karma.


When you add a new test always add it to CMake
----------------------------------------------

Otherwise the test will never be run as part of the default test suite.
``git grep add_test`` to see where tests are defined. For details, see below.


Strive for portability
----------------------

Avoid shell programming or symlinks in test scripts otherwise the tests are not
portable to Windows. Therefore do not use os.system() or os.symlink(). Do not
use explicit forward slashes for paths, instead use os.path.join().


Never add inputs to the test directories which are never run
------------------------------------------------------------

We want all inputs and outputs inside test/ to be accessible by the default test
suite. Otherwise we have no automatic way to detect that some inputs or outputs
have degraded. And degraded inputs and outputs are useless and confusing.


How can I run a single test?
----------------------------

The test runner is CTest. You can run a single test with::

  ctest -R mynewtest

where *mynewtest* is the name of test directory.

Alternatively you can run a test directly in the test directory.
For this check out::

  cd test/mynewtest
  ./test --help


Always test that the test really works
--------------------------------------

It is easy to make a mistake and create a test which is always "successful".
Test that your test catches mistakes. Verify whether it extracts the right
numbers.


How to add labels to tests
--------------------------

Individual tests can have extra descriptors (labels, tags) helping to categorize them.
Pick suitable label names characterizing given test.

Here is an example::

  dirac_test('dft_ac', 'short;dft')

If you give several labels, separate them with ";".

The advantage of labeling a test is that one can call group of tests according to their labels. For instance,
launching all tests containing label "short"::

  ctest -L short


How to add a new test
---------------------

1. First thing to add a new test is to add a python script called "test" in your 
test directory.

In most cases you probably want to use the runtest_dirac.py test library
since it gives you the framework to run DIRAC jobs and
to filter numbers with numerical tolerance.  But you don't have to. You have
full liberty to test whatever you want in whatever way you want. Important is
that you return an exit code. If the exit code is 0, the test was successful,
if non-0, then the test has failed. 

Please provide also README.rst description file to the test directory so that other developers
know what functionalities the test covers, what is the meaning of the test.

2. You have to add the new test also to CMake infrastructure,
otherwise the test will never be run as part of the selected test suite.
Read the following.

Tests categories in Dirac
-------------------------

Tests are grouped into several categories. Consider to choose suitable class for your test
depending on test's execution time, test's intention etc.

2.1. Default tests. These are always run as they cover at least released Dirac functionalities. 
Be sure they do take not too much time !

Modify the :download:`cmake/custom/test.cmake file <../../../cmake/custom/test.cmake>`.


2.2. Tutorial tests. These are added to default tests upon  **-D ENABLE_TUTORIALS=ON**.
Tutorial tests are directly interconnected with tutorial documentation files,
which refer also to corresponding input and output files of given tutorial test(s). 
Therefore these tests serve mainly to 
verify correctness of (ascii) files read and produced by Dirac.
These tests are not restricted by time as are default tests, but please consider rather shorter than longer examples.

See and modify the :download:`cmake/custom/tutorials.cmake file <../../../cmake/custom/tutorials.cmake>`.


2.3 Benchmark tests.
These tests serve to asses the Dirac execution performance.
These tests are activated with *ctest -C benchmarks -L benchmark*.

See and modify the :download:`cmake/custom/benchmarks.cmake file <../../../cmake/custom/benchmarks.cmake>`.


2.4 Unreleased and features restricted tests. 
In the DIRAC development version, unreleased and preprocessor restricted
features are accompanied with dedicated tests. 

How to include them, see
:download:`src/CMakeLists.tx file <../../../src/CMakeLists.txt>` for feature attached tests and the
:download:`cmake/custom/unreleased.cmake <../../../cmake/custom/unreleased.cmake>` file for testing of 
remaining unreleased features.
