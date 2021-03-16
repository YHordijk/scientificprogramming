:orphan:


Release preparation
===================


Feature freeze
--------------

Two months before the release (example: DIRAC16), a release branch (example:
release-16) is created from master. This is a feature freeze and from this
moment on the release branch ideally only receives cosmetics and bugfixes, no
major new features, no merges from master or other branches (except
cherry-picks; more about it later). You as developer are responsible to commit
or merge the to-be-released code to master, before the release branch is
created.


Can I commit bleeding edge code to master after the feature freeze?
-------------------------------------------------------------------

Yes you can! This is one of the reasons we have the release branch and exactly for
this reason we will never merge from master to the release branch, only from
the release branch to master.


How to check out the release branch (example release-16 for DIRAC16)
--------------------------------------------------------------------

Check out the release branch::

  $ git checkout release-16

Now you have it next to master (verify this)::

  $ git branch

You can switch back to master::

  $ git checkout master

And then switch back to the release branch::

  $ git checkout release-16


How to commit changes and bugfixes that are relevant for the released code
--------------------------------------------------------------------------

Commit all such changes to the release branch, not to master. This way no
commit or bugfix will get lost. Please do *not* transfer commits from the
release branch to master manually. This is not only unnecessary work but it is
harmful (conflicts)!
Again, do not commit to master, commit to the release branch::

  $ git checkout release-16
  $ git pull                      # update the local release-16 branch
  $ git commit                    # commit your modifications
  $ git push                      # push your changes to release-16 branch
  $ git checkout master           # switch to the master branch
  $ git merge --no-ff release-16  # merge changes to master, fix conflicts if any
  $ git push                      # push your changes to master

Note that by merging your changes to master you might get conflicts. Git points
them out clearly. In such cases open conflicting file(s) and fix discrepancies
inside manually. They are marked with "<<<<<<" and ">>>>>>" strings.

What if you accidentally committed something to master which belongs to
the release? In this case do *not* merge master to the release branch, but
rather cherry-pick the individual commit(s) to the release branch::

  $ git checkout release-16
  $ git cherry-pick [hash]        # with [hash] that corresponds
                                  # to the commit


How to exclude code from being released
---------------------------------------

You as the developer are responsible for removing your code from the
release branch that you do not wish to be released. This is *not* done by
actually removing the code but by using CMake in combination with CPP
statements. So again, *do not actually delete* code lines, use CPP statements
as demonstrated below:

Modules can be excluded from the release version by surrounding the appropriate
parts of code by preprocessor statements, e.g.::

  #ifdef MOD_KRMC
           CALL TKRMCORB
           CALL SETDC2(0)
           CALL PSIOPT('MCSCF',WFCONV,WORK,LWORK)
           IF(DOANA) CALL PAMANA(WORK,LWORK)
  #else
              CALL QUIT(
       &  'Second order KMCSCF optimization not included in this version')
  #endif

This means that the calls are included only if MOD_KRMC is defined.
In the above example, the
script would remove the actual MOD_KRMC code and what would be left is::

  CALL QUIT('Second order KMCSCF optimization not included in this version')

If you want to protect code from being released which is not part of any CMake
option, you can use MOD_UNRELEASED::

  #ifdef MOD_UNRELEASED
  !      secret code; this will not end up in the release tarball
  #endif

To summarize, you can remove code using either existing CPP filters or using
MOD_UNRELEASED. Do this on the release branch and merge to the master.


How does the preprocessor-based code removal work?
--------------------------------------------------

Consider this CMake code at the time of writing found in ``src/CMakeLists.txt``::

    cmake_dependent_option(ENABLE_INTEREST "Enable Interest" ON  "ENABLE_UNRELEASED" OFF)
    if(ENABLE_INTEREST)
        add_definitions(-DMOD_INTEREST)
    endif()

    ...

    if(ENABLE_INTEREST)
        add_subdirectory(${PROJECT_SOURCE_DIR}/src/interest)
        set(EXTERNAL_OBJECTS $<TARGET_OBJECTS:interest> ${EXTERNAL_OBJECTS})
        dirac_test(mdirac_closed-shell_dft "interest")
    else()
        set(CODE_REMOVAL_FLAGS ${CODE_REMOVAL_FLAGS} -UMOD_INTEREST)
    endif()
    message(STATUS "Interest library: ${ENABLE_INTEREST}")

This code adds a CMake option ``ENABLE_INTEREST`` to toggle the Interest
library on or off.  This is a dependent option which depends on
``ENABLE_UNRELEASED``. This means that Interest library is on by default unless
``ENABLE_UNRELEASED`` is selected.

We see that selecting ``ENABLE_INTEREST`` defines ``-DMOD_INTEREST`` (among
other things).  Further down we see that if ``ENABLE_INTEREST`` is false, then
``-UMOD_INTEREST`` is added to ``CODE_REMOVAL_FLAGS``. This variable is used in
``cmake/custom/cpack.cmake`` and passed on to the
``maintenance/remove_unreleased_code`` script which removes code based on these
flags.  Since ``maintenance/remove_unreleased_code`` works in-place,
``cmake/custom/cpack.cmake`` first copies the source tree to ``EXPORT_DIR`` to
make sure we don't delete code in our working directory.


How to "remove" entire directories?
-----------------------------------

See ``cmake/custom/cpack.cmake`` and search for ``rm -rf``.


How to "remove" tests?
----------------------

Currently you cannot. If you really want to remove it, list it explicitly
in ``cmake/custom/cpack.cmake``.


How to release or export the code
---------------------------------

Most important rule: treat unpublished code of others with respect.  Do not
release code of others without asking.

1. Verify copyright headers.
2. Verify logo (list of authors).
3. Bump the version (the version number is kept in file VERSION; but also git grep -i for "orphaned" version numbers in the source code).

As a final step create the tarball::

  $ mkdir build
  $ cd build
  $ cmake -DENABLE_UNRELEASED=OFF ..
  $ make release


Copyright statements and markers
--------------------------------

In the development source code you can find the following
copyright statement markers::

  !dirac_copyright_start
  !dirac_copyright_end
  !dalton_copyright_start
  !dalton_copyright_end
  \* dirac_copyright_start \*
  \* dirac_copyright_end \*
  \* dalton_copyright_start \*
  \* dalton_copyright_end \*

They are there so that we do not have to update
copyrights manually. You can update them with the
``maintenance/update_copyright.py`` script.
This script will overwrite whatever is between
these markers so it is a bad idea to modify text
between the markers manually.

These markers are removed by "make release"
which runs ``maintenance/remove_copyright_tags`` in the background
(see ``cmake/ConfigPackaging.cmake``).
