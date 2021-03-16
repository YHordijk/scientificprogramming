:orphan:
 

F.A.Q.
======

How can I get the AO density matrix?
------------------------------------

Use subroutine DENMAT dirac/dirden.F.


How can I get MO coefficients?
------------------------------

Use subroutine REACMO_NO_WORK
in dirac/dirgp.F.


How are density/Fock/... matrices blocked?
------------------------------------------

They are blocked according to NBBAS and IBBAS.


Why do all routines end with "RETURN" before "END"?
---------------------------------------------------

Most of the routines end with::

  subroutine foo()
  ...
  return
  end

The "return" is redundant. Simply write "end" without "return".  To make it
better, we recommend programmers to stick to the Fortran90/95 coding standards.
So the structure of the routine shall be as::

  subroutine foo()
  ...
  end subroutine


What are all the "Decks" before the subroutines good for?
---------------------------------------------------------

Historical reasons.
They have no meaning today.


Can I use zgemm or do I have to use the module matrix_defop?
------------------------------------------------------------

There is nothing wrong with zgemm. You can use either one - whichever you prefer.


I want my new module to be completely modular. Can I then use DIRAC infrastructure routines like QUIT or READT?
---------------------------------------------------------------------------------------------------------------

If you want to be 100% independent of the DIRAC infrastructure then you cannot
use such routines.
