:orphan:
 

star(MULPOP)

Mulliken population analysis is performed in AO basis. The analysis is
based on the concept of *labels*. Each basis function is labeled by its
functional type and center (and boson irrep when in SO basis, the
default). The labels are given in the output. A set of primitive labels
may be collected to group labels as specified by the user.


keyword(VECPOP)

For each fermion irrep, give an :ref:`orbital_strings` of orbitals to analyze.

*Default:* Analyze the occupied electronic solutions.


keyword(NETPOP)

Give Mulliken net/overlap populations in addition to gross populations.

keyword(LABEL)

Use pre-defined labels for use in Mulliken population analysis.
The following defintions are recognized by DIRAC:

::

    .LABEL
    ATOM

Provide labels for individual atoms. This is useful for getting Mulliken charges.

::

    .LABEL
    SHELL

Provide labels for individual orbital types.

keyword(LABDEF)

Define labels for use in Mulliken population analysis.

Give first the number of labels to define, then the label name (12 characters)
and :ref:`orbital_strings` for each fermion irrep.

*Example:*

::

    .LABDEF
    4 (total number of lines to follow)
    Re s        1
    Re p        2..4
    Re d        5..10
    Re f        11..20
     ...        21..35

*Pitfall:*

If you decide to redefine labels for the population analysis you have to
reassign **all** of them, otherwise Dirac is unable to continue. The
list of functions with the initial label definitions (*default:
SO-labels*) can be drawn from the DIRAC output starting at the section :

::

                                   GETLAB: SO-labels

       Large components:   35
        1  L Ag AR s        2  L Ag AR dxx      3  L Ag AR dyy      4  L Ag AR dzz      etc.


keyword(AOLAB)

Base definition of labels on AO basis.

*Default:* Base definition of labels on SO basis.


keyword(INDSML)

Use individual small component labels.

*Default:* Gather all small component functions belonging to a given
center and irrep, irrespective of function type.


keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

