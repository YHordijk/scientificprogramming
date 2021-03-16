C     File: dcbpop.h
C
C     Local common block for Mulliken Population analysis in dirana.F
C     (set in routine POPINP in dirrdn.F)
C
C     POPLAB - labels for population analysis
C     VECPOP(2) - orbital strings for population analysis
C     VECPOP_SAVE(2) - for backup copy of user input VECPOP
C
C     ADD_ALL   : collect large and small components for a given center
C                 into one label (AO-def): for DO_SCFPOP
C     ADDSML    : collect small components for a given center into one 
C                 label (AO-def)
C                 collect small components for given center combination
C                 into one label (SO-def)
C     DONETP    : also print net populations
C     LABDEF    : use user defined labels in POPLAB(*)
C     DO_SCFPOP : do atomic gross mull. populations for SCF iterations
C
C     IPRPOP - print level
C     ILABDF = 1: base label definitions on primitive AO-labels
C            = 2: base label definitions on primitive SO-labels
C     NPOPLAB   : number of labels for population analysis
C     IPOPLAB(1:PNLAB(0)) - pointer from AO/SO-label to POPLAB
C     ICHRG(1:NPOPLAB) - nuclear charge of each entry in POPLAB (for DO_SCFPOP)
C
#include <dcblab.h>
      CHARACTER POPLAB*12,VECPOP*72, VECPOP_SAVE*72
      COMMON/CBCPOP/POPLAB(MAXLAB),VECPOP(2),VECPOP_SAVE(2)
      LOGICAL ADD_ALL,ADDSML,DONETP,LABDEF, DO_SCFPOP
      COMMON/CBLPOP/ADD_ALL,ADDSML,DONETP,LABDEF,DO_SCFPOP
      COMMON/CBIPOP/IPRPOP,ILABDF,NPOPLAB,IPOPLAB(MAXLAB)
      REAL*8  CHRGNUC
      COMMON/CBRPOP/CHRGNUC(MAXLAB)
C -- end of dcbpop.h --
