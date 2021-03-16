!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DISKUN
      IMPLICIT REAL*8           (A-H,O-Z)
*
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      INTEGER   ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "files.inc"
#include "clunit.inc"
      INTEGER LUSCR_MPIL,LUDIA_MPIL,LUC_MPIL, LUHC_MPIL, LUSC123_MPIL
      INTEGER LUSC344_MPIL, LUDMMY_MPIL
      CHARACTER LUSCR_MPIN*11,NLUSCR_MPI*15,LUDMMY_MPIN*5, NLUDMMY_MPI*9
*
* Assign logical unit numbers for LUCIA:
*==================================================
*  Complete revision motivated by MCSCF development
*==================================================
*
*   1 -  9 :  Retained for I/O
*  ---------
*        5 :  Standard IN
*        6 :  Standard OUT
*
*  60 - 69 :  Important files retained (DIAG,C,HC etc.)
*  ---------
*       60 :    CI diagonal
*       61 :    CI reference vector
*       62 :    CI sigma vector
*
*  70 - 79 :  Scratch files
*  ---------
*
*  80 - 99 :  Files in use in different context
*             (e.g. standalone, INS, MOLCAS etc.)
*  ---------
*
*  Timo Fleig    Mar 28, 2002
*==================================================
*
*  Old setup:
*
* All file with some kind of input information  :  10 - 19
* All files containing final results            :  90 - 99
* Scratch files                                 :  30 - 40
* Internal files (retained through job)         :  20 - 29
*
* =========================
* Standard input and output
* =========================
*. Input file
      LUIN = 5
*. Special unit for input file from MOLCAS conversion routine
*    Timo Fleig
      LUCMOL = 91
*. Output file
      IF (MYPROC .EQ. MASTER) THEN
      LUOUT = 6
      END IF
* =================
* Input information
* =================
* Input file containing MO-AO transformation matrix
      LUMOIN = 80
*. Input file for CI-vectors
*. restart from file 21 is assumed
*. Input , two electron integrals - MOLCAS
      LUTINT = 81
*. Input , one electron integrals - MOLCAS
      LUOINT = 82
*. Input , property one-electron integral files
      LUPRP  = 83
*. Sirius interface file
      LUSIR1 = 84
*. File containing additional states for transition densities
      LUEXC = 85
*. Scratch file for master <--> slave communication
      LUSCR_MPI = 51
      LUSCR_MPIN ="SIGMA_PAR" 
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSCR_MPI,'(A9,A1,I1)') LUSCR_MPIN,'.',MYPROC
         LUSCR_MPIL=13
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSCR_MPI,'(A9,A1,I2)') LUSCR_MPIN,'.',MYPROC
         LUSCR_MPIL=14
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSCR_MPI,'(A9,A1,I3)') LUSCR_MPIN,'.',MYPROC
         LUSCR_MPIL=15
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCI*_NODE MODULE")
      ENDIF
#else
      NLUSCR_MPI=LUSCR_MPIN
      LUSCR_MPIL=13
#endif
      Open(Unit=LUSCR_MPI,File=NLUSCR_MPI(1:LUSCR_MPIL),
     &     Status='UNKNOWN',
     &     Form='UNFORMATTED')

* =================
* Internal files
* =================
*. CI diagonal
      LUDIA = 60
*
      LUDIA_MPIN ="CIDIA"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUDIA_MPI,'(A5,A1,I1)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=7
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUDIA_MPI,'(A5,A1,I2)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=8
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUDIA_MPI,'(A5,A1,I3)') LUDIA_MPIN,'.',MYPROC
         LUDIA_MPIL=9
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCI*_NODE MODULE")
      ENDIF
*#elif defined (VAR_MPI2)
      
#else
      NLUDIA_MPI=LUDIA_MPIN
      LUDIA_MPIL=5
#endif
*
      Open(Unit=LUDIA,File=NLUDIA_MPI(1:LUDIA_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. CI vector
      LUC = 61
      LUC_MPIN ="LUCVECT"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUC_MPI,'(A7,A1,I1)') LUC_MPIN,'.',MYPROC
         LUC_MPIL=9
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUC_MPI,'(A7,A1,I2)') LUC_MPIN,'.',MYPROC
         LUC_MPIL=10
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUC_MPI,'(A7,A1,I3)') LUC_MPIN,'.',MYPROC
         LUC_MPIL=11
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCI*_NODE MODULE")
      ENDIF
      ! master...
      IF(MYPROC == 0)THEN
        NLUC_MPI=LUC_MPIN
        LUC_MPIL=7
      END IF
#else
      NLUC_MPI=LUC_MPIN
      LUC_MPIL=7
#endif
      Open(Unit=LUC,File=NLUC_MPI(1:LUC_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. Sigma vector file
      LUHC = 62
      LUHC_MPIN ="HCFILE"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUHC_MPI,'(A6,A1,I1)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=8
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUHC_MPI,'(A6,A1,I2)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=9
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUHC_MPI,'(A6,A1,I3)') LUHC_MPIN,'.',MYPROC
         LUHC_MPIL=10
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCI*_NODE MODULE")
      ENDIF
#else
      NLUHC_MPI=LUHC_MPIN
      LUHC_MPIL=6
#endif
      Open(Unit=LUHC,File=NLUHC_MPI(1:LUHC_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
* =================
* Scratch files
* =================
      LUSC1 = 70
      LUSC1_MPIN ="LUSC1"
      LUSC2 = 71
      LUSC2_MPIN ="LUSC2"
      LUSC3 = 72
      LUSC3_MPIN ="LUSC3"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC1_MPI,'(A5,A1,I1)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A5,A1,I1)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A5,A1,I1)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=7
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC1_MPI,'(A5,A1,I2)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A5,A1,I2)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A5,A1,I2)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=8
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC1_MPI,'(A5,A1,I3)') LUSC1_MPIN,'.',MYPROC
         WRITE (NLUSC2_MPI,'(A5,A1,I3)') LUSC2_MPIN,'.',MYPROC
         WRITE (NLUSC3_MPI,'(A5,A1,I3)') LUSC3_MPIN,'.',MYPROC
         LUSC123_MPIL=9
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCI*_NODE MODULE")
      ENDIF
#else
      NLUSC1_MPI=LUSC1_MPIN
      NLUSC2_MPI=LUSC2_MPIN
      NLUSC3_MPI=LUSC3_MPIN
      LUSC123_MPIL=5
#endif
      Open(Unit=LUSC1,File=NLUSC1_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC2,File=NLUSC2_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
      Open(Unit=LUSC3,File=NLUSC3_MPI(1:LUSC123_MPIL),Status='UNKNOWN',
     &     Form='UNFORMATTED')
*. Scratch space for subspace handling
      LUSC34 = 73
      LUSC34_MPIN ="LUSC34"
      LUSC35 = 74
      LUSC35_MPIN ="LUSC35"
      LUSC36 = 75
      LUSC36_MPIN ="LUSC36"
      LUSC37 = 76
      LUSC37_MPIN ="LUSC37"
      LUSC38 = 77
      LUSC38_MPIN ="LUSC38"
      LUSC39 = 78
      LUSC39_MPIN ="LUSC39"
      LUSC40 = 79
      LUSC40_MPIN ="LUSC40"
      LUSC41 = 97
      LUSC41_MPIN ="LUSC41"
#if defined (VAR_MPI2)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUSC34_MPI,'(A6,A1,I1)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A6,A1,I1)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A6,A1,I1)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A6,A1,I1)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A6,A1,I1)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A6,A1,I1)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC40_MPI,'(A6,A1,I1)') LUSC40_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A6,A1,I1)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=8
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUSC34_MPI,'(A6,A1,I2)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A6,A1,I2)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A6,A1,I2)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A6,A1,I2)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A6,A1,I2)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A6,A1,I2)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC40_MPI,'(A6,A1,I2)') LUSC40_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A6,A1,I2)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=9
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUSC34_MPI,'(A6,A1,I3)') LUSC34_MPIN,'.',MYPROC
         WRITE (NLUSC35_MPI,'(A6,A1,I3)') LUSC35_MPIN,'.',MYPROC
         WRITE (NLUSC36_MPI,'(A6,A1,I3)') LUSC36_MPIN,'.',MYPROC
         WRITE (NLUSC37_MPI,'(A6,A1,I3)') LUSC37_MPIN,'.',MYPROC
         WRITE (NLUSC38_MPI,'(A6,A1,I3)') LUSC38_MPIN,'.',MYPROC
         WRITE (NLUSC39_MPI,'(A6,A1,I3)') LUSC39_MPIN,'.',MYPROC
         WRITE (NLUSC40_MPI,'(A6,A1,I3)') LUSC40_MPIN,'.',MYPROC
         WRITE (NLUSC41_MPI,'(A6,A1,I3)') LUSC41_MPIN,'.',MYPROC
         LUSC344_MPIL=10
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCI*_NODE MODULE")
      ENDIF
#else
      NLUSC34_MPI=LUSC34_MPIN
      NLUSC35_MPI=LUSC35_MPIN
      NLUSC36_MPI=LUSC36_MPIN
      NLUSC37_MPI=LUSC37_MPIN
      NLUSC38_MPI=LUSC38_MPIN
      NLUSC39_MPI=LUSC39_MPIN
      NLUSC40_MPI=LUSC40_MPIN
      NLUSC41_MPI=LUSC41_MPIN
      LUSC344_MPIL=6
#endif

      Open(Unit=LUSC34,File=NLUSC34_MPI(1:LUSC344_MPIL),
     &Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC35,File=NLUSC35_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC36,File=NLUSC36_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC37,File=NLUSC37_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC38,File=NLUSC38_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC39,File=NLUSC39_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC40,File=NLUSC40_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
      Open(Unit=LUSC41,File=NLUSC41_MPI(1:LUSC344_MPIL),
     &     Status='UNKNOWN',Form='UNFORMATTED')
* =================
* Output files
* =================
*. output file for CI-vectors
*. Not in use
      LUCIVO = 86
*. Natural orbitals in terms of input orbitals
*.
      LUMOUT = 87
*     Open(Unit=LUMOUT,File='CIAORB',Status='UNKNOWN',Form='FORMATTED')
*. Dumping 1- and 2- electron integrals in formatted form
* ( LU90 just defined here, it is not exported )
      LU90  = 88
*. Dumping symmmetry info, MO-AO expansion matrix and property integral
      LU91 = 69
*. CC amplitudes in formatted form
      LU_CCAMP = 90
*  LUCITA converted input unit for parallel calculations
      LU_LUWRT = 95
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      subroutine setunits
*
#include "clunit.inc"
#include "units.inc"
*
* Scratch file
      LUSC_INT = LU91
* Input for LUCITA, converted
      LUINP = LU_LUWRT
* COMFILE from MOLCAS-SCF
      LUCOM = LUSC36
* JOBIPH file from MOLCAS-RASSCF
      LUIPH = LUSC37
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TODSCN(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSC(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,0)
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TODSCNP(VEC,NREC,LREC,LBLK,LU)
*
* Write VEC as multiple record file accordin to NREC and LREC
* Only nonzero elements are written to DISC
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC(*)
      INTEGER LREC(NREC)
*
      IOFF = 1
      DO IREC = 1, NREC
C?      WRITE(6,*) ' TODSCN: IREC, LREC ',IREC,LREC(IREC)
C?      WRITE(6,*) ' Input record '
C?      CALL WRTMAT(VEC(IOFF),1,LREC(IREC),1,LREC(IREC))
        IF(LREC(IREC).GE.0) THEN
*. Normal complete record
          CALL ITODS(LREC(IREC),1,LBLK,LU)
          CALL TODSCP(VEC(IOFF),LREC(IREC),LBLK,LU)
          IOFF = IOFF + LREC(IREC)
        ELSE
*. zero record of length -LREC(IREC)
          CALL ITODS(-LREC(IREC),1,LBLK,LU)
          CALL ZERORC(IDUMMY,LU,1)
        END IF

*
      END DO
*
      RETURN
      END
