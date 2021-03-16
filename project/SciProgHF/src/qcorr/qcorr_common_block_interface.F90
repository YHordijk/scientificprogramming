       module qcorr_interface

       use qcorr_cfg

       implicit none

       public init_qcorr
       public exit_qcorr


contains

      subroutine init_qcorr(relQC, print_header)

      logical, intent(in)           :: relQC
      logical, intent(in), optional :: print_header
#include "priunit.h"

      wunit_q = lupri

      if(print_header) call print_header_qcorr(wunit_q)

      select case(relQC)

        case(.true.)
#ifdef PRG_DIRAC
          call init_qcorr_rel()
#else
          call quit('no relativistic CI in this program environment')
#endif
        case(.false.)
          call init_qcorr_nr()

      end select
      
      end subroutine init_qcorr

      subroutine exit_qcorr(relQC)

      logical, intent(in) :: relQC

      deallocate(KNSTSO_q)
      deallocate(KNSTSO2_q)
      deallocate(NELEC_q)
      deallocate(nash_q)
      if(relQC)then
        deallocate(MS2VAL_q)
        deallocate(NBLK_MS2_q)
        deallocate(IST_FOR_DT_q)
      end if
      
      end subroutine exit_qcorr

#ifdef PRG_DIRAC
      subroutine init_qcorr_rel()

 use symmetry_setup_krci
#include "dgroup.h"
#include "../luciarel/ipoist8.inc"
#include "../luciarel/mxpdim.inc"
#include "../luciarel/strbas.inc"
#include "../luciarel/cstate.inc"
#include "../luciarel/strinp.inc"
#include "../luciarel/stinf.inc"
#include "../luciarel/cgas.inc"
#include "../luciarel/gasstr.inc"

! integer block
      MXNSTR_q    = MXNSTR
      NMS2VAL_q   = NMS2VAL
      NACTEL_q    = real(NACTEL)
      NORB_q      = 0
      NSMST_q     = NSMST
      nfsym_q     = nfsym

      allocate(nash_q(nfsym_q))
      call helper_cb(nash_q,nfsym_q)

      allocate(KNSTSO_q(MXPSTT))
      KNSTSO_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
      
      allocate(KNSTSO2_q(MXPSTT))
      KNSTSO2_q(1:MXPSTT) = KNSTSO2(1:MXPSTT) 
     
      allocate(NELEC_q(MXPSTT))
      NELEC_q(1:MXPSTT) = NELEC(1:MXPSTT)

      allocate(MS2VAL_q(100))
      MS2VAL_q(1:100) = MS2VAL(1:100)

      allocate(NBLK_MS2_q(100))
      NBLK_MS2_q(1:100) = NBLK_MS2(1:100)

      allocate(IST_FOR_DT_q(2,MXPSTT))
      IST_FOR_DT_q(1:2,1:mxpstt) = IST_FOR_DT(1:2,1:MXPSTT)

      end subroutine init_qcorr_rel
#endif

      subroutine init_qcorr_nr()

#ifdef PRG_DIRAC
#include "dgroup.h"
#else
#include "infpar.h"
#include "inforb.h"
#endif
#include "../lucita/mxpdim.inc"

! contains NACTEL
#include "../lucita/lucinp.inc"

#include "../lucita/strbas.inc"
#include "../lucita/cstate.inc"
#include "../lucita/strinp.inc"
#include "../lucita/stinf.inc"
#include "../lucita/cgas.inc"
#include "../lucita/gasstr.inc"
#include "../lucita/orbinp.inc"
#include "../lucita/csm.inc"


! integer block
      MXNSTR_q    = MXNSTR
      NMS2VAL_q   = 1 
      NACTEL_q    = real(NACTEL)
      NORB_q      = ntoob
      NSMST_q     = NSMST
#ifdef PRG_DIRAC
      nfsym_q     = nfsym

#else
      nfsym_q     = nsym
#endif
      allocate(nash_q(nfsym_q))
      call helper_cb(nash_q,nfsym_q)

      allocate(KNSTSO_q(MXPSTT))
      KNSTSO_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
      
      allocate(KNSTSO2_q(MXPSTT))
      KNSTSO2_q(1:MXPSTT) = KNSTSO(1:MXPSTT) 
     
      allocate(NELEC_q(MXPSTT))
      NELEC_q(1:MXPSTT) = NELEC(1:MXPSTT)

      end subroutine init_qcorr_nr
 
      subroutine helper_cb(nash_q,nfsym_q)
#ifdef PRG_DIRAC
#include "dcborb.h"
#else
#include "inforb.h"
#endif
        integer, intent(out) :: nash_q(nfsym_q)
        integer, intent(in)  :: nfsym_q
        nash_q(1:nfsym_q) = nash(1:nfsym_q)
      end subroutine helper_cb

      subroutine print_header_qcorr(wunit)

      integer, intent(in) :: wunit

      write(wunit,'(/18x,a)') ' ********************************************************************'
      write(wunit,'(18x,a )') ' ***                        QCORR-CI                              ***'
      write(wunit,'(18x,a )') ' ***      a module for Davidson(-type) energy corrections of      ***'
      write(wunit,'(18x,a )') ' ***                  (MR)CISD wave functions                     ***'
      write(wunit,'(18x,a )') ' ***                                                              ***'
      write(wunit,'(18x,a )') ' ***  authors:         - Stefan Knecht                            ***'
      write(wunit,'(18x,a )') ' ***                   - Stefan Heinen                            ***'
      write(wunit,'(18x,a )') ' ***  contributor:     - Markus Reiher                            ***'
      write(wunit,'(18x,a )') ' ***                                                              ***'
      write(wunit,'(18x,a )') ' ***  please cite: work in progress, private communication with   ***'
      write(wunit,'(18x,a )') ' ***               S. Knecht, ETH Zuerich, 2014.                  ***'
      write(wunit,'(18x,a )') ' ***                                                              ***'
      write(wunit,'(18x,a )') ' ***  reference for +Q models: Chem. Rev. 112, 108-181 (2012)     ***'
      write(wunit,'(18x,a )') ' ***                                                              ***'
      write(wunit,'(18x,a/)') ' ********************************************************************'

      end subroutine print_header_qcorr

      end module

