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

module moltra_labeling

  integer, parameter, private ::  MAX_KRAMERS_PAIRS = 3200

  save

! Arrays needed to toggle between different indexing (TODO: figure out how to use parameters below)
  integer, public  :: moltra_isp (-MAX_KRAMERS_PAIRS:MAX_KRAMERS_PAIRS)
  integer, public  :: moltra_ikr (2*MAX_KRAMERS_PAIRS)
  integer, public  :: moltra_nkr
 
contains

      subroutine Make_Kramer_to_SpinorIndex (nfsym,nstr)
      integer ::    nfsym, nstr(nfsym), ii, jj, ifsym, i
      ii = 0
      jj = 0
      do ifsym = 1, nfsym
         do i = 1, nstr(ifsym)
            ii = ii + 1
            jj = jj + 1
            moltra_isp( ii) = jj
            moltra_isp(-ii) = jj + nstr(ifsym)
            moltra_ikr(jj)              =  ii
            moltra_ikr(jj+nstr(ifsym))  = -ii
         end do
         jj = jj + nstr(ifsym)
         if (jj > MAX_KRAMERS_PAIRS) call quit ('Update moltra_labeling')
      end do
      moltra_nkr = ii
      end subroutine

      integer function Kramer_to_SpinorIndex (ikr)
      Kramer_to_SpinorIndex =  moltra_isp(ikr)
      end function

      integer function Spinor_to_KramerIndex (isp)
      Spinor_to_KramerIndex =  moltra_ikr(isp)
      end function

end module
