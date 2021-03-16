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

      function ifind(i,arr,len)
C     Return the index of the first occurance of i in arr, or 0 if not
C     found.
      implicit none
      integer ifind, k,i,arr,len
      dimension arr(len)
      do k=1,len
         if (arr(k).eq.i) then
            ifind = k
            return
         endif
      enddo
      ifind = 0
      end

      LOGICAL FUNCTION ILL_NUMBER(x)
! Returns true if the input number is "ill", NaN,+Inf or -Inf.
! Calls internal C-function "is_ill_number(x)" in gp/gpc.c
      real*8, intent(in) :: x
      external :: is_ill_number
      if (is_ill_number(x) .ne. 0) then
          ILL_NUMBER = .True.
      else
          ILL_NUMBER = .False.
      endif
      return
      end 
