!======================================================================!
subroutine heapsorti(iarr,ndim,linvrt,imask,nmask)
!----------------------------------------------------------------------!
! HEAP SORT (non-recursive, in-place)
!
!   If you do not know anything about the INTEGER
!   vector you want to sort you can pick the HEAP SORT that scales 
!   with O(N log N) and is in-place, which means that it has no further 
!   storage requirements.
!     
! Benjamin Helmich-Paris
!----------------------------------------------------------------------!

 implicit none

! input:
 integer,intent(in) :: ndim

! output:
 integer,intent(inout) :: iarr(ndim)

! optional:
 logical,intent(in),optional :: linvrt
 integer,intent(in),optional :: nmask
 integer,intent(inout),optional :: imask(*)

! local scalar
 logical :: linvrt_loc
 integer :: nmsk,kdx,ival,imsk

 if (present(imask)) then
  if (present(nmask)) then
   nmsk = nmask
  else
   nmsk = 1
  end if
 else
  nmsk = 0
 end if

 if (present(linvrt)) then
  linvrt_loc = linvrt
 else
  linvrt_loc = .false.
 end if

 do kdx = ndim / 2, 1, -1
  call downheapi (iarr,imask,linvrt_loc,nmsk,ndim,kdx)
 end do

 kdx = ndim 
 do while (kdx > 1)
  ival = iarr(1)
  iarr(1) = iarr(kdx)
  iarr(kdx) = ival

  do imsk = 1,nmsk
   ival = imask(imsk)
   imask(imsk) = imask(imsk+(kdx-1)*nmsk)
   imask(imsk+(kdx-1)*nmsk) = ival
  end do

  kdx = kdx - 1
  call downheapi (iarr,imask,linvrt_loc,nmsk,kdx,1)
 end do

end subroutine heapsorti
!======================================================================!
!======================================================================!
subroutine downheapi(iarr,imask,linvrt,nmsk,ndim,kdx)

 implicit none

! dimensions:
 integer,intent(in)    :: nmsk,ndim,kdx

! input:
 logical,intent(in)    :: linvrt

! output:
 integer,intent(inout) :: iarr(ndim)
 integer,intent(inout) :: imask(nmsk,ndim)

! local:
 integer :: idx,jdx,imsk,jval,ival(nmsk), iswp(nmsk)

 idx = kdx
 jval = iarr(idx)
 do imsk = 1,nmsk
  ival(imsk) = imask(imsk,idx)
 end do

 if (linvrt) then

  do while (idx <= ndim / 2) 
   jdx = ishft(idx,1)
   if (jdx < ndim) then
    if (iarr(jdx) > iarr(jdx+1)) jdx = jdx + 1
   end if
   if (jval <= iarr(jdx)) exit
   iarr(idx) = iarr(jdx)
   do imsk = 1,nmsk
    iswp(imsk) = imask(imsk,jdx)
   end do
   do imsk = 1,nmsk
    imask(imsk,idx) = iswp(imsk)
   end do
   idx = jdx
  end do
  iarr(idx) = jval
  do imsk = 1,nmsk
   imask(imsk,idx) = ival(imsk)
  end do

 else

  do while (idx <= ndim / 2) 
   jdx = ishft(idx,1)
   if (jdx < ndim) then
    if (iarr(jdx) < iarr(jdx+1)) jdx = jdx + 1
   end if
   if (jval >= iarr(jdx)) exit
   iarr(idx) = iarr(jdx)
   do imsk = 1,nmsk
    iswp(imsk) = imask(imsk,jdx)
   end do
   do imsk = 1,nmsk
    imask(imsk,idx) = iswp(imsk)
   end do
   idx = jdx
  end do
  iarr(idx) = jval
  do imsk = 1,nmsk
   imask(imsk,idx) = ival(imsk)
  end do

 end if

end subroutine downheapi
!======================================================================!
!======================================================================!
subroutine heapsortd (darr,ndim,linvrt,imask,nmask)
!----------------------------------------------------------------------!
! HEAP SORT (non-recursive, in-place)
!
!   If you do not know anything about the REAL or DOUBLE PRECISION 
!   vector you want to sort you can pick the HEAP SORT that scales 
!   with O(N log N) and is in-place, which means that it has no further 
!   storage requirements.
!     
! Benjamin Helmich-Paris
!----------------------------------------------------------------------!

 implicit none

! input:
 integer,intent(in) :: ndim

! output:
 real(8),intent(inout) :: darr(ndim)

! optional:
 logical,intent(in),optional :: linvrt
 integer,intent(in),optional :: nmask
 integer,intent(inout),optional :: imask(*)

! local:
 logical :: linvrt_loc
 integer :: nmsk,kdx,ival,imsk
 real(8)  :: dval

 if (present(imask)) then
  if (present(nmask)) then
   nmsk = nmask
  else
   nmsk = 1
  end if
 else
  nmsk = 0
 end if

 if (present(linvrt)) then
  linvrt_loc = linvrt
 else
  linvrt_loc = .false.
 end if

 do kdx = ndim / 2, 1, -1
  call downheapd (darr,imask,linvrt_loc,kdx,ndim,nmsk)
 end do

 kdx = ndim 
 do while (kdx > 1)
  dval = darr(1)
  darr(1) = darr(kdx)
  darr(kdx) = dval

  do imsk = 1,nmsk
   ival = imask(imsk)
   imask(imsk) = imask(imsk+(kdx-1)*nmsk)
   imask(imsk+(kdx-1)*nmsk) = ival
  end do

  kdx = kdx - 1
  call downheapd (darr,imask,linvrt_loc,1,kdx,nmsk)
 end do

end subroutine heapsortd
!======================================================================!
!======================================================================!
subroutine downheapd (darr,imask,linvrt,kdx,ndim,nmsk)

 implicit none

! dimensions:
 integer,intent(in)    :: ndim, nmsk

! input:
 integer,intent(in)    :: kdx
 logical,intent(in)    :: linvrt

! output:
 real(8),intent(inout)  :: darr(ndim)
 integer,intent(inout) :: imask(nmsk,ndim)

 integer :: idx,jdx,ival(nmsk), imsk, iswp(nmsk)
 real(8)  :: dval

 idx = kdx
 dval = darr(idx)
 do imsk = 1,nmsk
  ival(imsk) = imask(imsk,idx)
 end do

 if (linvrt) then

  do while (idx <= ndim / 2) 
   jdx = ishft(idx,1)
   if (jdx < ndim) then
    if  (darr(jdx) > darr(jdx+1)) jdx = jdx + 1
   end if
   if (dval <= darr(jdx)) exit
   darr(idx) = darr(jdx)
   do imsk = 1,nmsk
    iswp(imsk) = imask(imsk,jdx)
   end do
   do imsk = 1,nmsk
    imask(imsk,idx) = iswp(imsk)
   end do
   idx = jdx
  end do
  darr(idx) = dval
  do imsk = 1,nmsk
   imask(imsk,idx) = ival(imsk)
  end do

 else

  do while (idx <= ndim / 2) 
   jdx = ishft(idx,1)
   if (jdx < ndim) then
    if (darr(jdx) < darr(jdx+1)) jdx = jdx + 1
   end if
   if (dval >= darr(jdx)) exit
   darr(idx) = darr(jdx)
   do imsk = 1,nmsk
    iswp(imsk) = imask(imsk,jdx)
   end do
   do imsk = 1,nmsk
    imask(imsk,idx) = iswp(imsk)
   end do
   idx = jdx
  end do
  darr(idx) = dval
  do imsk = 1,nmsk
   imask(imsk,idx) = ival(imsk)
  end do

 end if

end subroutine downheapd
!======================================================================!
