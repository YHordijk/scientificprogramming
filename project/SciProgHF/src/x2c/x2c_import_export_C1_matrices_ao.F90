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
!
!
! module for import/export of matrices from/to C1* symmetry.
!
! written by sknecht feb 2011
!
module x2c_import_exportC1
  use x2cmod_cfg
  use x2c_fio
  use x2c_utils, only: print_x2cmat

  implicit none

  public x2c_export_mat2C1
  public x2c_export_pctmat2C1

  private

contains

!**********************************************************************
  subroutine x2c_export_mat2C1(h1_saao,                                &
                               naosh_ls,                               &
                               nz,                                     &
                               ipq_off_in,                             &
                               iqp_off_in,                             &
                               boson_to_fermion_sym_operator,          &
                               boson_irrep_in_matrix,                  &
                               quaternion_multiplication_table,        &
                               quaternion_pointer,                     &
                               quaternion_pointer_operator,            &
                               do_read_write,                          &
                               fh,                                     &
                               print_lvl)                               
!----------------------------------------------------------------------
!
!    purpose: export (defining h1) matrix in symmetry (quaternion packed)
!             AO basis to AO basis in C1* double group symmetry.
!
!----------------------------------------------------------------------
     real(8), intent(inout)          :: h1_saao(*)
     integer, intent(in)             :: naosh_ls
     integer, intent(in)             :: nz
     integer, intent(in)             :: ipq_off_in(4,0:7)
     integer, intent(in)             :: iqp_off_in(4,0:7)
     integer, intent(in)             :: boson_to_fermion_sym_operator(0:7,1:2)
     integer, intent(in)             :: boson_irrep_in_matrix(4,0:7)
     integer, intent(in)             :: quaternion_multiplication_table(4,4,4)
     integer, intent(in)             :: quaternion_pointer(0:7,2)
     integer, intent(in)             :: quaternion_pointer_operator(0:7)
     integer, intent(in)             :: print_lvl
     logical, intent(in)             :: do_read_write
     integer, intent(in)             :: fh
!--------------------------------------------------------------
     integer                         :: i
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     character (len=12)              :: flabel
     real(8), allocatable            :: scr1_mat(:)
!**********************************************************************

!      initialize dimensions and pointers which might be used outside this module in dirac calls
       call setdc2(0)

!      open file with 4c-h1
       if(do_read_write)then
         open(fh,file='x2ch1dirty',status='unknown',form='unformatted', &
              access='sequential',action='readwrite',position='rewind')
         call dzero(h1_saao,nz*naosh_ls**2)
         i = 0
         write(flabel,'(a7,i4,i1)') 'h1_4cao',1,0
         call x2c_read(flabel,h1_saao,naosh_ls*naosh_ls*nz,fh)
       end if

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(h1_saao,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - h1_4c-SO-DIRAC',6)
       end if

!      a. insert quaternion phase factors
!      ----------------------------------
       if(nz < 4)then
         do iz = 1, nz
           iq = ipq_off_in(nz,0)
           call q2bphase('F',iq,1,h1_saao(1+(naosh_ls*naosh_ls*(iz-1))))
         end do
       end if

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(h1_saao,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - h1_4c-SO-DIRAC-phase',6)
       end if

!      b. transform hamiltonian elements from DIRAC-sorted AO to Hermit-sorted AO basis
!      --------------------------------------------------------------------------------
       call bstobu_no_work(h1_saao,nz)

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(h1_saao,naosh_ls,naosh_ls,nz,ipq_off_in(1,0),'x2c - h1_4c-SO-Hermit-sorted',6)
       end if

!      allocate scratch matrix
       allocate(scr1_mat(naosh_ls*naosh_ls*4))
       call dzero(scr1_mat,naosh_ls*naosh_ls*4)

!      c. transform from SO-AO basis to AO-basis (nosym)
!      -------------------------------------------------
       i = boson_to_fermion_sym_operator(0,1)
       do iz = 1, 4
         irepd = boson_irrep_in_matrix(iz,0)
         iq    = quaternion_multiplication_table(1,quaternion_pointer(irepd,i),iz)
         ipq   = iqp_off_in(iq,0)
         call mtsoao(h1_saao( 1+(naosh_ls*naosh_ls*(ipq-1))),    &
                     scr1_mat(1+(naosh_ls*naosh_ls*(iz-1))),     &
                     naosh_ls,                                   &
                     irepd,                                      &
                     print_lvl)
       end do

!      debug print
       if(print_lvl > 2)then
         call print_x2cmat(scr1_mat,naosh_ls,naosh_ls,4,ipq_off_in(1,0),'x2c - h1_4c-AO-nosym',6)
       end if
      
       if(do_read_write)then
!        d. write h4c (nosym) in AO basis to file 
!        ----------------------------------------
         call x2c_write('h1_4caoNOSYM',scr1_mat,naosh_ls*naosh_ls*4,x2c_fh1dirty)
         close(x2c_fh1dirty, status='keep')
       end if

       deallocate(scr1_mat)

  end subroutine x2c_export_mat2C1

!**********************************************************************
  subroutine x2c_export_pctmat2C1()
!----------------------------------------------------------------------
!
!    purpose: export picture-change transformation matrix U matrix in symmetry (quaternion packed)
!             AO basis to AO basis in C1* double group symmetry.
!
!
!    note: the structure of U in symmetry-ordered Dirac AO-basis (SO-Dirac-AO) is as follows:
!
!    do i = 1, nr_fsym
!
!       U( [L+S](i), L(i)) == U( nfbas(i,0), nfbas(i,1))
!      
!    end do
!
!    in other words:
!
!          nfsym = 1   nfsym = 2
!
!            L_1         L_2
!         __                  _ 
!        |   L_1               |  
!        |   +           0     |  
!    U = |   S_1               |  
!        |               L_2   |  
!        |   0           +     |  
!        |               S_2   |  
!        |_                   _| 
!
!    we thus most likely need: 
!    - to define new AO-basis offset pointers to take properly care of this 
!      non-quadratic matrix structure.
!    - "reorder" the matrix elements for NZ > 1 thus that the "array"-order for the matrix elements is: 
!       first ALL elements for NZ == 1 and then all elements for NZ == 2, ...
!
!    ... to be continued with hjj next week.
!                         
!-----------------------------------------------------------------------
     real(8), allocatable            :: pctmat(:)
     real(8), allocatable            :: scrmat(:)
     integer                         :: i
     integer                         :: iq
     integer                         :: ipq
     integer                         :: iz
     integer                         :: irepd
     integer                         :: ioff_pct_mat
     real(8)                         :: mydummy(2)
     logical                         :: fndlab12
     character (len=12)              :: flabel
     character (len=23)              :: debug_label
!**********************************************************************

!      open file with pctmatC1
       open(x2c_fh1dirty,file='x2cpctmatC1',status='unknown',form='unformatted', &
            access='sequential',action='readwrite',position='rewind')
       if(fndlab12('pctmtaoNOSYM',x2c_fh1dirty))then
          print *,' pctmatC1 present on file'
       else
          call x2c_write('dummy       ',mydummy,-1,x2c_fh1dirty)
       end if

!      initialize dimensions and pointers which might be used outside this module in dirac calls
       call setdc1(2)
       call setdc2(2)

!      allocate pct-matrix
       allocate(pctmat(nr_ao_total_x2c*nr_ao_large_x2c*nr_quat))
       pctmat = 0

!      initialize offset
       ioff_pct_mat = 1

       do i = 1, nr_fsym
         write(flabel,'(a7,i4,i1)') 'pctmtAO',1,i
         call x2c_read(flabel,pctmat(ioff_pct_mat),nr_ao_all(i)*nr_ao_l(i)*nr_quat,x2c_funit)
         ioff_pct_mat  = ioff_pct_mat + nr_ao_all(i) * nr_ao_l(i)* nr_quat
       end do

!      DEBUG purification!!!
       do i = 1, nr_ao_total_x2c*nr_ao_large_x2c*nr_quat       
         if(abs(pctmat(i)).lt.1.0d-14) pctmat(i) = 0.0d0
       end do

!      debug print
       if(x2c_prt > 2)then 
         ioff_pct_mat = 1
         do i = 1, nr_fsym
           write(debug_label,'(a21,i2)') 'x2c - pctmat-orig -->',i
           call print_x2cmat(pctmat(ioff_pct_mat),nr_ao_all(i),nr_ao_l(i),nr_quat,x2c_cb_pq_to_uq, &
                             debug_label,6)
           ioff_pct_mat  = ioff_pct_mat + nr_ao_all(i) * nr_ao_l(i)* nr_quat
         end do
       end if

!      0. "expand" pctmat to full dimension (L+S,L) combining nr_fsym 1 and 2 in
!      the row dimension - this will allow a simpler loop structure in the
!      following routines ("gug" means gerade-ungerade)
!      if(nr_fsym > 1)then
!        do iz = 1, nr_quat
!          call expand_nonquadmat_LS_L_gug()
!        end do
!      end if

!      a. insert quaternion phase factors
!      ----------------------------------
       print *,'total matrix size is ==> ',nr_ao_total_x2c*nr_ao_large_x2c
       if(nr_quat < 4)then
         do iz = 1, nr_quat
           iq = x2c_cb_pq_to_uq(iz,0)
           call q2bphase_nonquadmat_LS_L('F',iq,1,pctmat(1+(nr_ao_total_x2c*nr_ao_large_x2c*(iz-1))))
         end do
       end if

!      debug print
       if(x2c_prt > 2)then 
         ioff_pct_mat = 1
         do i = 1, nr_fsym
           write(debug_label,'(a21,i2)') 'x2c - pctmx-phase -->',i
           call print_x2cmat(pctmat(ioff_pct_mat),nr_ao_all(i),nr_ao_l(i),nr_quat,x2c_cb_pq_to_uq, &
                             debug_label,6)
           ioff_pct_mat  = ioff_pct_mat + nr_ao_all(i) * nr_ao_l(i)* nr_quat
         end do
       end if

!      b. transform pctmat elements from DIRAC-sorted AO to Hermit-sorted AO basis
!      ---------------------------------------------------------------------------
       call bstobu_no_work_nonquadmat_LS_L(pctmat,nr_quat)

       if(x2c_prt > 2)then 
         call print_x2cmat(pctmat,nr_ao_total_x2c,nr_ao_large_x2c,nr_quat,x2c_cb_pq_to_uq, & 
                           'x2c - pctmat-after-bstobu',6)
       end if

!      c. transform from SO-AO basis to AO-basis (nosym)
!      -------------------------------------------------
       allocate(scrmat(nr_ao_total_x2c*nr_ao_large_x2c*4))
       i = x2c_bs_to_fs(0,1)
       do iz = 1, 4
         irepd = x2c_bs_irrep_mat(iz,0)
         iq    = x2c_iqmult_trip_q(1,x2c_pointer_quat(irepd,i),iz)
         ipq   = x2c_cb_uq_to_pq(iq,0)
         print *,'ipq-1 and irepd are ',ipq-1, irepd

         call MTSOAO_nonquadmatL(pctmat(1+(nr_ao_total_x2c*nr_ao_large_x2c*(ipq-1))),      &
                                 scrmat(1+(nr_ao_total_x2c*nr_ao_large_x2c*(iz-1))),       &
                                 nr_ao_total_x2c,                                          &
                                 nr_ao_large_x2c,                                          &
                                 irepd,                                                    &
                                 10)
       end do
!      debug print
       if(x2c_prt > 2)then 
         call print_x2cmat(scrmat,nr_ao_total_x2c,nr_ao_large_x2c,4,x2c_cb_pq_to_uq, & 
                           'x2c - pctmat-in-C1',6)
       end if       

!      d. write pctmat (nosym) in AO basis to file 
!      -------------------------------------------
       call x2c_write('pctmtaoNOSYM',scrmat,nr_ao_total_x2c*nr_ao_large_x2c*4,x2c_fh1dirty)

       deallocate(scrmat)
       deallocate(pctmat)

       close(x2c_fh1dirty, status='keep')

  end subroutine x2c_export_pctmat2C1

end module
