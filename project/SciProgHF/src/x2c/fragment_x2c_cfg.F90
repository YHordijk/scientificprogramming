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
module fragment_x2c_cfg

! stefan: this module contains all required internal basis set, offset and 
!         pointer informations for the fragment-/or atomic-X2C approach.

  use x2c_fio

  implicit none

  public fragment_x2c
  public fragment_x2c_init
  public fragment_x2c_free

! type
  type fragment_x2c

!   double precision block
    real(8), allocatable ::      &
    pctmat(:)

!   integer block
    integer              ::      &
    charge,                      &
    funit
    integer, allocatable ::      &
    naosh_all(:),                &
    naosh_L(:)

!   logical block
    logical              ::                &
    fragment_approach_enabled,             &
    fragment_approach_ismolecule = .false.,&
    initialized
    
  end type

! fragment_x2c type
  type(fragment_x2c), public, save :: fragment_x2c_info
! ----------------------------------------------------------------------------

! parameters

  integer, parameter, public :: max_nr_symm_ind_centers = 2000

contains 

!*******************************************************************************
  subroutine fragment_x2c_init(A,                     & 
                               nfsym,                 &
                               nz,                    &
                               charge_cent,           &
                               naosh_all_cent,        &
                               naosh_L_cent)
!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------
    type(fragment_x2c)           :: A
    integer, intent(in)          :: nfsym
    integer, intent(in)          :: nz
    integer, intent(in)          :: charge_cent
    integer, intent(in)          :: naosh_all_cent(nfsym)
    integer, intent(in)          :: naosh_L_cent(nfsym)
!   ----------------------------------------------------------------------------
    integer                      :: pct_frag_dimension
    integer                      :: i
    character(len=3)             :: extension
    character(len=10)            :: filename
    character(len=12)            :: flabel
!   ----------------------------------------------------------------------------

!   reset old type information
    call fragment_x2c_free(A)

    A%initialized              = .true.

    allocate(A%naosh_all(nfsym))
    allocate(A%naosh_L(nfsym))

    A%naosh_all(1:nfsym) = -1
    A%naosh_L(1:nfsym)   = -1
    A%charge             = charge_cent
    A%funit              = 89

    pct_frag_dimension   = 0
    do i = 1, nfsym
      A%naosh_all(i)     = naosh_all_cent(i)
      A%naosh_L(i)       = naosh_L_cent(i)
      pct_frag_dimension = pct_frag_dimension + naosh_all_cent(i) * naosh_L_cent(i) * nz
    end do

    allocate(A%pctmat(pct_frag_dimension))

    A%pctmat         =  0

!   open file with fragment U matrix
    extension        = '000' 
    write(extension,'(i3)') A%charge

    if(extension(1:1) == ' ') write(extension(1:1),'(i1)') 0
    if(extension(2:2) == ' ') write(extension(2:2),'(i1)') 0
    write(filename,'(a7,a3)') 'X2CMAT.',extension

    open(A%funit,file=filename,status='old',form='unformatted',    &
         access='sequential',action='readwrite',position='rewind')

!   read atomic pctmat from file
    write(flabel,'(a7,i4,i1)') 'pctmtAO',1,1
    call x2c_read(flabel,                 &
                  A%pctmat,               &
                  A%naosh_all(1) *        &
                  A%naosh_L(1)   *        &
                  nz,                     &
                  A%funit)

  end subroutine fragment_x2c_init
!*******************************************************************************

  subroutine fragment_x2c_free(A)

!   ----------------------------------------------------------------------------
    type(fragment_x2c) :: A
!   ----------------------------------------------------------------------------

    if(.not. A%initialized) return

    close(A%funit, status='keep')

    A%initialized = .false.
    A%charge      = -1
    A%funit       = -1
    deallocate(A%naosh_all)
    deallocate(A%naosh_L)
    deallocate(A%pctmat)

  end subroutine fragment_x2c_free
!*******************************************************************************

end module fragment_x2c_cfg
