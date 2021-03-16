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

module adc_mat

! This module contains the matrix descriptor of the ADC matrix in
! the **current** run. All diagonalizers read the matrix information from
! this descriptor module.
! Additionally for the parallel run it contains local dimension descriptors

  implicit none

  save    !  if values of variables are changed they are still there after the next call 

  integer, public                  ::  reladc_md_iobase   = 0   !file handle for matrix handling
  integer, public                  ::  reladc_md_ionizl   = 0   !level od ionization for this matrix
  integer, public                  ::  reladc_md_ioldnew  = 0   !indicator for separate diag file
  integer, public                  ::  reladc_md_intbuf   = 0   !buffer length in matrix
  integer, public                  ::  reladc_md_desrep   = 0   !current symmetry
  integer, public                  ::  reladc_md_rcw      = 0   !real/complex matrix
  integer, public                  ::  reladc_md_lnzitr   = 0   !requested Krylov dimension
  integer, public                  ::  reladc_md_matdim   = 0   !matrix dimension
  integer, public                  ::  reladc_md_irecl    = 0   !record length of configuration data
  integer, public                  ::  reladc_md_nmain    = 0   !size of ADC matrix main block
  integer, public                  ::  reladc_md_nbufs    = 0   !number of buffers of len "intbuf"
  integer, public                  ::  reladc_md_neigenv  = 0   !number of requested full eigenvectors
  integer, public                  ::  reladc_md_davroots = 0   !number of requested Davidson roots
  integer, public                  ::  reladc_md_maxdavsp = 0   !maximum dimension of Davidson space
  integer, public                  ::  reladc_md_maxdavit = 0   !maximum number of Dav. iterations
  real*8,  public                  ::  reladc_md_convthr  = 0.0 !convergence of Davidson eigenvectors
  logical, public                  ::  reladc_md_davooc   = .true. !save memory in dav for large calculations
  character*6, public              ::  reladc_md_fileadc  = ''  !file name of ADC matrix
  character*6, public              ::  reladc_md_filediag = ''  !file name of diagonal of matrix, ioldnew=1
  character*6, public              ::  reladc_md_filecnf  = ''  !file name of configuration file
  character*5, public              ::  reladc_md_nmspec   = ''  !file name of spectral output file
  logical                          ::  reladc_md_isfano   = .false. !is this a fano run or not

  real, public, dimension(32,2)    ::  reladc_md_sip_eeigv= 0.0 !upper bound to eigenvalue for analysis
  real, public, dimension(32,2)    ::  reladc_md_dip_eeigv= 0.0 !upper bound to eigenvalue for analysis
  real, public                     ::  reladc_md_eeigv_lower = 0.0
  real, public                     ::  reladc_md_eeigv_upper = 0.0

end module

