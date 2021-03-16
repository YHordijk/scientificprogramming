      module qcorr_cfg

!     variable block for Q-correction module for LUCI* CI modules. The naming convention is as follows:
!     _qcorr: input parameters
!     _q    : internal variables transferred from each CI module in the interface routine
!     written by S. Heinen and S. Knecht, ETH Zuerich winter 2013/14.

      implicit none

! logical block
      logical, public                          :: luci_cfg_qcorr   = .false.
      logical, public                          :: set_nash_q_qcorr = .true.
! integer block
      integer, public                          :: cvorb_qcorr(1:2) =  0 ! to allow for core orbitals in non-rel CI we should use nbsym instead for the dimensioning of cvorb
      integer, public                          :: wunit_q          = -1
      integer, public                          :: nref_qcorr       = -1
      integer, public                          :: nref_e_qcorr     = -1
      integer, public                          :: print_qcorr      = -1
      integer, public                          :: MXNSTR_q         =  0
      integer, public                          :: NMS2VAL_q        =  0
      integer, public                          :: NORB_q           =  0
      integer, public                          :: NSMST_q          =  0
      integer, public                          :: nfsym_q          =  0
      integer*8, public, allocatable           :: KNSTSO_q(:)
      integer*8, public, allocatable           :: KNSTSO2_q(:)
      integer, public,   allocatable           :: nash_q(:)
      integer, public,   allocatable           :: nash_qcorr(:)
      integer, public,   allocatable           :: ncorbs_mc(:)
      integer, public,   allocatable           :: MS2VAL_q(:)
      integer, public,   allocatable           :: NBLK_MS2_q(:)
      integer, public,   allocatable           :: NELEC_q(:)
      integer, public,   allocatable           :: IST_FOR_DT_q(:,:)

! real*8 block
      real(8), public, allocatable             :: ref_wavefunction_coeff_qcorr(:)
      real(8), public, allocatable             :: reference_energy_qcorr(:)
      real(8), public                          :: NACTEL_q                  = 0.0d0

! character block
      character (len=240), public, allocatable :: ref_wavefunction_qcorr(:,:)

      end module 
