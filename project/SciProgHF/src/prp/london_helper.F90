#ifdef MOD_LAO_REARRANGED
module london_helper

  implicit none

  save

  logical, public :: lao_lr_rearrange       = .false.
  logical, public :: shielding_rearrange    = .false.

contains

  subroutine set_london_keywords(word)
     character :: word*6

!    shiel2 keyword:
!    ---------------
!    calculate NMR shielding tensor in LAO basis from reorganized equations
!    Eq. 40-42 in JCP 131, 124119 (2009)
     if (word(1:6) == 'shiel2') then
        shielding_rearrange = .true.
     end if

!    laomod keyword:
!    ---------------
!    reorganize LR equations when LAOs are used, needed for connection-independent formulation of shielding
!    (also, automatically, shielding is calculated from Eq. 40-42 in JCP 131, 124119 (2009))
     if (word(1:6) == 'laomod') then
        lao_lr_rearrange = .true.
     end if

  end subroutine



end module
#endif
