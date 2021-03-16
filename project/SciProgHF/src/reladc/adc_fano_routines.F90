MODULE adc_fano_routines

CONTAINS

  SUBROUTINE make_bigicra(aklrep,bigicra,iw)

!
!  Purpose: Write a more detailed icra table also containing absolute
!           spinor numbers for the 2h1p states
!
!           A
!           AREP
!           K
!           KREP
!           L
!           LREP
!           AABS
!           KABS
!           LABS
!

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)                        :: aklrep
    INTEGER, INTENT(INOUT), DIMENSION(9,NVOOT(aklrep)) :: bigicra!allocate outside
    INTEGER, INTENT(IN)                        :: iw

!
!  Data dictionary: Local Variables
!
    INTEGER    :: off, aoff, koff, loff     !offsets
    INTEGER    :: krep, lrep, arep          !representations
    INTEGER    :: k, l, a                   !relative spinor numbers
    INTEGER    :: kabs, labs, aabs          !absolute spinor numbers
    INTEGER    :: klrep                     !gamma_k x gamma_l
    INTEGER    :: kmin

!
!  Code
!
    bigicra = 0

!                                                                          
!  loop over irrep structure and form energy array & lookup table          
! 
    off = 1
    do_klrep:DO klrep = 1, nrep
      arep = multb(aklrep,klrep+nrep,2)
      aoff = iv(arep) + io(nrep+1)
      !WRITE(iw,*) arep! correct
      !WRITE(iw,*) 'iv(arep) = ', iv(arep)!correct
      do_lrep:DO lrep = 1, nrep
        loff = io(lrep)
        krep = multb(lrep,klrep+nrep,2)
        !WRITE(IW,*) 'lrep = ', lrep
        !WRITE(IW,*) 'krep = ', krep
        IF (krep < lrep) THEN
          CYCLE do_lrep
        END IF
        koff = io(krep)
        !WRITE(iw,*) 'koff = ', koff
        DO l = 1, no(lrep)
          kmin = 1
          IF(krep == lrep) THEN
            kmin = l + 1
          END IF
          DO k = kmin, no(krep)
            DO a = 1, nv(arep)
              bigicra(1,off) = a
              bigicra(2,off) = arep
              bigicra(3,off) = k
              bigicra(4,off) = krep
              bigicra(5,off) = l
              bigicra(6,off) = lrep
              bigicra(7,off) = aoff + a  !aabs
              bigicra(8,off) = koff + k  !kabs
              bigicra(9,off) = loff + l  !labs

              off = off + 1
            END DO
          END DO
        END DO
      END DO do_lrep
    END DO do_klrep
    
    IF (off-1 /= nvoot(aklrep)) THEN
      WRITE(IW,*) 
      WRITE(IW,*) ' *** Expected value of counter:',NVOOT(AKLREP)
      WRITE(IW,*) ' *** Actual value of counter:  ',OFF-1
      CALL QUIT('Inconsistency in make_bigicra')
    END IF

    WRITE(IW,*)
    WRITE(IW,*) 'Construct the big ICRA lookup table'
!    WRITE(IW,123) bigicra
!    123 FORMAT (' ',9I8)
    WRITE(IW,*)

  END SUBROUTINE make_bigicra




!
!--------------------------------------------------------------------
!
  INTEGER FUNCTION det_nr_configs(start,length,bigicra,krep,  &
                                  logvar,iw)
!
!  Purpose: Count the number of initial or final state configurations
!
    use adc_cfg

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)                           :: start, length
    INTEGER, INTENT(IN)                           :: krep
    INTEGER, INTENT(IN), DIMENSION(9,nvoot(krep)) :: bigicra
    INTEGER, INTENT(IN)                           :: logvar! 1=is_in, 2=is_fin
    INTEGER, INTENT(IN)                           :: iw

!
!  Data dictionary: Local Variables
!
    LOGICAL  :: is_in, is_coupl, is_fin
    INTEGER  :: aabs, kabs, labs
    INTEGER  :: asabs, ksabs, lsabs
    INTEGER  :: row
    INTEGER  :: inline
    INTEGER  :: first, second
    INTEGER  :: counter

!
!  Code
!

    counter = 0
     
    SELECT CASE (logvar)
      CASE (1) ! We want to now how many initial state configs we have
        counter = start - 1
        DO row = 1, (length - start + 1)
          aabs   = bigicra(7,row)
          kabs   = bigicra(8,row)
          labs   = bigicra(9,row)
          each_channel1:DO inline = 1, reladc_fano_nrchannels
            first  = reladc_fano_channels(inline,1)
            second = reladc_fano_channels(inline,2)
            IF (first == kabs) THEN
              IF (second == labs) THEN
                is_in = .FALSE. !We hit a final state configuration
              ELSE
                is_in = .TRUE.
              END IF
            ELSE
              is_in = .TRUE.
            END IF

            IF (is_in .EQV. .FALSE.) THEN
              EXIT each_channel1
            END IF
          END DO each_channel1
          IF (is_in .EQV. .TRUE.) THEN
            counter = counter + 1
          END IF
        END DO
      CASE (2) 
        DO row = 1, (length - start + 1)
          aabs   = bigicra(7,row)
          kabs   = bigicra(8,row)
          labs   = bigicra(9,row)
          each_channel2:DO inline = 1, reladc_fano_nrchannels
            first  = reladc_fano_channels(inline,1)
            second = reladc_fano_channels(inline,2)
            IF (first == kabs) THEN
              IF (second == labs) THEN
                is_fin = .TRUE. !We hit a final state configuration
              ELSE
                is_fin = .FALSE.
              END IF
            ELSE
              is_fin = .FALSE.
            END IF

            IF (is_fin .EQV. .TRUE.) THEN
              EXIT each_channel2
            END IF
          END DO each_channel2
          IF (is_fin .EQV. .TRUE.) THEN
            counter = counter + 1
            !WRITE(iw,*) counter !debug
            !WRITE(iw,*) 'aabs = ', aabs
            !WRITE(iw,*) 'kabs = ', kabs
            !WRITE(iw,*) 'labs = ', labs
          END IF
        END DO
      CASE DEFAULT ! Programmer has introduced a bug, since users can't
                   ! change this variable
        CALL QUIT ('state can only be initial or final state')
    END SELECT

    det_nr_configs = counter



  END FUNCTION det_nr_configs
    

!
!--------------------------------------------------------------------
!
  SUBROUTINE determine_block(row,col,bigicra,krep,is_in,  &
                             is_coupl,is_fin,block)
!
!  Purpose: to determine whether an matrix element belongs to an initial
!           or final state configuration. Hereby sort the 1h/2h1p space.
!           "block" is telling the program which part of the matrix is
!           currently being built.
!
    use adc_cfg

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)                           :: row, col
    INTEGER, INTENT(IN)                           :: krep
    INTEGER, INTENT(IN), DIMENSION(9,nvoot(krep)) :: bigicra
    LOGICAL, INTENT(INOUT)                          :: is_in, is_coupl, is_fin
    INTEGER, INTENT(IN)                           :: block ! 1=1h/2h1p, 2=2h1p/2h1p

!
!  Data dictionary: Local Variables
!
    INTEGER  :: aabs, kabs, labs
    INTEGER  :: asabs, ksabs, lsabs
    INTEGER  :: inline
    INTEGER  :: first, second

!
!  Code
!
    aabs   = bigicra(7,row)
    kabs   = bigicra(8,row)
    labs   = bigicra(9,row)

    asabs  = bigicra(7,col)
    ksabs  = bigicra(8,col)
    lsabs  = bigicra(9,col)


!
!  We now know the configuration of our current matrix element :-)
!  Now we need to decide what to do with this knowledge
!
    SELECT CASE (block)
      CASE (1) ! We are dealing with the h/2h1p block
        each_channel:DO inline = 1, reladc_fano_nrchannels
          first  = reladc_fano_channels(inline,1)
          second = reladc_fano_channels(inline,2)
          IF (first == kabs) THEN
            IF (second == labs) THEN
              is_in = .FALSE. !We hit a final state configuration
            ELSE
              is_in = .TRUE.
            END IF
          ELSE
            is_in = .TRUE.
          END IF

          IF (is_in .EQV. .FALSE.) THEN
            EXIT each_channel
          END IF
        END DO each_channel

!      CASE (2) ! We are dealing with the 2h1p/2h1p block
      CASE DEFAULT ! Programmer has introduced a bug, since users can't
                   ! change this variable
        CALL QUIT ('block in determine_block can only be 1 or 2')
    END SELECT

  END SUBROUTINE determine_block

!
!----------------------------------------------------------------------
!
  SUBROUTINE make_finalpos(krep,bigicra,finalpos,fin_max,iw)
!
!  Purpose: Create an array, that holds the positions of final state
!           configurations within the 2h1p part
!
    use adc_cfg
    use adc_fano_exchange
    use memory_allocator

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)                           :: krep
    INTEGER, INTENT(IN), DIMENSION(9,nvoot(krep)) :: bigicra
    INTEGER, INTENT(IN)                           :: fin_max
    INTEGER, INTENT(INOUT), DIMENSION(fin_max)    :: finalpos
    INTEGER, INTENT(IN)                           :: iw

!
!  Data dictionary: Local Variables
!
    INTEGER   :: row, inline
    INTEGER   :: aabs, kabs, labs
    INTEGER   :: first, second
    INTEGER   :: counter
    INTEGER   :: count2h1p
    INTEGER   :: channel
    INTEGER   :: gfirst, glast, off
    LOGICAL   :: is_fin

!
!  Code
!
    WRITE(iw,*) 'Enter make_finalpos'

    counter = 0
    off = 0

    CALL alloc(nr2h1p,reladc_fano_nrgroups)

    each_group:DO channel = 1, reladc_fano_nrgroups
      gfirst = off + 1
      glast  = off + reladc_fano_groups(channel)
      count2h1p = 0

      each_channel:DO inline = gfirst, glast
        first  = reladc_fano_channels(inline,1)
        second = reladc_fano_channels(inline,2)
        DO row = 1, nvoot(krep)
          aabs   = bigicra(7,row)
          kabs   = bigicra(8,row)
          labs   = bigicra(9,row)
          IF (first == kabs) THEN
            IF (second == labs) THEN
              is_fin            = .TRUE. !We hit a final state configuration
              counter           = counter + 1
              finalpos(counter) = row
              !WRITE(iw,*) 'position = ', row
              count2h1p = count2h1p + 1
            ELSE
              is_fin = .FALSE.
            END IF
          ELSE
            is_fin = .FALSE.
          END IF

          !IF (is_fin .EQV. .TRUE.) THEN
          !  CYCLE each_channel
          !END IF
        END DO 
      END DO each_channel

      nr2h1p(channel) = count2h1p
      WRITE(iw,*) 'Number of 2h1p states in this group', count2h1p
      off    = glast

    END DO each_group

  END SUBROUTINE make_finalpos



END MODULE adc_fano_routines
