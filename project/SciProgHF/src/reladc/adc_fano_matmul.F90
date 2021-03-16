MODULE adc_fano_matmul

CONTAINS

  SUBROUTINE fano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                               fhmat,fhdiag,intbuf,e2h1p,t_mom,nbufs, &
                               do_partial,assign_fin_channel,iw)
!
! Purpose: Do the final multiplication of the initial wave function
!          with the coupling block and the final state wavefunctions.
!
    use adc_fano_exchange
    use memory_allocator
    use adc_cfg
#ifdef HAS_STIELTJES
    use stieltjesmod
#endif

    IMPLICIT NONE
#include "pi.h"
!
! Calling Variables
!
    CHARACTER(6), INTENT(IN)                        :: fn_adcdiag,fn_adcmat
    INTEGER, INTENT(IN)                             :: ladc
    INTEGER, INTENT(IN)                             :: nr_in, nr_fin
    INTEGER, INTENT(IN)                             :: intbuf
    INTEGER, INTENT(IN)                             :: fhmat, fhdiag
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: e2h1p
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: t_mom
    INTEGER, INTENT(IN)                             :: nbufs
    LOGICAL, INTENT(IN)                             :: do_partial
    INTEGER, INTENT(IN), DIMENSION(nr_fin)          :: assign_fin_channel
    INTEGER, INTENT(IN)                             :: iw

!
! Local Variables
!
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: diag
    REAL(8)                               :: autoev  = 27.2113957d0
    INTEGER                               :: ixx
!    INTEGER                               :: nbufs
    INTEGER                               :: nact
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: buf
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: ioi, ioj
    INTEGER                               :: idummy
    REAL(8)                               :: buffer
    INTEGER                               :: irow, jcol
    INTEGER                               :: i, j, k
    INTEGER                               :: counter
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: vec_start
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: vec1, vec2

    CHARACTER(3), PARAMETER               :: nmspec = 'FPS'
    CHARACTER(11)                         :: fn_trans
    INTEGER                               :: fh_trans
    CHARACTER(3), PARAMETER               :: nmstie = 'STI'
    CHARACTER(21)                         :: fn_stieltjes
    INTEGER                               :: ierror
    INTEGER                               :: channel
    INTEGER                               :: first, last, off

    REAL(8)                               :: e_init
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: total_gamma
    REAL(8), ALLOCATABLE, DIMENSION(:,:)  :: partial_gamma
    REAL(8)                               :: gammae !result from stieltjes
    REAL(8)                               :: sum_partial
    REAL(8)                               :: norm


!
! Code
!

    CALL alloc(buf,intbuf)
    CALL alloc(ioi,intbuf)
    CALL alloc(ioj,intbuf)

    CALL alloc(vec_start,ladc) !to be deallocated
    CALL alloc(vec1,ladc)
    CALL alloc(vec2,ladc)

    CALL alloc(e2h1p,nr_in_evecs*nr_fin)
    CALL alloc(t_mom,nr_in_evecs*nr_fin)

    CALL alloc(total_gamma,nr_in_evecs)
    CALL alloc(partial_gamma,reladc_fano_nrgroups,nr_in_evecs)

    WRITE(iw,*) 'Number of Initial States: ', nr_in_evecs

    vec_start = 0

    OPEN(fhmat,FILE=fn_adcmat,FORM='unformatted',STATUS='unknown')
!    OPEN(fhdiag,FILE=fn_adcdiag,FORM='unformatted',STATUS='unknown')

    DO i = 1, nr_in_evecs
      counter = 0
      vec1 = 0
      vec_start(1:nr_in) = in_evecs(1:nr_in,i)
!      WRITE(iw,*) 'starting vector'
!      WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)

! Choose file name of the pseudo-spectrum
      IF(i.GT.9) THEN
        WRITE(fn_trans,'(a3,a1,i2)') nmspec,'.',i
      ELSE
        WRITE(fn_trans,'(a3,a2,i1)') nmspec,'.0',i
      END IF
! Choose file name for Stieltjes output
      IF(i.GT.9) THEN
        WRITE(fn_stieltjes,'(a3,a1,i2)') nmstie,'.',i
      ELSE
        WRITE(fn_stieltjes,'(a3,a2,i1)') nmstie,'.0',i
      END IF


!      CALL alloc(diag,ladc)
!      REWIND(fhdiag)
      REWIND(fhmat)
!      READ(fhdiag) (diag(ixx), ixx=1, ladc), nbufs
!      DO j = 1, ladc
!        vec1(j) = vec1(j) + diag(j) * vec_start(j)
!      END DO
!      CALL dealloc(diag)
!      !WRITE(iw,*) 'Vector with diagonal'
!      !WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
!      WRITE(iw,*) 'nbufs =', nbufs

      DO k = 1, nbufs
        READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                    (ioi(ixx),ixx=1,intbuf), &
                    (ioj(ixx),ixx=1,intbuf), &
                    nact, idummy
!        WRITE(iw,*) 'nact = ', nact
        DO j = 1, nact
          irow   = ioi(j)
          jcol   = ioj(j)
          buffer = buf(j)
!          WRITE(iw,*) 'buffer = ', buffer
          !WRITE(iw,*) 'irow = ', irow, ' jcol = ', jcol

          vec1(irow) = vec1(irow) + buffer * vec_start(jcol)
          vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
        END DO
      END DO
!      WRITE(iw,*) 'Vector with total sorted matrix'
!      WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
      
!
! Now we have the initial state vector multiplied with the sorted ADC matrix.
! We will further multiply this vector to each and every final state vector.
!
      DO j = 1, nr_fin
        vec2 = 0
        counter = counter + 1
        vec2(ladc-nr_fin+1:ladc) = fin_evecs(1:nr_fin,j)
        !WRITE(iw,*) 'Final state vector'
        !WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)
        e2h1p(counter) = fin_energies(j)
        ! This part is ony necessary, when the initial and final states
        ! are non-orthogonal, which they are not at the moment!
        !t_mom(counter) = DOT_PRODUCT(vec1,vec2) &
        !                 - DOT_PRODUCT(vec_start,vec2) * in_energy(i)
        t_mom(counter) = 2*pi*DOT_PRODUCT(vec1,vec2)*DOT_PRODUCT(vec1,vec2)
      END DO

      WRITE(iw,*)
      WRITE(iw,*) '-------------------------------------------------------------'
      WRITE(iw,'(1X,A45,I4)') 'Pseudo-spectrum for initial state number ', i
      WRITE(iw,'(1X,A10,F10.6,A3,A20,F10.6)') 'energy: ', in_energy(1,i)*autoev,&
                                              'eV', 'polestrength: ', in_energy(2,i)
      WRITE(iw,*) '-------------------------------------------------------------'
      WRITE(iw,*)
      WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'
      IF (counter > 20) THEN
        WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,20)
      ELSE
        WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
      END IF

      WRITE(iw,*)
      WRITE(iw,*) 'Shifting pseudo-spectrum'
      FORALL (ixx = 1:counter)
        e2h1p(ixx) = e2h1p(ixx) - in_energy(1,i)
        !t_mom(ixx) = 2 * pi * t_mom(ixx) * t_mom(ixx)
      END FORALL

      WRITE(iw,*)
      WRITE(iw,*) 'Calling Stieltjes.'
      !e_init = in_energy(1,i)
      e_init = 0.0
#ifdef HAS_STIELTJES
      CALL stieltjes(counter,e2h1p,t_mom,e_init,fn_stieltjes,&
                     gammae,2004,0)
#else
      call quit('stieltjes not supported in this version')
#endif
!
! printflag 1 for no external output(dirac.out), 2 for output
!

      WRITE(iw,*) 'Gamma = ', gammae*autoev, 'eV'
      total_gamma(i) = gammae * autoev
!!
!!
      fh_trans = fhdiag
      OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    &
           ACTION='WRITE',iostat=ierror)
      WRITE(fh_trans,*) in_energy(1,i)
      WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
      CLOSE(fh_trans)

    END DO !loop over intial states


!
! Starting the partial decay widths
!
!

    IF (do_partial) THEN
      WRITE(iw,*) 
      WRITE(iw,*) ' Starting creation of partial decay width pseudo-spectra'
      WRITE(iw,*)

      DO i = 1, nr_in_evecs
        last = nr_in
        vec_start(1:nr_in) = in_evecs(1:nr_in,i)
!        WRITE(iw,*) 'starting vector'
!        WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)

        DO channel = 1, reladc_fano_nrgroups
          counter = 0
          vec1  = 0
          first = last + 1
          last  = first + nr2h1p(channel) - 1
          !WRITE(iw,*) 'first = ', first, ' last = ', last

          REWIND(fhmat)

! Ch  oose file name of the pseudo-spectrum
          IF(i.GT.9) THEN
            WRITE(fn_trans,'(a3,a1,i2,a1,a4)') nmspec,'.',i,'.',reladc_fano_labels(channel)
          ELSE
            WRITE(fn_trans,'(a3,a2,i1,a1,a4)') nmspec,'.0',i,'.',reladc_fano_labels(channel)
          END IF
! Choose file name of Stieltjes output
          IF(i.GT.9) THEN
            WRITE(fn_trans,'(a3,a1,i2,a1,a4)') nmstie,'.',i,'.',reladc_fano_labels(channel)
          ELSE
            WRITE(fn_trans,'(a3,a2,i1,a1,a4)') nmstie,'.0',i,'.',reladc_fano_labels(channel)
          END IF


          DO k = 1, nbufs
            READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                        (ioi(ixx),ixx=1,intbuf), &
                        (ioj(ixx),ixx=1,intbuf), &
                        nact, idummy
!            WRITE(iw,*) 'nact = ', nact
            DO j = 1, nact
              irow   = ioi(j)
              jcol   = ioj(j)
              buffer = buf(j)
!              WRITE(iw,*) 'buffer = ', buffer

              IF (first <= irow .AND. irow <= last) THEN
                vec1(irow) = vec1(irow) + buffer * vec_start(jcol)
                vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
              END IF
            END DO
          END DO
!          IF (channel == 2) THEN
!            WRITE(iw,*) 'Vector with total sorted matrix'
!            WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
!          END IF
        
!
! No  w we have the initial state vector multiplied with the sorted ADC matrix.
! We   will further multiply this vector to each and every final state vector.
!
          mul_fin:DO j = 1, nr_fin
              vec2 = 0
              counter = counter + 1
              vec2(ladc-nr_fin+1:ladc) = fin_evecs(1:nr_fin,j)
              !IF (j == 1) THEN
              !  WRITE(iw,*) 'Final state vector'
              !  WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)
              !END IF
              e2h1p(counter) = fin_energies(j)
              ! This part is ony necessary, when the initial and final states
              ! are non-orthogonal, which they are not at the moment!
              !t_mom(counter) = DOT_PRODUCT(vec1,vec2) &
              !                 - DOT_PRODUCT(vec_start,vec2) * in_energy(i)
              t_mom(counter) = 2*pi*DOT_PRODUCT(vec1,vec2)*DOT_PRODUCT(vec1,vec2)
!              WRITE(iw,*) 't_mom: ', t_mom(counter)
          END DO mul_fin

! This if clause has to be out around the multiplication of final
! states in case of an selection of final states
!            IF (channel == assign_fin_channel(j)) THEN
          
          WRITE(iw,*)
          WRITE(iw,*) 'name: ', fn_trans
          WRITE(iw,*) '-------------------------------------------------------------'
          WRITE(iw,'(1X,A45,I4)') 'Pseudo-spectrum for initial state number ', i
          WRITE(iw,'(1X,A10,F10.6,A3,A20,F10.6)') 'energy: ', in_energy(1,i)*autoev,&
                                                  'eV', 'polestrength: ', in_energy(2,i)
          WRITE(iw,*) '-------------------------------------------------------------'
          WRITE(iw,*)
          WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'
          IF (counter > 20) THEN
            WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,20)
          ELSE
            WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
          END IF
          WRITE(iw,*)

!
! This procedure only makes sense, when those final vectors not
! associated with the chosen final states are not considered in
! the calculation of the transition moments.
!
          FORALL (ixx = 1:counter)
            e2h1p(ixx) = e2h1p(ixx) - in_energy(1,i)
            !t_mom(ixx) = 2 * pi * t_mom(ixx) * t_mom(ixx)
          END FORALL

          WRITE(iw,*) 'Channel #', channel
          WRITE(iw,*)
          WRITE(iw,*) 'Calling Stieltjes.'
          WRITE(iw,*)
          WRITE(iw,*)
          !e_init = in_energy(1,i)
          e_init = 0.0
#ifdef HAS_STIELTJES
          CALL stieltjes(counter,e2h1p,t_mom,e_init,fn_stieltjes,&
                         gammae,1000+channel,0)
#else
      call quit('stieltjes not supported in this version')
#endif
!
! printflag 1 for no external output(dirac.out), 2 for output
!
          WRITE(iw,*) 'Gamma = ', gammae*autoev, 'eV'
          partial_gamma(channel,i) = gammae*autoev
!!
!!
          OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    &
               ACTION='WRITE',iostat=ierror)
          WRITE(fh_trans,*) in_energy(1,i)
          WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
          CLOSE(fh_trans)

        END DO ! end of loop over channels
      END DO   ! end of loop over initial states

    END IF

    CALL dealloc(buf)
    CALL dealloc(ioi)
    CALL dealloc(ioj)

    CALL dealloc(vec_start)
    CALL dealloc(vec1)
    CALL dealloc(vec2)

    CALL dealloc(in_evecs)
    CALL dealloc(fin_evecs)
    CALL dealloc(fin_energies)

    CALL dealloc(e2h1p)
    CALL dealloc(t_mom)

!
! Beautiful Summary of the calculated decay widths
!

    WRITE(iw,*)
    WRITE(iw,*)

    DO i = 1, nr_in_evecs
      sum_partial = 0.0
      DO j = 1, reladc_fano_nrgroups
        sum_partial = sum_partial + partial_gamma(j,i)
      END DO
      norm = sum_partial / total_gamma(i)
      partial_gamma(1:reladc_fano_nrgroups,i) = partial_gamma(1:reladc_fano_nrgroups,i) &
                                                / norm
      WRITE(iw,149) 'Initial state # ', i, 'energy: ', in_energy(1,i)*autoev, 'eV', &
                  'polestr.: ', in_energy(2,i)
      WRITE(iw,*) '--------------------------------------------------------------'
      WRITE(iw,*) 'Total decay width: ', total_gamma(i), 'eV'
      WRITE(iw,*) 'Sum of partial decay widths: ', sum_partial, 'eV'
      WRITE(iw,*) 'Renormalized partial decay widths in eV:'
      WRITE(iw,*) '--------------------------------------------------------------'
      DO channel = 1, reladc_fano_nrgroups
        WRITE(iw,150) 'Channel #', channel, partial_gamma(channel,i)
      END DO
      WRITE(iw,*) '--------------------------------------------------------------'
      WRITE(iw,*)
      WRITE(iw,*)
    END DO

    149 FORMAT(' ',A16,2X,I3,3X,A8,2X,F14.5,1X,A2,3X,A10,2X,F7.5)
    150 FORMAT(' ',A9,2X,I3,3X,ES14.5)

    CALL dealloc(total_gamma)
    CALL dealloc(partial_gamma)
    CALL dealloc(in_energy)

  END SUBROUTINE fano_create_t_mom



!
!-----------------------------------------------------------------
!
  SUBROUTINE cfano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                               fhmat,fhdiag,intbuf,e2h1p,t_mom,nbufs, &
                               do_partial,assign_fin_channel,iw)
!
! Purpose: Do the final multiplication of the initial wave function
!          with the coupling block and the final state wavefunctions.
!

    use adc_fano_exchange
    use memory_allocator
    use adc_cfg
#ifdef HAS_STIELTJES
    use stieltjesmod
#endif

    IMPLICIT NONE
#include "pi.h"

!
! Calling Variables
!
    CHARACTER(6), INTENT(IN)                        :: fn_adcdiag,fn_adcmat
    INTEGER, INTENT(IN)                             :: ladc
    INTEGER, INTENT(IN)                             :: nr_in, nr_fin
    INTEGER, INTENT(IN)                             :: intbuf
    INTEGER, INTENT(IN)                             :: fhmat, fhdiag
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: e2h1p
    REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: t_mom
    INTEGER, INTENT(IN)                             :: nbufs
    LOGICAL, INTENT(IN)                             :: do_partial
    INTEGER, INTENT(IN), DIMENSION(nr_fin)          :: assign_fin_channel
    INTEGER, INTENT(IN)                             :: iw

!
! Local Variables
!
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: diag
    REAL(8)                               :: autoev  = 27.2113957d0
    INTEGER                               :: ixx
    INTEGER                               :: nact
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: buf
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: ioi, ioj
    INTEGER                               :: idummy
    COMPLEX(8)                               :: buffer
    INTEGER                               :: irow, jcol
    INTEGER                               :: i, j, k
    INTEGER                               :: counter
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: vec_start
    COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: vec1, vec2

    CHARACTER(3),PARAMETER                :: nmspec = 'FPS'
    CHARACTER(11)                         :: fn_trans
    INTEGER                               :: fh_trans
    CHARACTER(3), PARAMETER               :: nmstie = 'STI'               
    CHARACTER(21)                         :: fn_stieltjes
    INTEGER                               :: ierror

    INTEGER                               :: channel                      
    INTEGER                               :: first, last, off             
                                                                          
    REAL(8)                               :: e_init                       
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: total_gamma                  
    REAL(8), ALLOCATABLE, DIMENSION(:,:)  :: partial_gamma                
    REAL(8)                               :: gammae !result from stieltjes
    REAL(8)                               :: sum_partial                  
    REAL(8)                               :: norm


!
! Code
!

    CALL alloc(buf,intbuf)
    CALL alloc(ioi,intbuf)
    CALL alloc(ioj,intbuf)

    CALL alloc(vec_start,ladc) !to be deallocated
    CALL alloc(vec1,ladc)
    CALL alloc(vec2,ladc)

    CALL alloc(e2h1p,nr_in_evecs*nr_fin)
    CALL alloc(t_mom,nr_in_evecs*nr_fin)

    CALL alloc(total_gamma,nr_in_evecs)                                   
    CALL alloc(partial_gamma,reladc_fano_nrgroups,nr_in_evecs)

    WRITE(iw,*) 'Number of Initial States: ', nr_in_evecs

    vec_start = 0

    OPEN(fhmat,FILE=fn_adcmat,FORM='unformatted',STATUS='unknown')
!    OPEN(fhdiag,FILE=fn_adcdiag,FORM='unformatted',STATUS='unknown')

    DO i = 1, nr_in_evecs
      counter = 0
      vec1 = 0
      vec_start(1:nr_in) = cin_evecs(1:nr_in,i)
!      WRITE(iw,*) 'starting vector'
!      WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)

! Choose file name of the pseudo-spectrum                                 
      IF(i.GT.9) THEN                                                     
        WRITE(fn_trans,'(a3,a1,i2)') nmspec,'.',i                         
      ELSE                                                                
        WRITE(fn_trans,'(a3,a2,i1)') nmspec,'.0',i                        
      END IF                                                              
! Choose file name for Stieltjes output                                   
      IF(i.GT.9) THEN                                                     
        WRITE(fn_stieltjes,'(a3,a1,i2)') nmstie,'.',i                     
      ELSE
        WRITE(fn_stieltjes,'(a3,a2,i1)') nmstie,'.0',i
      END IF

!      CALL alloc(diag,ladc)
!      REWIND(fhdiag)
      REWIND(fhmat)
!      READ(fhdiag) (diag(ixx), ixx=1, ladc), nbufs
!      DO j = 1, ladc
!        vec1(j) = vec1(j) + diag(j) * vec_start(j)
!      END DO
!      CALL dealloc(diag)
!      !WRITE(iw,*) 'Vector with diagonal'
!      !WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
!      WRITE(iw,*) 'nbufs =', nbufs

      DO k = 1, nbufs
        READ(fhmat) (buf(ixx),ixx=1,intbuf), &
                    (ioi(ixx),ixx=1,intbuf), &
                    (ioj(ixx),ixx=1,intbuf), &
                    nact, idummy
!        WRITE(iw,*) 'nact = ', nact
        DO j = 1, nact
          irow   = ioi(j)
          jcol   = ioj(j)
          buffer = buf(j) ! lower triangle is stored
!          WRITE(iw,*) 'buffer = ', buffer

          vec1(irow) = vec1(irow) + CONJG(buffer) * vec_start(jcol)
          vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)
        END DO
      END DO
!      WRITE(iw,*) 'Vector with total sorted matrix'
!      WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)
      
!
! Now we have the initial state vector multiplied with the sorted ADC matrix.
! We will further multiply this vector to each and every final state vector.
!
      DO j = 1, nr_fin
        vec2 = 0
        counter = counter + 1
        vec2(ladc-nr_fin+1:ladc) = cfin_evecs(1:nr_fin,j)
        !WRITE(iw,*) 'Final state vector'
        !WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)
        e2h1p(counter) = fin_energies(j)
        ! This part is ony necessary, when the initial and final states
        ! are non-orthogonal, which they are not at the moment!
        !t_mom(counter) = DOT_PRODUCT(vec1,vec2) &
        !                 - DOT_PRODUCT(vec_start,vec2) * in_energy(i)
        t_mom(counter) = 2*pi*ABS(DOT_PRODUCT(vec1,vec2)) * ABS(DOT_PRODUCT(vec1,vec2))
      END DO

      WRITE(iw,*)                                                         
      WRITE(iw,*) '-------------------------------------------------------------'       
      WRITE(iw,'(1X,A45,I4)') 'Pseudo-spectrum for initial state number ', i
      WRITE(iw,'(1X,A10,F10.6,A3,A20,F10.6)') 'energy: ', in_energy(1,i)*autoev,&       
                                              'eV', 'polestrength: ', in_energy(2,i)    
      WRITE(iw,*) '-------------------------------------------------------------'       
      WRITE(iw,*)                                                         
      WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'              
      IF (counter > 20) THEN                                              
        WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,20)    
      ELSE                                                                
        WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
      END IF

      WRITE(iw,*)                                                         
      WRITE(iw,*) 'Shifting pseudo-spectrum'                              
      FORALL (ixx = 1:counter)                                            
        e2h1p(ixx) = e2h1p(ixx) - in_energy(1,i)                          
       ! t_mom(ixx) = 2 * pi * t_mom(ixx) * t_mom(ixx)                     
      END FORALL                                                          
                                                                          
      WRITE(iw,*)                                                         
      WRITE(iw,*) 'Calling Stieltjes.'                                    
      !e_init = in_energy(1,i)                                            
      e_init = 0.0                                                        
#ifdef HAS_STIELTJES
      CALL stieltjes(counter,e2h1p,t_mom,e_init,fn_stieltjes,&            
                     gammae,2004,0)                                       
#else
      call quit('stieltjes not supported in this version')
#endif
!                                                                         
! printflag 1 for no external output(dirac.out), 2 for output             
!                                                                         
                                                                          
      WRITE(iw,*) 'Gamma = ', gammae*autoev, 'eV'                         
      total_gamma(i) = gammae * autoev                                    
!!                                                                        
!!                                                                        
      fh_trans = fhdiag                                                   
      OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    & 
           ACTION='WRITE',iostat=ierror)                                  
      WRITE(fh_trans,*) in_energy(1,i)                                    
      WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
      CLOSE(fh_trans)

    END DO !loop over initial states

!                                                                         
! Starting the partial decay widths                                       
!                                                                         
!                                                                         
    IF (do_partial) THEN                                                  
      WRITE(iw,*)                                                         
      WRITE(iw,*) ' Starting creation of partial decay width pseudo-spectra'
      WRITE(iw,*)                                                         
                                                                          
      DO i = 1, nr_in_evecs                                               
        last = nr_in                                                      
        vec_start(1:nr_in) = in_evecs(1:nr_in,i)                          
!        WRITE(iw,*) 'starting vector'                                    
!        WRITE(iw,'(1X,ES14.7)') (vec_start(ixx), ixx=1,ladc)             
                                                                          
        DO channel = 1, reladc_fano_nrgroups                              
          counter = 0                                                     
          vec1  = 0                                                       
          first = last + 1                                                
          last  = first + nr2h1p(channel) - 1                             
          !WRITE(iw,*) 'first = ', first, ' last = ', last                
                                                                          
          REWIND(fhmat)                                                   
                                                                          
! Ch  oose file name of the pseudo-spectrum                               
          IF(i.GT.9) THEN                                                 
            WRITE(fn_trans,'(a3,a1,i2,a1,a4)') nmspec,'.',i,'.',reladc_fano_labels(channel)
          ELSE                                                            
            WRITE(fn_trans,'(a3,a2,i1,a1,a4)') nmspec,'.0',i,'.',reladc_fano_labels(channel)
          END IF
! Choose file name of Stieltjes output                                    
          IF(i.GT.9) THEN                                                 
            WRITE(fn_trans,'(a3,a1,i2,a1,a4)') nmstie,'.',i,'.',reladc_fano_labels(channel)
          ELSE                                                            
            WRITE(fn_trans,'(a3,a2,i1,a1,a4)') nmstie,'.0',i,'.',reladc_fano_labels(channel)
          END IF                                                          
                                                                          
                                                                          
          DO k = 1, nbufs                                                 
            READ(fhmat) (buf(ixx),ixx=1,intbuf), &                        
                        (ioi(ixx),ixx=1,intbuf), &                        
                        (ioj(ixx),ixx=1,intbuf), &
                        nact, idummy                                      
!            WRITE(iw,*) 'nact = ', nact                                  
            DO j = 1, nact                                                
              irow   = ioi(j)                                             
              jcol   = ioj(j)                                             
              buffer = buf(j)                                             
!              WRITE(iw,*) 'buffer = ', buffer                            
                                                                          
              IF (first <= irow .AND. irow <= last) THEN                  
                vec1(irow) = vec1(irow) + CONJG(buffer) * vec_start(jcol)        
                vec1(jcol) = vec1(jcol) + buffer * vec_start(irow)        
              END IF                                                      
            END DO                                                        
          END DO                                                          
!          IF (channel == 2) THEN                                         
!            WRITE(iw,*) 'Vector with total sorted matrix'                
!            WRITE(iw,'(1X,ES14.7)') (vec1(ixx), ixx=1,ladc)              
!          END IF                                                         
                                                                          
!                                                                         
! No  w we have the initial state vector multiplied with the sorted ADC matrix.         
! We   will further multiply this vector to each and every final state vector.
!                 
          mul_fin:DO j = 1, nr_fin                                        
              vec2 = 0                                                    
              counter = counter + 1                                       
              vec2(ladc-nr_fin+1:ladc) = fin_evecs(1:nr_fin,j)            
              !IF (j == 1) THEN                                           
              !  WRITE(iw,*) 'Final state vector'                         
              !  WRITE(iw,'(1X,ES14.7)') (vec2(ixx), ixx=1,ladc)          
              !END IF                                                     
              e2h1p(counter) = fin_energies(j)                            
              t_mom(counter) = 2*pi*ABS(DOT_PRODUCT(vec1,vec2))*ABS(DOT_PRODUCT(vec1,vec2))                     
!              WRITE(iw,*) 't_mom: ', t_mom(counter)                      
          END DO mul_fin                                                  
                                                                          
! This if clause has to be out around the multiplication of final         
! states in case of an selection of final states                          
!            IF (channel == assign_fin_channel(j)) THEN                   
                                                                          
          WRITE(iw,*)                                                     
          WRITE(iw,*) 'name: ', fn_trans                                  
          WRITE(iw,*) '-------------------------------------------------------------'   
          WRITE(iw,'(1X,A45,I4)') 'Pseudo-spectrum for initial state number ', i        
          WRITE(iw,'(1X,A10,F10.6,A3,A20,F10.6)') 'energy: ', in_energy(1,i)*autoev,&   
                                                  'eV', 'polestrength: ', in_energy(2,i)
          WRITE(iw,*) '-------------------------------------------------------------'   
          WRITE(iw,*)                                                     
          WRITE(iw,'(1X,2A25)') 'Energies', 'Transition moments'          
          IF (counter > 20) THEN                                          
            WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,20)
          ELSE                                                            
            WRITE(iw,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
          END IF                                                          
          WRITE(iw,*)
!                                                                         
! This procedure only makes sense, when those final vectors not           
! associated with the chosen final states are not considered in           
! the calculation of the transition moments.                              
!                                                                         
          FORALL (ixx = 1:counter)                                        
            e2h1p(ixx) = e2h1p(ixx) - in_energy(1,i)                      
          END FORALL                                                      
                                                                          
          WRITE(iw,*) 'Channel #', channel                                
          WRITE(iw,*)                                                     
          WRITE(iw,*) 'Calling Stieltjes.'                                
          WRITE(iw,*)                                                     
          WRITE(iw,*)                                                     
          !e_init = in_energy(1,i)                                        
          e_init = 0.0                                                    
#ifdef HAS_STIELTJES
          CALL stieltjes(counter,e2h1p,t_mom,e_init,fn_stieltjes,&        
                         gammae,1000+channel,0)                           
#else
      call quit('stieltjes not supported in this version')
#endif
!                                                                         
! printflag 1 for no external output(dirac.out), 2 for output             
!                                                                         
          WRITE(iw,*) 'Gamma = ', gammae*autoev, 'eV'                     
          partial_gamma(channel,i) = gammae*autoev                        
!!                                                                        
!!                                                                        
          OPEN(fh_trans,FILE=fn_trans,STATUS='UNKNOWN',FORM='FORMATTED',    &
               ACTION='WRITE',iostat=ierror)                              
          WRITE(fh_trans,*) in_energy(1,i)                                
          WRITE(fh_trans,'(1X,2ES25.15)') (e2h1p(ixx), t_mom(ixx), ixx=1,counter)
          CLOSE(fh_trans)                                                 
                                                                          
        END DO ! end of loop over channels                                
      END DO   ! end of loop over initial states                          
                                                                          
    END IF 


    CALL dealloc(buf)
    CALL dealloc(ioi)
    CALL dealloc(ioj)

    CALL dealloc(vec_start)
    CALL dealloc(vec1)
    CALL dealloc(vec2)

    CALL dealloc(in_evecs)
    CALL dealloc(fin_evecs)
    CALL dealloc(fin_energies)

    CALL dealloc(e2h1p)
    CALL dealloc(t_mom)

!                                                                         
! Beautiful Summary of the calculated decay widths                        
!                                                                         
                                                                          
    WRITE(iw,*)                                                           
    WRITE(iw,*)                                                           
                                                                          
    DO i = 1, nr_in_evecs                                                 
      sum_partial = 0.0                                                   
      DO j = 1, reladc_fano_nrgroups                                      
        sum_partial = sum_partial + partial_gamma(j,i)                    
      END DO                                                              
      norm = sum_partial / total_gamma(i)                                 
      partial_gamma(1:reladc_fano_nrgroups,i) = partial_gamma(1:reladc_fano_nrgroups,i) &
                                                / norm                    
      WRITE(iw,149) 'Initial state # ', i, 'energy: ', in_energy(1,i)*autoev, 'eV', &   
                  'polestr.: ', in_energy(2,i)                            
      WRITE(iw,*) '--------------------------------------------------------------'      
      WRITE(iw,*) 'Total decay width: ', total_gamma(i), 'eV'             
      WRITE(iw,*) 'Sum of partial decay widths: ', sum_partial, 'eV'      
      WRITE(iw,*) 'Renormalized partial decay widths in eV:'              
      WRITE(iw,*) '--------------------------------------------------------------'      
      DO channel = 1, reladc_fano_nrgroups                                
        WRITE(iw,150) 'Channel #', channel, partial_gamma(channel,i)      
      END DO                                                              
      WRITE(iw,*) '--------------------------------------------------------------'
      WRITE(iw,*)                                                         
      WRITE(iw,*)                                                         
    END DO                                                                
                                                                          
    149 FORMAT(' ',A16,2X,I3,3X,A8,2X,F14.5,1X,A2,3X,A10,2X,F7.5)         
    150 FORMAT(' ',A9,2X,I3,3X,ES14.5)                                    
                                                                          
    CALL dealloc(total_gamma)                                             
    CALL dealloc(partial_gamma)                                           
    CALL dealloc(in_energy)

  END SUBROUTINE cfano_create_t_mom
END MODULE adc_fano_matmul
