MODULE adc_fano_diag

CONTAINS

  SUBROUTINE fano_wcondat_s(LUN,FILENAME,IRECL,DESREP,lenmat, &
                            finalpos,fin_max,iw)

    IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!   Writes explicit irep names and spinor numbers for the various
!   final state configurations belonging to symmetry DESREP.
!   Avoiding too complicated things we fix IRECL in this routine
!   no matter which value is transferred to this routine. Therefore
!   IRECL is an OUTPUT parameter.
!   Configurations of the final states are kicked out, they won't
!   won't be used in FanoADC.
!
!---------------Calling variables--------------------------------------
!
    CHARACTER(6), INTENT(IN)                 :: FILENAME
    INTEGER, INTENT(IN)                     :: LUN,DESREP,lenmat
    INTEGER, INTENT(OUT)                    :: irecl
    INTEGER, INTENT(IN)                     :: fin_max
    INTEGER, INTENT(IN), DIMENSION(fin_max) :: finalpos
    INTEGER, INTENT(IN)                     :: iw
!
!---------------Common Blocks--------------------------------------
!
#include "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!                                                                 
    INTEGER       :: lcount                                              
    INTEGER       :: satcount, fincount
    CHARACTER*44  :: FIELD44                                        
    INTEGER       :: bl_conf                                             
    CHARACTER*4   :: bl_irna,bl_abs                                  

    INTEGER       :: a, k , l
    INTEGER       :: arep, krep, lrep, klrep
    INTEGER       :: aoff, koff, loff
    INTEGER       :: kmin
                                                                  
!                                                                 
!---------------Executable code--------------------------------------
!                                                                 
!  Code for ** single ionization **:  1h/2h1p                     
!                                                                 
! fix record length for the direct access file and set some       
! blank variables. Then open file and initialize configuration    
! counter.                                                        
! *ATT* Is has turned out advantageous that storing the *absolute*
! spinor numbers is superior to storing the relative ones for the 
! subsequent postprocessing!                                      
!                                                                 
! this will be followed by the single ionization case for the moment
! (1h and 2h1p) configurations                       
!                                                                 
    IRECL = 44                                                  
    BL_CONF = 0                                                 
    BL_IRNA = '    '                                            
    BL_ABS  = 'abs.'                                            
    OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,    &
   &     STATUS='UNKNOWN')                                      
    LCOUNT = 1                                                  
    satcount = 0
    !fincount = 1
                                                                
!_________________________ SIPS ______________________________    
!|                                                                
!|                                                                
!|                                                                
      WRITE(iw,*)
      WRITE(IW,*) 'Writing FanoSIP configs to ',FILENAME            
      WRITE(IW,*) 'Final state symmetry: ',DESREP               
      WRITE(iw,*)
                                                                  
! ** start with 1h configurations                                 
                                                                  
      DO k=1,no(DESREP)                                         
        koff=io(DESREP)+k                                       
        WRITE(FIELD44,381)                                &      
                k,repna(DESREP),                          &      
                koff,repna(DESREP),                       &      
                bl_conf,bl_irna,                          &      
                bl_conf,bl_irna                                 
        WRITE(LUN,REC=lcount) FIELD44                           
        lcount = lcount + 1                                     
      ENDDO                                                     
                                                                
! ** continue with 2h1p configurations                            
! ** write absolute spinor numbers                                
                                                                  
      klloop:DO  KLREP=1,NREP                                        
        AREP=MULTB(DESREP,KLREP+NREP,2)                         
        lloop:DO LREP=1,NREP                                       
          KREP=MULTB(LREP,KLREP+NREP,2)                          
          IF(KREP.LT.LREP) CYCLE lloop                           
          DO L=1,NO(LREP)                                        
            LOFF=IO(LREP)+L                                      
            KMIN=1                                               
            IF(KREP.EQ.LREP) KMIN = L + 1                        
            DO K=KMIN,NO(KREP)                                   
              KOFF=IO(KREP)+K                                    
              aloop:DO A=1,NV(AREP)                                    
                AOFF=IO(NREP+1)+IV(AREP)+A                       
                satcount = satcount + 1
                IF (ANY(finalpos == satcount)) THEN
                  CYCLE aloop
                END IF
!                WRITE(iw,*) 'satcount = ', satcount
!                 WRITE(FIELD44,381)                               
!    &              K,REPNA(KREP),                                 
!    &              L,REPNA(LREP),                                 
!    &              A,REPNA(AREP),                                 
!    &              BL_CONF,BL_IRNA                                
                WRITE(FIELD44,381)                         &      
   &              KOFF,REPNA(KREP),                        &      
   &              LOFF,REPNA(LREP),                        &      
   &              AOFF,REPNA(AREP),                        &      
   &              BL_CONF,BL_IRNA                                
                WRITE(LUN,REC=lcount) FIELD44                    
                LCOUNT=LCOUNT+1                                  
              END DO aloop                                             
            ENDDO                                                
          ENDDO                                                  
        END DO lloop
      END DO klloop
                                                                  
! ** file is open only within this subroutine                     
                                                                  
      CLOSE(LUN)                                                
                                                                  
! ** correct configuration counter                                
                                                                  
      LCOUNT = LCOUNT - 1                                       
      IF(LCOUNT.NE.lenmat) THEN                                
        WRITE(IW,*) 'LCOUNT,lenmat:',LCOUNT,lenmat            
        CALL QUIT('Error in fano_wcondat_s for single ionizations')  
      ENDIF                                                     
!|                                                                
!|                                                                
!|____________________________________________________________    
                                                                  
  381  FORMAT(4(I3,'(',A4,')  '))                                  
                                                                  
  END SUBROUTINE fano_wcondat_s             


!
!--------------------------------------------------------------
!
  SUBROUTINE fano_find_in_vec(n,evl,z,positions,nr_in_vecs,iw)
!
!  Purpose: To find the eigenvectors of the Lanczos matrix having
!           a non-negliable contribution to the Fano initial state
!           and write the positions within Z to an array

    use adc_cfg
    use memory_allocator
    use adc_fano_exchange

    IMPLICIT NONE
!
! Calling Variables
!
    INTEGER, INTENT(IN)                 :: n
    REAL(8), INTENT(IN), DIMENSION(N)   :: evl
    REAL(8), INTENT(IN), DIMENSION(N,N) :: z
    INTEGER, INTENT(OUT), DIMENSION(N)  :: positions
    INTEGER, INTENT(OUT)                :: nr_in_vecs
    INTEGER, INTENT(IN)                 :: iw

!
! Local Variables
!
    REAL(8)                          :: contrib
    REAL, PARAMETER                  :: thresh=0.15
    REAL(8), PARAMETER               :: ovlthresh=1.0D-2
    REAL(8), PARAMETER               :: ethresh=1.0D-12
    INTEGER                          :: initial
    INTEGER                          :: i,j,k
    INTEGER,ALLOCATABLE,DIMENSION(:) :: selected
    INTEGER                          :: nr_in
    REAL(8),ALLOCATABLE,DIMENSION(:) :: vec1, vec2 
    REAL(8)                          :: ovl
    REAL(8)                          :: energy1, energy2

!
! Code
!
    positions  = 0
    nr_in      = 0
    nr_in_vecs = 0
    initial    = reladc_fano_inrelsp

    CALL alloc(selected,n,id="selected eigenvectors in Fano Lanczos")
    selected   = 0

    DO j=1,n
      contrib = z(initial,j) * z(initial,j)
      IF (contrib >= thresh) THEN
        nr_in      = nr_in + 1
        selected(nr_in) = j
      END IF
    END DO

    SELECT CASE (reladc_fano_ovl_in)
      CASE (.TRUE.)
        !CALL alloc(vec1,n)
        !CALL alloc(vec2,n)
        !DO i = 1, nr_in
        !  vec1(1:n) = z(1:n,selected(i))
        !  DO j = 1, n
        !    vec2(1:n) = z(1:n,j)
        !    ovl = dot_product(vec1, vec2)
        !    WRITE(iw,*) 'j = ', j
        !    WRITE(iw,*) 'Overlap = ', ovl
        !    IF (1.0D0 - ABS(ovl) < ovlthresh) THEN
        !      nr_in_vecs = nr_in_vecs + 1
        !      positions(nr_in_vecs) = j
        !    END IF
        !  END DO
        !END DO
        !CALL dealloc(vec1)
        !CALL dealloc(vec2)
        WRITE(iw,*) 'Selection by overlap does not work yet.'

      CASE (.FALSE.)
        DO i = 1, nr_in
          energy1      = evl(selected(i))
          DO j = 1, n
            energy2 = evl(j)
            IF (ABS(energy1 - energy2) < ethresh) THEN
              nr_in_vecs = nr_in_vecs + 1
              positions(nr_in_vecs) = j
            END IF
          END DO
        END DO
    END SELECT

    CALL dealloc(selected)

    WRITE(iw,*)
    WRITE(iw,*) '--------------------------------------------------------------'
    WRITE(iw,*) 'Selection of Lanczos eigenvectors as Fano initial state(s)'
    WRITE(iw,*) '--------------------------------------------------------------'
    WRITE(iw,*) 'Energy threshold: ', ethresh
    WRITE(iw,*) 'Overlap threshold: ', ovlthresh
    WRITE(iw,*) 'Selection method overlap(T), energy(F): ', reladc_fano_ovl_in
    WRITE(iw,*) 'Number of contributions: ', nr_in
    WRITE(iw,*) 'Number of selected eigenvectors: ', nr_in_vecs
    WRITE(iw,*) 'Selected positions:'
    WRITE(iw,'(1X,5I4)') (positions(i), i=1, nr_in_vecs)
    WRITE(iw,*)


  END SUBROUTINE fano_find_in_vec




!
!--------------------------------------------------------------
!
  SUBROUTINE cfano_find_in_vec(n,evl,z,positions,nr_in_vecs,iw)
!
!  Purpose: To find the eigenvectors of the Lanczos matrix having
!           a non-negliable contribution to the Fano initial state
!           and write the positions within Z to an array

    use adc_cfg
    use memory_allocator
    use adc_fano_exchange

    IMPLICIT NONE
!
! Calling Variables
!
    INTEGER, INTENT(IN)                    :: n
    REAL(8), INTENT(IN), DIMENSION(N)      :: evl
    COMPLEX(8), INTENT(IN), DIMENSION(N,N) :: z
    INTEGER, INTENT(OUT), DIMENSION(N)     :: positions
    INTEGER, INTENT(OUT)                   :: nr_in_vecs
    INTEGER, INTENT(IN)                    :: iw

!
! Local Variables
!
    REAL(8)                             :: contrib
    REAL, PARAMETER                     :: thresh=0.05
    REAL(8), PARAMETER                  :: ovlthresh=1.0D-2
    REAL(8), PARAMETER                  :: ethresh=1.0D-12
    INTEGER                             :: initial
    INTEGER                             :: i,j,k
    INTEGER,ALLOCATABLE,DIMENSION(:)    :: selected
    INTEGER                             :: nr_in
    COMPLEX(8),ALLOCATABLE,DIMENSION(:) :: vec1, vec2 
    REAL(8)                             :: ovl
    REAL(8)                             :: energy1, energy2

!
! Code
!
    positions  = 0
    nr_in      = 0
    nr_in_vecs = 0
    initial    = reladc_fano_inrelsp

    CALL alloc(selected,n,id="selected eigenvectors in Fano Lanczos")
    selected   = 0

    DO j=1,n
      contrib = CDABS(z(initial,j)) * CDABS(z(initial,j))
      IF (contrib >= thresh) THEN
        nr_in      = nr_in + 1
        selected(nr_in) = j
      END IF
    END DO

    SELECT CASE (reladc_fano_ovl_in)
      CASE (.TRUE.)
        !CALL alloc(vec1,n)
        !CALL alloc(vec2,n)
        !DO i = 1, nr_in
        !  vec1(1:n) = z(1:n,selected(i))
        !  DO j = 1, n
        !    vec2(1:n) = z(1:n,j)
        !    ovl = dot_product(vec1, vec2)
        !    WRITE(iw,*) 'j = ', j
        !    WRITE(iw,*) 'Overlap = ', ovl
        !    IF (1.0D0 - ABS(ovl) < ovlthresh) THEN
        !      nr_in_vecs = nr_in_vecs + 1
        !      positions(nr_in_vecs) = j
        !    END IF
        !  END DO
        !END DO
        !CALL dealloc(vec1)
        !CALL dealloc(vec2)
        WRITE(iw,*) 'Selection by overlap does not work yet.'

      CASE (.FALSE.)
        DO i = 1, nr_in
          energy1      = evl(selected(i))
          DO j = 1, n
            energy2 = evl(j)
            IF (ABS(energy1 - energy2) < ethresh) THEN
              nr_in_vecs = nr_in_vecs + 1
              positions(nr_in_vecs) = j
            END IF
          END DO
        END DO
    END SELECT

    CALL dealloc(selected)

    WRITE(iw,*)
    WRITE(iw,*) '--------------------------------------------------------------'
    WRITE(iw,*) 'Selection of Lanczos eigenvectors as Fano initial state(s)'
    WRITE(iw,*) '--------------------------------------------------------------'
    WRITE(iw,*) 'Energy threshold: ', ethresh
    WRITE(iw,*) 'Overlap threshold: ', ovlthresh
    WRITE(iw,*) 'Selection method overlap(T), energy(F): ', reladc_fano_ovl_in
    WRITE(iw,*) 'Number of contributions: ', nr_in
    WRITE(iw,*) 'Number of selected eigenvectors: ', nr_in_vecs
    WRITE(iw,*) 'Selected positions:'
    WRITE(iw,'(1X,5I4)') (positions(i), i=1, nr_in_vecs)
    WRITE(iw,*)


  END SUBROUTINE cfano_find_in_vec


END MODULE adc_fano_diag
