MODULE adc_fano_real_routines

!
!  Data dictionary: common variables and arrays for the Fano ADC
!  
  INTEGER, PARAMETER                   :: fano_intbuf = 5 * 1024 * 1024 !as in ADC
  INTEGER                              :: itapfano !same as itapadc, written in a module
                                                   !file handle basis
  INTEGER                              :: fin_max  !number of final state configs
  INTEGER                              :: ierror = 0
  
!  Variables for the big, sorted ADC matrix
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: adcbuf   !alloc in fano_write_hhp
                                                   !dealloc in fano_make_hphp
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: adcioi, adcioj !"
  INTEGER                              :: adcnbufs !number of buffers for ADC matrix
  INTEGER                              :: adcib    !book keeping number for adcbuf
  INTEGER                              :: adccol, adcrow
  INTEGER                              :: adcfh    ! file handle for ADC matrix
  INTEGER                              :: adcdiagfh

!  Variables for the intial state matrix
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: inbuf
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: inioi, inioj !"
  INTEGER                              :: innbufs
  INTEGER                              :: inib
  INTEGER                              :: incol, inrow
  INTEGER                              :: infh
  INTEGER                              :: indiagfh

!  Variables for the final state matrix
  INTEGER                              :: finnbufs
  INTEGER                              :: finfh, findiagfh


!  File names
  CHARACTER(6), PARAMETER              :: fn_indiag  = 'INDIAG'
  CHARACTER(6), PARAMETER              :: fn_inmat   = 'IN_MAT'
  CHARACTER(6), PARAMETER              :: fn_adcdiag = 'FADCDG'
  CHARACTER(6), PARAMETER              :: fn_adcmat  = 'ADCMAT'
  CHARACTER(6), PARAMETER              :: fn_finmat  = 'FINMAT'

                                               
  
  

CONTAINS

  SUBROUTINE fano_write_hhp(ckks,bufhp,ladc,crep,bigicra,itapadc,iw, &
                            finalpos)

!
!  Purpose: Calculate the coupling block of the ADC matrix and sort
!           the matrix elements according to their description of
!           initial or final states. Write out the lower triangle
!           of the matrix for the initial state description and write
!           the coupling block of the sorted ADC matrix.
!
!           -----                -----
!     1h    |   |             1h |   |
!           -----                -----
!           |   |   ------>      |   |
!     fin   ~~~~~                |   |
!           |   |            fin ~~~~~
!           -----                -----
!
!  Author: Elke Fasshauer
!

    use adc_cfg
    use memory_allocator

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)    :: crep              ! current representation
    REAL(8), INTENT(IN), DIMENSION(no(crep)*no(crep))    :: ckks
    REAL(8), INTENT(IN), DIMENSION(no(crep)*nvoot(crep)) :: bufhp
    INTEGER, INTENT(IN)    :: ladc              ! length of adc matrix
    INTEGER, DIMENSION(9,nvoot(crep)), INTENT(IN) :: bigicra
    INTEGER, INTENT(IN)    :: itapadc           ! file handling number of ADC matrix
    INTEGER, INTENT(IN)    :: iw                ! file handling number of output file
    INTEGER, INTENT(IN), DIMENSION(fin_max)              :: finalpos

!
!  Data dictionary: Local Variables
!
    REAL(8),ALLOCATABLE,DIMENSION(:)   :: single_col   !single column of ADC matrix
    REAL(8)                            :: temp
    INTEGER                            :: temp_row     
    INTEGER                            :: finrow       !row of final states
    INTEGER                            :: k, ixx, fincount ! counters
    LOGICAL                            :: is_in  = .TRUE.
    LOGICAL                            :: is_coupl = .FALSE.
    LOGICAL                            :: is_fin = .FALSE.
    INTEGER, PARAMETER                 :: jdummy = 0

!
!  Data dictionary: Variables taken from WRITE_HHP, not initialized there
!
    INTEGER                            :: off1, off2, off3
    INTEGER                            :: nok, najl
    INTEGER                            :: dummy

!
!  allocate arrays
!
    CALL alloc(single_col,ladc,id='single column in Fano hhp')
    CALL alloc(inbuf,fano_intbuf,id='buffer of initial state matrix')
    CALL alloc(inioi,fano_intbuf,id='row buffer initial state')
    CALL alloc(inioj,fano_intbuf,id='col buffer initial state')
    CALL alloc(adcbuf,fano_intbuf,id='buffer of sorted ADC matrix')
    CALL alloc(adcioi,fano_intbuf,id='row buffer sorted ADC matrix')
    CALL alloc(adcioj,fano_intbuf,id='col buffer sorted ADC matrix')

! ---------------------------------
!  Code
! ---------------------------------

!  clear output buffers
    single_col = 0.0D0
    inbuf      = 0.0D0
    inioi      = 0
    inioj      = 0
    adcbuf     = 0.0D0
    adcioi     = 0
    adcioj     = 0

    inrow      = 1
    finrow     = 0
    fincount   = 0

    itapfano   = itapadc
!  remember: at this moment, itapadc+4 and itapadc+5
!            are open due to adc_sngl_main
    adcfh      = itapfano + 7 !Files need to be opened
    infh       = itapfano + 9

    ! open the matrix files of the sorted ADC matrix and the matrix of the
    ! initial states
    OPEN(adcfh,FILE=fn_adcmat,STATUS='UNKNOWN',FORM='UNFORMATTED',  &
         ACTION='WRITE',iostat=ierror)
    OPEN(infh,FILE=fn_inmat,STATUS='UNKNOWN',FORM='UNFORMATTED',    &
         ACTION='WRITE',iostat=ierror)


!                                                                       
!  create h/h and h/2h1p column by column
!  write out h/h and h/2h1p buffers                                     
!
    nok  = no(crep)
    WRITE(iw,*) 'nok = ', nok
    najl = nvoot(crep)

    inib  = 0
    adcib = 0

    off1 = 1
    off3 = 1

    one_column:DO k = 1, nok
      adccol = k
      off2   = 1
      CALL XCOPY(nok,ckks(off1),1,single_col(off2),1)
      off2   = off2 + nok
      CALL XCOPY(najl,bufhp(off3),1,single_col(off2),1)

      rowloop:DO adcrow = adccol + 1, ladc
        incol = adccol
        IF (ANY(finalpos + nok == adcrow)) THEN
          !WRITE(iw,*) 'adcrow hhp = ', adcrow
          CYCLE rowloop

        ELSE
          inrow = inrow + 1
          !WRITE(iw,*) 'inrow = ', inrow
          IF (single_col(adcrow) /= 0.0D0) THEN
            ! Handle the matrix of initial states
            inib   = inib + 1
            inbuf(inib) = single_col(adcrow)
            inioi(inib) = inrow
            inioj(inib) = adccol
!            WRITE(iw,*) 'inioj = ', inioj(inib), 'inioi = ', inioi(inib)
            IF (inib == fano_intbuf) THEN
              CALL rwrite_to_disk(inbuf,inioi,inioj,fano_intbuf, &
                                  innbufs,inib,infh)
            END IF
          END IF
        END IF
      END DO rowloop
!!
!!  One column has been processed, now the final states are added to
!!  the end of the corresponding column
!!
      DO temp_row = 1, fin_max
        temp = single_col(finalpos(temp_row) + nok)
        !WRITE(iw,*) 'processing: ', finalpos(temp_row)
        IF (temp /= 0.0D0) THEN
          adcib = adcib + 1
          adcbuf(adcib) = temp
          adcioi(adcib) = temp_row + inrow
          adcioj(adcib) = adccol
          IF (adcib == fano_intbuf) THEN
            CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf,adcnbufs,  &
                  adcib,adcfh)
          END IF
        END IF
      END DO 

      inrow = adccol + 1
      finrow = 0
      fincount = 0

      off1 = off1 + nok
      off3 = off3 + najl

    END DO one_column

!    DO ixx=1, najl + nok
!      WRITE(iw,'(1X,2(ES14.5,I4,I4))') -inbuf(ixx),inioi(ixx),inioj(ixx),&
!                                     -adcbuf(ixx),adcioi(ixx),adcioj(ixx)
!    END DO

!
!  deallocate local arrays
!
    CALL dealloc(single_col)


  END SUBROUTINE fano_write_hhp




!
!---------------------------------------------------------------------
!
  SUBROUTINE fano_makehphp(ckks,eajl,oooo,vovo,bigicra,ladc,krep,  &
                           finalpos,iw,mxno)

!
!  Purpose: Create all elements of the ADC Sat Block and write them to the
!           different Fano matrices
!

    use adc_cfg
    use memory_allocator

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER,INTENT(IN)                               :: krep
    REAL(8),INTENT(IN),DIMENSION(no(krep)*no(krep))  :: ckks
    REAL(8),INTENT(IN)                               :: eajl(*)
    REAL(8),INTENT(IN)                               :: oooo(*), vovo(*)
    INTEGER,INTENT(IN),DIMENSION(9,nvoot(krep))      :: bigicra
    INTEGER,INTENT(IN)                               :: ladc
    INTEGER,INTENT(IN),DIMENSION(fin_max)            :: finalpos
    INTEGER,INTENT(IN)                               :: iw
    INTEGER,INTENT(IN)                               :: mxno
    

!
!  Data dictionary: Local Variables
!
    INTEGER                                :: nok, najl
    INTEGER                                :: temp
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: tempe !temporary array for diag
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: temp_fin
    INTEGER                                :: irep
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:)   :: oolt
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: diag
    INTEGER                                :: fincount
    INTEGER                                :: finrow, fincol
    REAL(8),ALLOCATABLE,DIMENSION(:)       :: finbuf
    INTEGER,ALLOCATABLE,DIMENSION(:)       :: finioi, finioj
    INTEGER                                :: finib
    INTEGER                                :: irow, jcol, i, j, k
    INTEGER                                :: coupcol, couprow
    INTEGER                                :: finalrow, finalcol
    INTEGER                                :: counter
    INTEGER                                :: off1, ioff 
    REAL(8)                                :: s
    INTEGER                                :: zeros, ixx

!
!  Code
!
    WRITE(iw,*) 'FanoADC: Entering ADC Satellite block'
    nok  = no(krep) ! /= mxno
    najl = nvoot(krep)

!
! Open all files used only in this routine
!
    adcdiagfh = itapfano + 6
    indiagfh  = itapfano + 8
    finfh     = itapfano + 10

    OPEN(adcdiagfh,FILE=fn_adcdiag,STATUS='UNKNOWN',FORM='UNFORMATTED',  &
         ACTION='WRITE',iostat=ierror)
    OPEN(indiagfh,FILE=fn_indiag,STATUS='UNKNOWN',FORM='UNFORMATTED',  &
         ACTION='WRITE',iostat=ierror)
    OPEN(finfh,FILE=fn_finmat,STATUS='UNKNOWN',FORM='UNFORMATTED',  &
         ACTION='WRITE',iostat=ierror)


    CALL alloc(diag,nok+najl, id='FanoADC: Array of diagonal')

    ! Fill in the diagonal of the h/h block
    off1 = 1
    DO k = 1, nok
      diag(k) = ckks(off1 + k - 1)
      off1 = off1 + nok
    END DO
    
    ! Take care of the diagonal of the 2h1p/2h1p block
    CALL alloc(tempe,najl, id='FanoADC: temporary array of diagonal elements')
    CALL alloc(temp_fin, fin_max)
    tempe = 0
    counter = 0
    ! add the 2h1p diag sorted elements
    diag_2h1p:DO k = 1, najl
      IF (ANY(finalpos == k)) THEN
        CYCLE diag_2h1p
      END IF  
      counter = counter + 1
      tempe(counter) = eajl(k)
    END DO diag_2h1p

    DO fincount = 1, fin_max
      temp = finalpos(fincount)
      temp_fin(fincount) = eajl(temp)
    END DO

    WRITE(iw,*) 'Diagonal sorted.'

    ! add final state diagonal
    tempe(najl-fin_max+1:najl) = temp_fin(1:fin_max)
    diag(nok+1:ladc)           = tempe(1:najl)

    CALL dealloc(temp_fin)
    CALL dealloc(tempe)

    WRITE(iw,*) 'Diagonal is written to temporary array'

    ! In case of strict ADC2 the offdiagonal elements are zero
    ! and no sorting and/or writing is necessary.
    offdiag:IF(reladc_adclevel >= 2) THEN
      CALL GETOOOO(oooo)
      CALL GETVOVO(vovo)

      ! construct oolt lookup table for KCOCC
      ALLOCATE(oolt(mxno,mxno,nrep))
      oolt = 0
      DO irep = 1, nrep
        ioff = 1
        DO j = 1, no(irep)
          DO i = j + 1, no(irep)
            oolt(i,j,irep) = ioff
            ioff = ioff + 1
          END DO
        END DO
      END DO

      WRITE(iw,*) 'OOLT constructed'

!  
!    loop over the original structure of the ADC matrix and write the matrix
!    elements sorted into different buffers
!  
      fincol = 1
      zeros = 0

      colloop:DO jcol = 1, najl !the loop is over the original ADC matrix cols
        IF (ANY(finalpos == jcol)) THEN
          coupcol = adccol
!     couplings between final states that are not taken care of in the column
!     wise handling are taken care of here (row < col).
          couploop:DO irow = jcol+1, najl 
            IF (ANY(finalpos == irow)) THEN
              CYCLE couploop
            END IF

            s = fano_kcocc(irow,jcol,oooo,vovo,bigicra,najl,oolt,mxno)

            det_row:DO i = 1, fin_max
              IF (finalpos(i) == jcol) THEN
                couprow = ladc - fin_max + i
                EXIT det_row
              END IF
            END DO det_row
            !couprow = ladc - fin_max + fincol
            coupcol = coupcol + 1

            IF (s /= 0.0D0) THEN
              adcib = adcib + 1
              adcbuf(adcib) = s
              adcioi(adcib) = couprow
              adcioj(adcib) = coupcol
              IF (adcib == fano_intbuf) THEN
                CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf, &
                                    adcnbufs,adcib,adcfh)
              END IF
            END IF
            
          END DO couploop
          fincol  = fincol + 1
          CYCLE colloop
        END IF
         
        ! determine the column of the initial and total matrix
        adccol = adccol + 1
        incol  = incol  + 1

        finrow = 1
        adcrow = adccol - 1
        inrow  = incol - 1

        rowloop:DO irow = jcol, najl
          IF (ANY(finalpos == irow)) THEN
            CYCLE rowloop
          END IF

          ! determine the row of the initial and total matrix
          adcrow = adcrow + 1
          inrow  = inrow  + 1

          s = fano_kcocc(irow,jcol,oooo,vovo,bigicra,najl,oolt,mxno)
          !WRITE(iw,'(1X,A4,ES15.4)') 's = ', s
          IF (s == 0.0D0) THEN
            zeros = zeros + 1
          END IF

          isdiag:IF (irow == jcol) THEN
            diag(adcrow) = diag(adcrow) + s
          ELSE isdiag
            IF (s /= 0.0D0) THEN
              inib = inib + 1
              inbuf(inib) = s
              inioi(inib) = inrow
              inioj(inib) = incol
              IF (inib == fano_intbuf) THEN
                CALL rwrite_to_disk(inbuf,inioi,inioj,fano_intbuf, &
                                    innbufs,inib,infh)
              END IF
              ! now total adc-matrix
              ! we skipt that, since we only need the interaction matrix
              !adcib = adcib + 1
              !adcbuf(adcib) = s
              !adcioi(adcib) = adcrow
              !adcioj(adcib) = adccol
              !IF (adcib == fano_intbuf) THEN
              !  CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf, &
              !                      adcnbufs,adcib,adcfh)
              !END IF
            END IF
          END IF isdiag

        END DO rowloop
!   construct coupling part between initial and final states
        DO fincount = 1, fin_max
          irow = finalpos(fincount)
          adcrow = adcrow + 1

          s = fano_kcocc(irow,jcol,oooo,vovo,bigicra,najl,oolt,mxno)
          IF (s /= 0.0D0) THEN
            adcib = adcib + 1
            adcbuf(adcib) = s
            adcioi(adcib) = adcrow
            adcioj(adcib) = adccol
            IF (adcib == fano_intbuf) THEN
              CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf, &
                                  adcnbufs,adcib,adcfh)
            END IF
          ELSE
            zeros = zeros + 1
          END IF
        END DO

      END DO colloop

      WRITE(iw,*) 'Number of zeros = ', zeros

      CALL rwrite_rest_to_disk(inbuf,inioi,inioj,fano_intbuf,innbufs, &
                               inib,infh)

      CALL dealloc(inbuf)
      CALL dealloc(inioi)
      CALL dealloc(inioj)
      

      ! take care of the block of final states
      fincol = 0
      finib = 0
      finnbufs = 0
      CALL alloc(finbuf, fano_intbuf,id='buffer of final state')
      CALL alloc(finioi, fano_intbuf,id='row buffer final state')
      CALL alloc(finioj, fano_intbuf,id='col buffer final state')
      DO finalcol = 1, fin_max
        adccol = adccol + 1
        fincol = fincol + 1
        jcol   = finalpos(finalcol)
        adcrow = adccol - 1
        finrow = fincol - 1
        DO finalrow = finalcol, fin_max
          adcrow = adcrow + 1
          finrow = finrow + 1
          irow   = finalpos(finalrow)

          s = fano_kcocc(irow,jcol,oooo,vovo,bigicra,najl,oolt,mxno)
          !WRITE(iw,'(1X,A4,ES15.4)') 's = ', s
          isfindiag:IF (adcrow == adccol) THEN
            diag(adcrow) = diag(adcrow) + s
            finib = finib + 1
            finbuf(finib) = diag(adcrow)
            finioi(finib) = finrow
            finioj(finib) = fincol
          ELSE isfindiag
            IF (s /= 0.0D0) THEN
              ! only save interaction block
              !adcib = adcib + 1
              !adcbuf(adcib) = s
              !adcioi(adcib) = adcrow
              !adcioj(adcib) = adccol
              !!WRITE(iw,*) 'adccolsat = ', adccol
              !!WRITE(iw,*) 'adcrowsat = ', adcrow
              !IF (adcib == fano_intbuf) THEN
              !  CALL rwrite_to_disk(adcbuf,adcioi,adcioj,fano_intbuf, &
              !                      adcnbufs,adcib,adcfh)
              !END IF
              finib = finib + 1
              finbuf(finib) = s
              finioi(finib) = finrow
              finioj(finib) = fincol
              !WRITE(iw,*) 'fincol = ', fincol
              !WRITE(iw,*) 'finrow = ', finrow
            END IF
          END IF isfindiag


          IF (finib == fano_intbuf) THEN
            CALL rwrite_to_disk(finbuf,finioi,finioj,fano_intbuf, &
                                finnbufs,finib,finfh)
          END IF
        END DO
      END DO

      WRITE(iw,*) 'ADCib = ', adcib
!      WRITE(iw,*) 'Finib = ', finib

      DEALLOCATE(oolt)    

      WRITE(iw,*) 'Writing final state contributions for level>1'
      CALL rwrite_rest_to_disk(adcbuf,adcioi,adcioj,fano_intbuf,adcnbufs, &
                               adcib,adcfh)

      CALL dealloc(adcbuf)
      CALL dealloc(adcioi)
      CALL dealloc(adcioj)

      CALL rwrite_rest_to_disk(finbuf,finioi,finioj,fano_intbuf,finnbufs, &
                               finib,finfh)

    ELSE offdiag
      ! deallocate all no longer used arrays
      CALL rwrite_rest_to_disk(inbuf,inioi,inioj,fano_intbuf,innbufs, &
                               inib,infh)
      CALL rwrite_rest_to_disk(adcbuf,adcioi,adcioj,fano_intbuf,adcnbufs, &
                               adcib,adcfh)
      CALL dealloc(inbuf)
      CALL dealloc(inioi)
      CALL dealloc(inioj)
      CALL dealloc(adcbuf)
      CALL dealloc(adcioi)
      CALL dealloc(adcioj)
      ! allocate final buffers
      CALL alloc(finbuf, fano_intbuf,id='buffer final state')
      CALL alloc(finioi, fano_intbuf,id='row buffer final state')
      CALL alloc(finioj, fano_intbuf,id='col buffer final state')
      finib = 0
      !WRITE(iw,*) 'finib = ', finib
      !WRITE(iw,*) 'end   = ', ladc-fin_max+1
      !WRITE(iw,*) 'ladc  = ', ladc
      DO fincount = ladc-fin_max+1,ladc
        finib = finib + 1
        finbuf(finib) = diag(fincount)
        finioi(finib) = finib
        finioj(finib) = finib
        IF (finib == fano_intbuf) THEN
          CALL rwrite_to_disk(finbuf,finioi,finioj,fano_intbuf, &
                              finnbufs,finib,finfh)
        END IF
      END DO
      CALL rwrite_rest_to_disk(finbuf,finioi,finioj,fano_intbuf,finnbufs, &
                               finib,finfh)

    END IF offdiag

    CALL dealloc(finbuf)
    CALL dealloc(finioi)
    CALL dealloc(finioj)

    ! write the diagonal to separate files
    !Coupling block does not contain diagonal elements
    !CALL rwrite_diag(diag,ladc,1,ladc,adcdiagfh,adcnbufs,iw)
    !WRITE(iw,'(1X,ES14.5)') (diag(fincount), fincount = 1,ladc)
    CALL rwrite_diag(diag,ladc,1,ladc-fin_max,indiagfh,innbufs,iw)
    !CALL rwrite_diag(diag,ladc,ladc-fin_max+1,ladc,findiagfh,adcnbufs,iw)

    CALL dealloc(diag)

!
!  close all files used only in this routine
!
    CLOSE(adcdiagfh)
    CLOSE(indiagfh)
    CLOSE(finfh)
    CLOSE(adcfh)
    CLOSE(infh)



  END SUBROUTINE fano_makehphp


!
!---------------------------------------------------------------------
!
  SUBROUTINE rwrite_to_disk(buf,ioi,ioj,intbuf,nbufs,ib,fh)
!
!  Purpose:  to write the three arrays buf, ioi, ioj of length intbuf
!            to disk
!
    IMPLICIT NONE
!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)                     :: intbuf, fh
    INTEGER, INTENT(INOUT)                  :: nbufs, ib
    REAL(8), DIMENSION(intbuf), INTENT(IN)  :: buf
    INTEGER, DIMENSION(intbuf), INTENT(IN)  :: ioi, ioj
!
!  Data dictionary: Local Variables
!
    INTEGER :: ixx
    INTEGER :: jdummy

!
!  Code
!

    jdummy = 0

    nbufs = nbufs + 1

    WRITE(fh) (-buf(ixx), ixx = 1, intbuf),    &
              (ioi(ixx), ixx = 1, intbuf),     &
              (ioj(ixx), ixx = 1, intbuf),     &
              intbuf, jdummy
    ib = 0

  END SUBROUTINE rwrite_to_disk


!
!--------------------------------------------------------------------
!
  SUBROUTINE rwrite_rest_to_disk(buf,ioi,ioj,intbuf,nbufs,ib,fh)
!
!  Purpose:  write the rest of the matrix to disk
!
    IMPLICIT NONE
!
!  Data dictionary: Calling variables
!
    INTEGER, INTENT(IN)                     :: intbuf, fh
    INTEGER, INTENT(INOUT)                  :: nbufs, ib
    REAL(8), DIMENSION(intbuf), INTENT(IN)  :: buf
    INTEGER, DIMENSION(intbuf), INTENT(IN)  :: ioi, ioj
!
!  Data dictionary: Local Variables
!
    INTEGER :: ixx
    INTEGER :: jdummy
    INTEGER :: iw=6

!
!  Code
!

    jdummy = 0

    nbufs = nbufs + 1

    WRITE(fh) (-buf(ixx), ixx = 1, intbuf),   &
              (ioi(ixx), ixx = 1, intbuf),    &
              (ioj(ixx), ixx = 1, intbuf),    &
              ib, jdummy
    ib = 0

  END SUBROUTINE rwrite_rest_to_disk

!
!---------------------------------------------------------------------
!
  SUBROUTINE rwrite_diag(diag,ladc,first,last,fh,nbufs,iw)
!
!  Purpose: to write the diagonal elements to disk
!
    IMPLICIT NONE
!
!  Data dictionary: Calling variables
!
    INTEGER, INTENT(IN)                 :: ladc
    INTEGER, INTENT(IN)                 :: first, last
    INTEGER, INTENT(IN)                 :: fh
    INTEGER, INTENT(IN)                 :: nbufs
    REAL(8), INTENT(IN), DIMENSION(ladc):: diag
    INTEGER, INTENT(IN)                 :: iw
!
!  Data dictionary: Local variables
!
    INTEGER :: ixx
!
!  Code
!

    !WRITE(iw,*) 'first = ', first
    !WRITE(iw,*) 'last  = ', last
    WRITE(fh) (-diag(ixx), ixx= first, last), nbufs


  END SUBROUTINE rwrite_diag
    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  FUNCTION fano_kcocc(IROW,JCOL,OOOO,VOVO,ICRA,najl,OOLT,MXNO)

     IMPLICIT NONE
!
!---------------Description--------------------------------------------
!
!     This function returns the C(1)_akl,a'k'l' element of the 2h1p/2h1p
!     block in the (N-1) space. The integrals OOOO and VOVO have to be
!     provided and ICRA/OOLT contain translation/lookup tables referring
!     to the irrep under consideration (generated by the caller).
!     The matrix is of dimension NVOOT x NVOOT and the proper range
!     of row I and column J has to be maintained by the caller !
!
!---------------Last modified------------------------------------------
!
!     Author : MP
!
!---------------Calling variables--------------------------------------
!
      INTEGER, INTENT(IN) :: najl
      INTEGER IROW,JCOL,MXNO,ICRA(9,najl),OOLT(MXNO,MXNO,*)
      REAL(8) OOOO(*),VOVO(*)
      REAL(8) fano_kcocc
!
!---------------Common Blocks--------------------------------------
!
#include  "../relccsd/symm.inc"
!
!---------------Local variables--------------------------------------
!
      REAL*8 S
      INTEGER :: A,AS,K,KS,L,LS
      INTEGER :: AREP,ASREP,KREP,KSREP,LREP,LSREP
      INTEGER :: ASKREP, ASLREP, KLREP
      INTEGER :: IOFF
      INTEGER :: tcol, trow
!
!---------------Executable code--------------------------------------
!
      IF (irow < jcol) THEN
        tcol = irow
        trow = jcol
      ELSE
        tcol = jcol
        trow = irow
      END IF

      AS    = ICRA(1,TCOL)
      ASREP = ICRA(2,TCOL)
      KS    = ICRA(3,TCOL)
      KSREP = ICRA(4,TCOL)
      LS    = ICRA(5,TCOL)
      LSREP = ICRA(6,TCOL)


      A     = ICRA(1,TROW)                                        
      AREP  = ICRA(2,TROW)
      K     = ICRA(3,TROW)
      KREP  = ICRA(4,TROW)                                        
      L     = ICRA(5,TROW)                                        
      LREP  = ICRA(6,TROW)                                        
                                                                  
      S=0.0D0                                                     
                                                                  
!                                                                 
!   ..... Part A                                                  
!                                                                 
      IF(AREP.EQ.ASREP) THEN                                      
        IF(A.EQ.AS) THEN                                          
          KLREP=MULTB(KREP,LREP,1)                                
          IF(KREP.EQ.LREP) THEN                                   
            IOFF=IIOOT(KREP,LREP) + OOLT(K,L,KREP) - 1            
          ELSE                                                    
            IOFF=IIOOT(KREP,LREP) + (L-1)*NO(KREP) + K - 1        
          END IF                                                   
          IOFF=IOFF * NOOT(KLREP) + IOOOOTT(KLREP)                
          IF(KSREP.EQ.LSREP) THEN                                 
            IOFF=IOFF + IIOOT(KSREP,LSREP) + OOLT(KS,LS,KSREP)    
          ELSE                                                    
            IOFF=IOFF + IIOOT(KSREP,LSREP) +                      &
     &        (LS-1)*NO(KSREP) + KS                               
          END IF                                                   
          S = S - OOOO(IOFF)                                      
        END IF                                                     
      END IF                                   
!                                                                 
!   ..... Part B                                                  
!                                                                 
      IF(KREP.EQ.KSREP) THEN                                      
        IF(K.EQ.KS) THEN                                          
          ASLREP=MULTB(ASREP,LREP,1)                              
          IOFF=IIVO(ASREP,LREP) + (L-1)*NV(ASREP) + AS - 1        
          IOFF=IOFF*NVO(ASLREP) + IVOVO(ASLREP)                   
          IOFF=IOFF + IIVO(AREP,LSREP) +                          &
     &        (LS-1)*NV(AREP) + A                                 
          S = S + VOVO(IOFF)                                      
        END IF                                                     
      END IF                                                       
!                                                                 
!   ..... Part C                                                  
!                                                                 
      IF(LREP.EQ.LSREP) THEN                                      
        IF(L.EQ.LS) THEN                                          
          ASKREP=MULTB(ASREP,KREP,1)                              
          IOFF=IIVO(ASREP,KREP) + (K-1)*NV(ASREP) + AS - 1        
          IOFF=IOFF*NVO(ASKREP) + IVOVO(ASKREP)                   
          IOFF=IOFF + IIVO(AREP,KSREP) +                          &
     &        (KS-1)*NV(AREP) + A                                 
          S = S + VOVO(IOFF)                                      
        END IF                                                     
      END IF                                                       
!                                                                 
!   ..... Part D                                                  
!                                                                 
      IF(LREP.EQ.KSREP) THEN                                      
        IF(L.EQ.KS) THEN                                          
          ASKREP=MULTB(ASREP,KREP,1)                              
          IOFF=IIVO(ASREP,KREP) + (K-1)*NV(ASREP) + AS - 1        
          IOFF=IOFF*NVO(ASKREP) + IVOVO(ASKREP)                   
          IOFF=IOFF + IIVO(AREP,LSREP) +                          &
     &        (LS-1)*NV(AREP) + A                                 
          S = S - VOVO(IOFF)                                      
        END IF                                                     
      END IF               
!                                                                 
!   ..... Part E                                                  
!                                                                 
      IF(KREP.EQ.LSREP) THEN                                      
        IF(K.EQ.LS) THEN                                          
          ASLREP=MULTB(ASREP,LREP,1)                              
          IOFF=IIVO(ASREP,LREP) + (L-1)*NV(ASREP) + AS - 1        
          IOFF=IOFF*NVO(ASLREP) + IVOVO(ASLREP)                   
          IOFF=IOFF + IIVO(AREP,KSREP) +                          &
     &        (KS-1)*NV(AREP) + A                                 
          S = S - VOVO(IOFF)                                      
        END IF                                                     
      END IF                                                       
                                                                  
      fano_kcocc = S                                                   
                                                                  
!      RETURN                                                      
      END FUNCTION     



END MODULE adc_fano_real_routines
