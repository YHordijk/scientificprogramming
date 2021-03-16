C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C      *****************************************************
C      ****        D H F - A D C (2)  M O D U L E       ****
C      ****           TWO-ELECTRON PROPAGATOR           ****
C      ****     D O U B L E   I O N I Z A T I O N S     ****
C      *****************************************************
C
C                       AUTHOR: M.Pernpointner
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ADCDBLE_MAIN (MXCORR,NSP,IREPSPI,EPS)
      use adc_mat
C
C---------------Description-----------------------------------------------
C
C    Setup for the two-electron propagator.
C
C---------------Implicit definition all integers--------------------------
C-------- this is vital, because buffer names are in fact integer
C-------- pointers to the next free space in the huge array CC
C
      IMPLICIT INTEGER(A-Z)
C
C---------------Calling variables--------------------------------------
C
      INTEGER                             ::  MXCORR
      INTEGER                             ::  NSP
      INTEGER                             ::  IREPSPI(NSP,16,2)
      REAL*8                              ::  EPS(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/files.inc"
#include  "adcinpt.inc"
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
      COMMON/RDMATND/NBUFS
C
C
C---------------Local variables memory requirements------------------
C
      LOGICAL ERRFLG
      INTEGER ICOLLEN(MXREP),INTBUF
      REAL(8), PARAMETER  :: autoev = 27.2113957D0
c file names for the Lanczos
      CHARACTER*6 NAME2,NAME3
      PARAMETER(NAME2='ADC2PF',NAME3='ADCCNF')
      CHARACTER*5 NMSPECF
      PARAMETER(NMSPECF='DSPEC')
c variables for the lanczos
      REAL*8 RDUMMY
      INTEGER IONIZLEVEL
      INTEGER IRECL
!
!  declare the buffer arrays needed for the calculations accordingly
!
      integer, allocatable, dimension(:)    :: BFIRCD, BIROOO,
     &         BLKUPO, BLKUPV, BSPTRO, BUFIOI, BUFIOJ

      real*8, allocatable, dimension(:)     :: BAOCOC, BFOOOO, BFVOVO, 
     &    BGIJKL, BIXOOO, BOOOO1, BOOOO2, BORIJK, BSATEV, BUFBUF,
     &    BUFCOL, BUFENE, BVOOO1, BVOOO2, BVVOO1, BVVOO2, E_IJKL,
     &    ERSMIK
C
C---------------Executable code--------------------------------------
C
C
C
C  _________________________ MEMORY ALLOCATION SECTION ______________
C |
C |
C |
      MXCORR=-1
      IONIZLEVEL=2
      WRITE(IW,*) 'Entering 2P propagator module'
C
C determine maximum dimensions for the buffers needed
C
      NDIMT4 = 0
      NOOTMAX = 0
      DO KREP = 1,NREP
        MX1 = NVVOT(KREP)*NO(KREP)*NO(KREP)
        IF(MX1.GT.NDIMT4) NDIMT4 = MX1
        MX2 = NOOT(KREP)
        IF(MX2.GT.NOOTMAX) NOOTMAX = MX2
      ENDDO

      NDIM = IVVOOTT(NREP+1)*RCW
      NDIM = MAX0(NDIM,KVVOOT(NREP+1)*RCW)
      allocate(BVVOO1(NDIM))
      allocate(BVVOO2(NDIM))

      NDIM = IOOOOTT(NREP+1)*RCW
      allocate(BOOOO1(NDIM))
      allocate(BOOOO2(NDIM))

      NDIM = NV1*RCW
      allocate(BGIJKL(NDIM))

      NDIM = IOO(NREP+1)*RCW
      allocate(BAOCOC(NDIM))

      NDIM = NDIMT4
      allocate(ERSMIK(NDIM))

      NDIM = IOOOOTT(NREP+1)
      allocate(E_IJKL(NDIM))
C
      NDIM = MAX0(IVOOOT(NREP+1),LOVOOT(NREP+1)) * RCW
      allocate(BVOOO1(NDIM))
      allocate(BVOOO2(NDIM))
C
      NDIM = IO(NREP+1)*IO(NREP+1)*IO(NREP+1)
      allocate(BIXOOO(NDIM))
      allocate(BORIJK(NDIM))
C
      NDIM = IO(NREP+1)*2
      allocate(BIROOO(NDIM))
C
      M = IO(NREP+1)
      NDIM = (M*M*M*IV(NREP+1) + NOOTMAX)*RCW
      allocate(BUFCOL(NDIM))
      NDIM = (M*M*M*IV(NREP+1) + NOOTMAX)
      allocate(BUFENE(NDIM))
!
! |
! |
! |
! |____________E_N_D______A_L_L_O_C_A_T_I_O_N______S_E_C_T_I_O_N____
!
!
!  _________________________ INTRO SECTION __________________________
! |
! |
! |
      CALL PST('Parameters for the double ionization run:+')

      write(IW,'(10X,A,T50,I3)') 'Order of DIP always ADC-2x !'
      IF(ADCTHR.gt.0.0d0) THEN
        WRITE(IW,*) 'A threshold of',ADCTHR,' will be applied to'
        WRITE(IW,*) 'the 2P ADC matrix.'
      ENDIF
      write(IW,'(10X,A,T50,F18.6)') 'Threshold for screen output:',
     &      dipprnt
      write(IW,'(10X,A)') 'Requested DIP symmetries (dipreps):'
      n=0
      do i=1,32
        n=n+dipreps(i)
      enddo
      if(n.eq.0) then
        write(IW,'(10X,A)') 'None. Program will do it automatically.'
      else
        write(IW,'(10X,16I4)') (dipreps(i),i=1,16)
        write(IW,'(10X,16I4)') (dipreps(i),i=17,32)
      endif
! This output printing became obsolete
!
!      write(IW,'(10X,A)') 'Requested DIP eigenvectors:'
!      n=0
!      do i=1,32
!        n=n+dipeigv(i)
!      enddo
!      if(n.eq.0) then
!        write(IW,'(10X,A)') 'None.'
!      else
!        write(IW,'(10X,16I4)') (dipeigv(i),i=1,16)
!        write(IW,'(10X,16I4)') (dipeigv(i),i=17,32)
!      endif
!
! determine symmetries to be calculated
!
      n = 0
      do i = 1,nrep
        n = n + dipreps(i)
      enddo
      write(iw,'(10X,A)') 'Program calculates the following symmetries:'
      write(iw,*)
      if(n.eq.0) then
        icount=1
        do i=1,nrep
          dipreps(icount) = i
          write(iw,*) '        No: ',icount,'  symmetry: ',i
          icount = icount + 1
        enddo
        write(iw,*)
      else
        icount=1
        do i=1,nrep
          if(dipreps(i).ne.0) then
            write(iw,*) '        No: ',icount,'  symmetry: ',dipreps(i)
            icount = icount + 1
          endif
        enddo
      endif
      write(iw,*)
C |
C |
C |
C |________________E_N_D______S_E_C_T_I_O_N__________________________
C

      INTBUF = 5*1024*1024  !output buffer is now 5 MWORDS (40 MByte)
C
C  __________________  BUILD 2H-2H BLOCK  ___________________________
C |
C |
C |
C   all symmetries in the 2H/2H block are treated simultaneously.
C
      CALL BUILD_2H2H (BGIJKL, BVVOO1, BVVOO2, 
     &                 BOOOO1, BOOOO2, BAOCOC, 
     &                 EPS, ERSMIK, E_IJKL)
C
C   the 2h/2h block is now available in BGIJKL and has
C   OOT||OOT structure with the BKC <|<| |>|>
C |
C |
C |________________E_N_D______S_E_C_T_I_O_N__________________________
C
C
C  __________________  BUILD 3H1P-2H BLOCK  ___________________________
C |
C |   THIS IS TO BE DONE OVER INDIVIDUAL SYMMETRIES.
C |   THE CORRESPONDING ADC-MATRIX BLOCKS ARE STORED PER SYMMETRY.
C |   THE Blocks are of the 2h+3h1p/2h structure. During the
C |   index array setup in 3H1P the satellite diagonal epsilon
C |   array is also calculated and stored in an intermediate file.
C |

c*****************************************************************
c**** master symmetry loop for the coupling/satellite blocks  ****
c**** the 3h1p/2h matrices are written out in this routine    ****
c**** For each symmetry the length of the satellite block     ****
c**** is determined and stored in ICOLLEN (needed for the     ****
c**** SAT block routine. The cc of the coupling block         ****
c**** is implemented. We therefore need to cc the result.     ****
c**** Important: Symmetries where no final states occur       ****
c**** are skipped throughout!                                 ****
c*****************************************************************

      IF(CARITH) THEN
        CALL PST('Building complex 2h/3h1p coupling block+')
      ELSE
        CALL PST('Building real 2h/3h1p coupling block+')
      ENDIF
      DO 200 IRPLOOP = 1,NREP
        IF(dipreps(IRPLOOP).EQ.0) GOTO 200
        IF(dipreps(IRPLOOP).GT.NREP) THEN
           WRITE(IW,*) 'Chosen Symmetry',dipreps(IRPLOOP),
     &     ' does not exist and is skipped!'
           GOTO 200
        ENDIF

        DESREP = dipreps(IRPLOOP)

        ICOLLEN(DESREP) = 0
        IF(NOOT(DESREP).GT.0) THEN
          WRITE(IW,'(A,I4)') 
     &       ' *** Found 2H final states in symm. ',DESREP
          IF(CARITH) THEN
            CALL BUILD_3H1P_2H_C(IO(NREP+1),BIXOOO,BIROOO,
     &                           BORIJK,BVOOO1,BVOOO2,
     &                           BUFCOL,BUFENE,BGIJKL,
     &                           EPS,DESREP,LENSATBL)
          ELSE
            CALL BUILD_3H1P_2H_R(IO(NREP+1),BIXOOO,BIROOO,
     &                           BORIJK,BVOOO1,BVOOO2,
     &                           BUFCOL,BGIJKL,EPS,DESREP,
     &                           LENSATBL)
          ENDIF
          ICOLLEN(DESREP) = LENSATBL
        ENDIF
 200  CONTINUE

c
c  releasing dynamically allocated memory used so far
c
      deallocate(BVVOO1)
      deallocate(BVVOO2)
      deallocate(BOOOO1)
      deallocate(BOOOO2)
      deallocate(BGIJKL)
      deallocate(BAOCOC)
      deallocate(ERSMIK)
      deallocate(E_IJKL)
      deallocate(BVOOO1)
      deallocate(BVOOO2)
      deallocate(BIXOOO)
      deallocate(BORIJK)
      deallocate(BIROOO)
      deallocate(BUFCOL)
      deallocate(BUFENE)
C |
C |
C |
C |________________E_N_D______S_E_C_T_I_O_N__________________________
C
C  __________________  BUILD 3H1P-3H1P (SAT) BLOCK  _________________
C |
C |   This is done over individual symmetries. The packed 3h1p/2h
C |   matrices are already on disk. For each symmetry
C |   the required memory is determined and freed again after
C |   termination.

      IF(CARITH) THEN
        CALL PST('Building 1. order complex satellite block+')
      ELSE
        CALL PST('Building 1. order real satellite block+')
      ENDIF

c
c  allocate space 
c  - for OOOO and VOVO integrals (needed in all symmetries)
c  - for the spinor rel/abs transposition arrays
c  - and for the master 4-index lookup maps OOOO and VOVO
c  ATT:  the array for the OOOO and VOVO integrals is now increased
c        by one entry in order to hold the OOOO(0)=0.0d0 value.
C this makes offset addressing particularly simple but may cause a
C slight computational overhead by adding zero to a floating number.
c
      NDIM=(IOOOOTT(NREP+1)+1)*RCW
      allocate(BFOOOO(NDIM))
      NDIM=(IVOVO(NREP+1)+1)*RCW
      allocate(BFVOVO(NDIM))
      NDIM=IO(NREP+1)*2
      allocate(BSPTRO(NDIM))
      NDIM=IO(NREP+1)
      NDIM=NDIM*NDIM*NDIM*NDIM
      allocate(BLKUPO(NDIM))
      NDIM=IV(NREP+1)*IO(NREP+1)
      NDIM=NDIM*NDIM
      allocate(BLKUPV(NDIM))
c
c  initialize transposition arrays
c
      CALL MAKESPTRA(BSPTRO,IO(NREP+1))
c
c  initialize lookup tables for OOOO and VOVO
c
      CALL MAKELOOKUPO(BLKUPO,IO(NREP+1))
      CALL MAKELOOKUPV(BLKUPV,IO(NREP+1),IV(NREP+1))
c
c  grab buffer space for sparse matrix storage
c
      allocate(BUFIOI(INTBUF))
      allocate(BUFIOJ(INTBUF))
      allocate(BUFBUF(INTBUF*RCW))
!
!
!___________________________________________________________
!| loop over the symmetries
!| att: the additional, symmetry-specific memory is allocated
!| within this loop and immediately released afterwards.
!| att2: it was observed that the ADC2 matrices get very large and
!| therefore a diagonalization is required immediately after the
!| SAT block construction!
!|     
!
      DO 300 IRPLOOP = 1,NREP
        IF(dipreps(IRPLOOP).EQ.0) GOTO 300
        IF(dipreps(IRPLOOP).GT.NREP) THEN
           WRITE(IW,*) 'Chosen Symmetry',dipreps(IRPLOOP),
     &     ' does not exist and is skipped!'
           GOTO 300
        ENDIF

        DESREP = dipreps(IRPLOOP)

        IF(NOOT(DESREP).EQ.0) GOTO 300
        WRITE(IW,*) ' *** Found final states in symm.',DESREP

        LENSATBL = ICOLLEN(DESREP)
        NDIM = LENSATBL
        allocate(BSATEV(NDIM))
        NDIM = LENSATBL*4
        allocate(BFIRCD(NDIM))

        IF(CARITH) THEN
          CALL BUILD_SAT_C(BFOOOO,BFVOVO,BSPTRO,IO(NREP+1),
     &                     IV(NREP+1),BSATEV,BFIRCD,LENSATBL,
     &                     BLKUPO,BLKUPV,DESREP,ADCTHR,
     &                     BUFIOI,BUFIOJ,BUFBUF,INTBUF)
        ELSE
          CALL BUILD_SAT_R(BFOOOO,BFVOVO,BSPTRO,IO(NREP+1),
     &                     IV(NREP+1),BSATEV,BFIRCD,LENSATBL,
     &                     BLKUPO,BLKUPV,DESREP,ADCTHR,
     &                     BUFIOI,BUFIOJ,BUFBUF,INTBUF)
        ENDIF
C
C  call 2h/3h1p configuration creator, file handle is available
C  afterwards
C
        LADC = NOOT(DESREP) + LENSATBL
        CALL WCONDAT_D(ITAPADC,NAME3,IRECL,DESREP,LADC)
        WRITE(IW,*) 'Record length as determined by WCONDAT:',IRECL
c
c Diagonalizers are decoupled and reallocate memory individually.
c
c
C  set matrix descriptor and perform Lanczos diagonalization
C
        reladc_md_iobase       = ITAPADC
        reladc_md_ionizl       = IONIZLEVEL
        reladc_md_ioldnew      = 2
        reladc_md_intbuf       = INTBUF
        reladc_md_desrep       = DESREP
        reladc_md_rcw          = RCW
        reladc_md_lnzitr       = DIPITER
        reladc_md_matdim       = LADC
        reladc_md_irecl        = IRECL
        reladc_md_nmain        = NOOT(DESREP)
        reladc_md_nbufs        = NBUFS
        reladc_md_eeigv_lower  = reladc_md_dip_eeigv(IRPLOOP,1)/autoev
        reladc_md_eeigv_upper  = reladc_md_dip_eeigv(IRPLOOP,2)/autoev
        reladc_md_fileadc      = NAME2
        reladc_md_filediag     = '      '
        reladc_md_filecnf      = NAME3
        reladc_md_nmspec       = NMSPECF


        DIPPRNT = 30.0d0
c
c  call full diagonalizer for new matrix format if requested.
c
        IF(DOFULL) THEN
          IF(CARITH) THEN
            CALL FULLDIAC2(ITAPADC, NAME2, INTBUF, NBUFS, LADC, DESREP)
          ELSE
            CALL FULLDIAR2(ITAPADC, NAME2, INTBUF, NBUFS, LADC, DESREP)
          ENDIF
        ENDIF
c
c  call Lanczos diagonalizer.
c  memory handling of Lanczos is self-contained.
c
        CALL DIAG_LANC(IW,DOINCORE)
c
c  we do Mulliken-type ADC population analysis after diagonalization
c  since in this part of the program we have access to all the symmetry
c  information we need. Since we need only the 2h part we can make use
c  of the TMATVEC matrix, the eigenvectors of the Lanczos matrix.
c
c  in the corresponding symmetry DESREP
c
c  ATT in this F90 routine an include symm.inc does not work!
c  therefore we have to transfer the relevant arrays to TWOHPOP
c  especially MULTB(64,64,2), NOOT(*) and NO(*)
c
        IF(DOADCPOP) THEN
          IF(CARITH) THEN
            CALL TWOHPOP_C(IW,ITAPADC,DESREP,NREP,MULTB,NO,NOOT)
          ELSE
            CALL TWOHPOP_R(IW,ITAPADC,DESREP,NREP,MULTB,NO,NOOT)
          ENDIF
        ENDIF
! 
!   deallocate symmetry-dependent buffers reserved in each symmetry cycle
!
        deallocate(BSATEV)
        deallocate(BFIRCD)
!|
!|
 300  CONTINUE
!|
!|                        symmetry loop end
!|                        deallocate remaining buffers
!|_____________________________________________________________

      deallocate(BFOOOO)
      deallocate(BFVOVO)
      deallocate(BSPTRO)
      deallocate(BLKUPO)
      deallocate(BLKUPV)
      deallocate(BUFIOI)
      deallocate(BUFIOJ)
      deallocate(BUFBUF)
c
c   perform a final DNM matrix check before we transfer back control
c
      IF(DOADCPOP) THEN
        call DNMMATCHK(IW,ITAPADC,ERRFLG)
        IF(ERRFLG) then
          call quit('dnm matrix check failed!')
        ELSE
          write(iw,*) '--------  DNM matrix check passed.  -------'
        ENDIF
      ENDIF

C |
C |
C |
C |________________E_N_D______S_E_C_T_I_O_N__________________________
c
      RETURN
      END
C |
C |
C |
C |________________E_N_D___M A I N___P R O G R A M __________________

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_2H2H (GIJKL,VVOO1,VVOO2,OOOO1,OOOO2,
     &                       AIK,EPS,EABLKKS,EIJKL)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 2H/2H block of the Double ionization ADC Matrix
C     The result goes into GIJKL
C
C---------------Calling variables--------------------------------------
C
      REAL*8 GIJKL(*), VVOO1(*), VVOO2(*), OOOO1(*), OOOO2(*)
      REAL*8 AIK(*), EPS(*), EABLKKS(*), EIJKL(*)
      REAL*8 ASYM
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
      CALL PST('Building 2h/2h block+')
c
c  -------------------------------
c  Start calculating  TERM 5A    |
c  RESULT in OOOO1               |
c  -------------------------------
c
c  Form the v_rs,kl (VVOO2) and v*_rs,ij (VVOO1)
c  then multiply v_rs,kl by -(e_r + e_s) to give
c  Y_rs,kl (VVOO2). Att: the factor 0.5 is absorbed because
c  the sum over r,s is replaced by the sum r>s !
c  then contract v*_rs,ij and Y_rs,kl to C_ij,kl (OOOO1)
c
      NUM = IOOOOTT(NREP+1)
      CALL XCOPY(NUM,A0,0,OOOO1,1)
      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL XCOPY(NV3,VVOO1,1,VVOO2,1)
      IF(CARITH) CALL CONJUGA(NV3,VVOO1,1)
      CALL DENOMY (EPS,VVOO2,VVOO2)
      CALL CNTRCT('T','N',NOOT,NOOT,NVVT,A1,VVOO1,VVOO2,A0,OOOO1,NREP)
c
c  -------------------------------
c  Start calculating  TERM 5B    |
c  RESULT in OOOO2               |
c  -------------------------------
c
c  energy array 0.5*(ei + ej + ek + el) in EIJKL (always real)
c  now v_rs,kl is in VVOO1 and v*_rs,ij is in VVOO2
c  form B_ij,kl  BKC: <|<| |>|>  being in OOOO2
c  then EIJKL has to be multiplied *elementwise* with B_ij,kl
c  not contracted ! result is final B_ij,kl (in OOOO2)
c
      NUM = IOOOOTT(NREP+1)
      CALL XCOPY(NUM,A0,0,OOOO2,1)
      CALL MAKE_EIJKL(EPS,EIJKL)
      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL XCOPY(NV3,VVOO1,1,VVOO2,1)
      IF(CARITH) CALL CONJUGA(NV3,VVOO2,1)
      CALL CNTRCT('T','N',NOOT,NOOT,NVVT,A1,VVOO2,VVOO1,A0,OOOO2,NREP)
      NUM=IOOOOTT(NREP+1)
      IF(CARITH) THEN
         CALL EWARRMULC(EIJKL,OOOO2,NUM)
      ELSE
         CALL EWARRMULR(EIJKL,OOOO2,NUM)
      ENDIF
c
c  -------------------------------
c  Start calculating  AIK block  |
c  RESULT in AIK(*)              |
c  -------------------------------
c
c compute the AIK block (needed for term 1-4) for all irreps.
c This is exactly the same as the ADC3 C(2)_kk' block !
c But we include it as a separate routine because MAKE_HH2 does
c additional things which to cancel would be too much fiddling around.
c Att. AIK has to be cleared because MAKE_HH2 *ADDS* on it !!!!!
c at the same time determine dimension of gamma (as in the main program)
c
      DIMGM = 0
      AOFFS = 1
      DO KREP = 1,NREP
        NUM = NO(KREP) * NO(KREP)
        CALL XCOPY(NUM,A0,0,AIK(AOFFS),1)
        CALL MAKE_AIK(AIK(AOFFS),VVOO1,VVOO2,EPS,EABLKKS,KREP)
        CALL CKKASYM(AIK(AOFFS),AIK(AOFFS),NO(KREP),ASYM)

        AOFFS = AOFFS + NUM*RCW
        DIMGM = DIMGM + NOO(KREP)*NOO(KREP)
      ENDDO
c
c  before processing the A_ik and diagonal (Term 6) expressions we
c  collect terms 0/5a/ and 5b. All have  OOT/OOT structure!
C  att: the GIJKL array is a bit larger since it has to keep the full
C  matrix!
c
c clear Gamma_ijkl array (2h/2h block)
c
      CALL XCOPY (NV1,A0,0,GIJKL,1)
c
c Term 0: fill in the  <OO||OO> integrals
c remember: the OOT/OOT space is smaller than the GIJKL space.
c
      CALL GETOOOO(GIJKL)
c
c add term 5A and term 5B to term 0
c
      CALL XAXPY (NV1,A1,OOOO1,1,GIJKL,1)
      CALL XAXPY (NV1,A1,OOOO2,1,GIJKL,1)
c
c Term 1-4 and Term 6, done by AWEAVER
c Term 6 is the mere diagonal element and contains only spinor energies
c Remember: the antisymmetry with respect to the A_ik contributions is
c accounted for in AWEAVER!
c
      IF(CARITH) THEN
        CALL AWEAVERC(EPS,GIJKL,AIK)
      ELSE
        CALL AWEAVERR(EPS,GIJKL,AIK)
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_3H1P_2H_R (MAXOCCSP,IXOOO,IROOO,ORIJK,
     &                            VOOO1,VOOO2,COL,GIJKL,EPS,DESREP,
     &                            LEN_RKLM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/2H block of the Double ionization ADC Matrix
C     This is the real number implementation.
c
C---------------Calling variables--------------------------------------
C
      INTEGER MAXOCCSP
      INTEGER IXOOO(MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER IROOO(MAXOCCSP,2)
      INTEGER ORIJK(MAXOCCSP*MAXOCCSP*MAXOCCSP)
      REAL*8 VOOO1(*),VOOO2(*),COL(*),GIJKL(*),EPS(*)
      INTEGER DESREP,LEN_RKLM
c
c  MAXOCCSP is the total number of occupied spinors over all ireps.
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUX
      CHARACTER*6 NAME3,NAME4
      PARAMETER(NAME3='ADC2PE',NAME4='ADC2PP')
      CHARACTER*9 FNAME3,FNAME4
      REAL*8 ENERGY_I,ENERGY_J,ENERGY_K,ENERGY_R
C
C---------------Executable code--------------------------------------
c
c  desrep is the destination symmetry of the two-hole final state
c  the packed matrix for the coupling block just consists of a
c  sequential packed file which will be reorganized in the routine
c  for the satellite block and transformed to its final form.
c
c_______  INITIALIZE ARRAYS FOR THE COUPL. BLOCK ____________________
c                                                                    |
c                                                                    |
c
c  all occupied spinors are counted contiguously through the
c  irreps. Assign the corresponding irep to this spinor number.
c  in the code the spinors are accessed through this unique number.
c  this unique number can be interpreted as *absolute* spinor number
c  in contrast to the *relative* spinor number per irep.
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
c
c  set up triple index array for totally antisymmetric storage of 
c  the Gamma_rklm,ij
c  The gamma array is totally antisymmetric in klm!
c  This means that I>J>K must be maintained throughout and
c  no I=J, J=K or I=K combination must occur 
c
      IX=1
      DO K=1,MAXOCCSP
        DO J=1,MAXOCCSP
          DO I=1,MAXOCCSP
            IXOOO(I,J,K) = 0
            ORIJK(IX)=0
            IX = IX + 1
          ENDDO
        ENDDO
      ENDDO
c
c  initialize column offset array.
c  In due course also calculate the energies and use the COL(*) 
c  buffer for storage. Later on the COL(*) buffer holds the PHHH 
c  column.
c  For the chosen storage mode I,J,K with I>J>K we additionally
c  know the corresponding IREPS and can then access the one-particle
c  spinor energies accordingly.
c
      IX = 1
      IEX = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            IXOOO(I,J,K)=IX
            IX = IX + 1
            RRP = MULTB(DESREP+NREP,IJKRP,2)
c
c  we have determined G_r such that G_r* x G_ijk = G_dest
c  attention: DESREP is bosonic !
c
            ORIJK(IX) = ORIJK(IX-1) + NV(RRP)
c
c  energies. from here I,J,K and IRP,JRP,KRP,RRP are known.
c  We loop over R for the different particle indices and construct the
c  energy array.
c
            ENERGY_I = EPS(IO(IRP) + IROOO(I,1))
            ENERGY_J = EPS(IO(JRP) + IROOO(J,1))
            ENERGY_K = EPS(IO(KRP) + IROOO(K,1))
            AUX = ENERGY_I + ENERGY_J + ENERGY_K
            DO 17 IXR = 1,NV(RRP)
              ENERGY_R = EPS(IO(NREP+1) + IV(RRP) + IXR)
              COL(IEX) = ENERGY_R - AUX
              IEX = IEX + 1
 17         CONTINUE

          ENDDO
        ENDDO
      ENDDO
      LEN_RKLM=ORIJK(IX)   
      WRITE(IW,'(4X,A,I3,A,I8)') 'Length of coupling block in symmetry',
     &      DESREP,':',LEN_RKLM
c                                                                    |
c____________________________________________________________________|
c_____________________________________________________________________
c                                                                    |
c            OPEN CONTROL/MATRIX/ENERGY FILES FOR THE ADC-2          |
c                                                                    |
      IF(DESREP.GT.9) THEN
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
        WRITE(FNAME4,'(A6,A1,I2)') NAME4,'.',DESREP
      ELSE
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
        WRITE(FNAME4,'(A6,A2,I1)') NAME4,'.0',DESREP
      ENDIF

c     WRITE(IW,*) '2P-ADC block data (packed) in       ',FNAME4
c     WRITE(IW,*) '2P-ADC matrix satellite energies in ',FNAME3

      OPEN(ITAPADC+3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAPADC+4,FILE=FNAME4,FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c       write satellite block diagonal energies immediately
c
      WRITE(ITAPADC+3) (COL(IX),IX=1,LEN_RKLM)
      CLOSE(ITAPADC+3)
c                                                                    |
c____________________________________________________________________|
c
c  read VOOO integrals and sort to VOO,O:  (I K>L, J)
c
      CALL GETVOOO(VOOO1)
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           VOOO1,VOOO2)
c
C *****************************************************
C ************  MASTER SYMMETRY LOOP OVER THE 2H SYMS *
C *****************************************************
C
c
c  the externally given symmetry is DESREP. Of course, there are
c  multiple combinations of gamma_j and gamma_i to fulfill
c  gamma_i x gamma_j = DESREP. They all have to be included and
c  form the final array in the genuine order below.
c  this genuine order also defines the order of columns in GIJKL
c
      GCOL=0
      DO 90 JXREP = 1, NREP
        IXREP = MULTB(JXREP,DESREP+NREP,2)
        IF (IXREP.LT.JXREP) GOTO 90
        DO 85 JX = 1, NO(JXREP)
          IXMIN = 1
          IF (IXREP.EQ.JXREP) IXMIN = JX + 1
          DO 80 IX = IXMIN, NO(IXREP)

c
c  the following initializations and updates *HAVE* to take place
c  for *EACH* new column, determined by Gamma_i, Gamma_j, IX and JX:
c
c  - GCOL:  increase column counter by one.
c  - CLEAR intermediate column COL(*). Contributions are added on top
c          therefore it must be set to zero for each new column. Once
c          a col is finished it is written out immediately.
c
            GCOL = GCOL + 1
            DO IY=1,LEN_RKLM
              COL(IY) = 0.0D0 
            ENDDO
c
c  from here we know the ireps IREP,JREP and the corresponding spinor
c  indices i and j. Now we rumble through the sorted VOOOs and distribute the
c  integrals
c_______________________________________________________________
c|
c|      internal VOOO Loop start
c|  
c|  first treat the delta_i? cases
c|
      OFFOVOO = LOVOOT(JXREP) + (JX-1)*NVOOT(JXREP)
      ICNT = 1
      DO 55 LMREP=1,NREP
        DO 50 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 50
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(JXREP,LMREP+NREP,2)

c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset.
c
c  *** Attention:  IX, JX, L and M are relative spinor numbers
c  but the IXOOO array needs absolute ones ! 
c  the R spinor index stays relative...
c
c  the rklm values are created by looping through the VOO,O
c  integrals only once. Successively the loop variables
c  L and M are interpreted as L,M/K,M and K,L !
c  
               DO 45 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 AUX = VOOO2(OFFS)
 
                 XXK = IO(IXREP) + IX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - AUX
                 ENDIF
                   
                 XXK = IO(LREP)  + L
                 XXL = IO(IXREP) + IX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + AUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(IXREP) + IX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - AUX
                 ENDIF

 45            CONTINUE
   
               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 50     CONTINUE
 55   CONTINUE
c|
c|  second treat the delta_j? case
c|
      OFFOVOO = LOVOOT(IXREP) + (IX-1)*NVOOT(IXREP)
      ICNT = 1
      DO 75 LMREP=1,NREP
        DO 70 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 70
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(IXREP,LMREP+NREP,2)
c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset
c
               DO 65 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 AUX = VOOO2(OFFS)

                 XXK = IO(JXREP) + JX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + AUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(JXREP) + JX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - AUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(JXREP) + JX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + AUX
                 ENDIF

 65            CONTINUE

               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 70     CONTINUE
 75   CONTINUE
c|
c|      internal VOOO Loop end
c|  
c|______________________________________________________________
c
c  from here one column of the coupling block matrix is finished
c  and ranges from 1 to LEN_RKLM
c  Next the corresponding 2h column
c  will be placed in front of it and then written to disk.
c_______________________________________________________________
c   FORM COMPLETE COLUMN AND WRITE TO DISK IN PACKED FORM.     |
c  The packed matrix storage implies column counting starting  |
c  from 1 and the next column is indicated by a row index of   |
c  -1. EOF is indicated by a row index of -255.                |
c                                                              |
          JOFF = NOOT(DESREP)
          GOFFSET = IOOOOTT(DESREP) + (GCOL-1)*JOFF
          DO JXX = LEN_RKLM,1,-1
            COL(JXX+JOFF) = COL(JXX)
          ENDDO
          DO JXX = 1,JOFF
            COL(JXX) = GIJKL(GOFFSET + JXX)
          ENDDO
          JXXM = LEN_RKLM + JOFF
          DO JXX = GCOL, JXXM
            IF(COL(JXX).NE.0.0D0) THEN
              WRITE(ITAPADC+4) JXX,COL(JXX)
            ENDIF
          ENDDO
          JXX = -1
          WRITE(ITAPADC+4) JXX,0.0d0
c                                                              |
c______________________________________________________________|
c
c
  80      CONTINUE
  85    CONTINUE
  90  CONTINUE
c
c  close adc matrix file (packed and unpacked)
c  write end of file descriptor.
c
      JXX = -255
      WRITE(ITAPADC+4) JXX,0.0d0
      CLOSE(ITAPADC+4)

      RETURN
      END
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_3H1P_2H_C(MAXOCCSP,IXOOO,IROOO,ORIJK,
     &                           VOOO1,VOOO2,COL,ECOL,GIJKL,
     &                           EPS,DESREP,LEN_RKLM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/2H block of the Double ionization ADC Matrix
C     This is the complex number implementation.
C
C---------------Calling variables--------------------------------------
C
      INTEGER MAXOCCSP
      INTEGER IXOOO(MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER IROOO(MAXOCCSP,2)
      INTEGER ORIJK(MAXOCCSP*MAXOCCSP*MAXOCCSP)
      COMPLEX*16 VOOO1(*),VOOO2(*),COL(*),GIJKL(*)
      REAL*8 ECOL(*),EPS(*)
      INTEGER DESREP,LEN_RKLM
c
c  MAXOCCSP is the total number of occupied spinors over all ireps.
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 AUX
      COMPLEX*16 CAUX
      CHARACTER*6 NAME3,NAME4
      PARAMETER(NAME3='ADC2PE',NAME4='ADC2PP')
      CHARACTER*9 FNAME3,FNAME4
      REAL*8 ENERGY_I,ENERGY_J,ENERGY_K,ENERGY_R
C
C---------------Executable code--------------------------------------
c
c  desrep is the destination symmetry of the two-hole final state
c  the packed matrix for the coupling block just consists of a
c  sequential packed file which will be reorganized in the routine
c  for the satellite block and transformed to its final form.
c
c_______  INITIALIZE ARRAYS FOR THE COUPL. BLOCK ____________________
c                                                                    |
c                                                                    |
c
c  all occupied spinors are counted contiguously through the
c  irreps. Assign the corresponding irep to this spinor number.
c  in the code the spinors are accessed through this unique number.
c  this unique number can be interpreted as *absolute* spinor number
c  in contrast to the *relative* spinor number.
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
c
c  set up triple index array for totally antisymmetric storage of 
c  the Gamma_rklm,ij
c  The gamma array is totally antisymmetric in klm!
c  This means that I>J>K must be maintained throughout!
c  no I=J, J=K or I=K combination must occur 
c
      IX=1
      DO K=1,MAXOCCSP
        DO J=1,MAXOCCSP
          DO I=1,MAXOCCSP
            IXOOO(I,J,K) = 0
            ORIJK(IX)=0
            IX = IX + 1
          ENDDO
        ENDDO
      ENDDO
c
c  initialize column offset array.
c  In due course also calculate the energies and use the COL(*) 
c  buffer for storage. Later on the COL(*) buffer holds the PHHH 
c  column.
c  For the chosen storage mode I,J,K with I>J>K we additionally
c  know the corresponding IREPS and can then access the one-particle
c  spinor energies accordingly.
c
      IX = 1
      IEX = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            IXOOO(I,J,K)=IX
            IX = IX + 1
            RRP = MULTB(DESREP+NREP,IJKRP,2)
c
c  we have determined G_r such that G_r* x G_ijk = G_dest
c  attention: DESREP is bosonic !
c
            ORIJK(IX) = ORIJK(IX-1) + NV(RRP)
c
c  energies. from here I,J,K and IRP,JRP,KRP,RRP are known.
c  We loop over R for the different particle indices and construct the
c  energy array.
c
            ENERGY_I = EPS(IO(IRP) + IROOO(I,1))
            ENERGY_J = EPS(IO(JRP) + IROOO(J,1))
            ENERGY_K = EPS(IO(KRP) + IROOO(K,1))
            AUX = ENERGY_I + ENERGY_J + ENERGY_K
            DO 17 IXR = 1,NV(RRP)
              ENERGY_R = EPS(IO(NREP+1) + IV(RRP) + IXR)
              ECOL(IEX) = ENERGY_R - AUX
              IEX = IEX + 1
 17         CONTINUE

          ENDDO
        ENDDO
      ENDDO
      LEN_RKLM=ORIJK(IX)   
      WRITE(IW,'(4X,A,I3,A,I8)') 'Length of coupling block in symmetry',
     &      DESREP,':',LEN_RKLM
c                                                                    |
c____________________________________________________________________|
c_____________________________________________________________________
c                                                                    |
c            OPEN CONTROL/MATRIX/ENERGY FILES FOR THE ADC-2          |
c                                                                    |
      IF(DESREP.GT.9) THEN
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
        WRITE(FNAME4,'(A6,A1,I2)') NAME4,'.',DESREP
      ELSE
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
        WRITE(FNAME4,'(A6,A2,I1)') NAME4,'.0',DESREP
      ENDIF

      OPEN(ITAPADC+3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAPADC+4,FILE=FNAME4,FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c write satellite block diagonal energies immediately
c
      WRITE(ITAPADC+3) (ECOL(IX),IX=1,LEN_RKLM)
      CLOSE(ITAPADC+3)
c                                                                    |
c____________________________________________________________________|
c
c  read VOOO integrals and sort to VOO,O:  (I K>L, J)
c
      CALL GETVOOO(VOOO1)
      CALL SRT6 (NREP,MULTB,LFA,NVO,NV,NO,NOOT,NVOOT,LOVOOT,LLVOOT,
     &           VOOO1,VOOO2)
C
C *****************************************************
C ************  MASTER SYMMETRY LOOP OVER THE 2H SYMS *
C *****************************************************
C
c
c  the externally given symmetry is DESREP. Of course, there are
c  multiple combinations of gamma_j and gamma_i to fulfill
c  gamma_i x gamma_j = DESREP. They all have to be included and
c  form the final array in the genuine order below.
c  this genuine order also defines the order of columns in GIJKL
c
      GCOL=0
      DO 90 JXREP = 1, NREP
        IXREP = MULTB(JXREP,DESREP+NREP,2)
        IF (IXREP.LT.JXREP) GOTO 90
        DO 85 JX = 1, NO(JXREP)
          IXMIN = 1
          IF (IXREP.EQ.JXREP) IXMIN = JX + 1
          DO 80 IX = IXMIN, NO(IXREP)

c
c  the following initializations and updates *HAVE* to take place
c  for *EACH* new column, determined by Gamma_i, Gamma_j, IX and JX:
c
c  - GCOL:  increase column counter by one.
c  - CLEAR intermediate column COL(*). Contributions are added on top
c          therefore it must be set to zero for each new column. Once
c          a col is finished it is written out immediately.
c
            GCOL = GCOL + 1
            DO IY=1,LEN_RKLM
              COL(IY) = (0.0D0,0.0D0) 
            ENDDO
c
c  from here we know the ireps IREP,JREP and the corresponding spinor
c  indices i and j. Now we rumble through the sorted VOOOs and distribute the
c  integrals
c_______________________________________________________________
c|
c|      internal VOOO Loop start
c|  
c|  first treat the delta_i? cases
c|
      OFFOVOO = LOVOOT(JXREP) + (JX-1)*NVOOT(JXREP)
      ICNT = 1
      DO 55 LMREP=1,NREP
        DO 50 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 50
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(JXREP,LMREP+NREP,2)

c       write(iw,'(A,4I3)') 'In LM LOOP:',LREP,MREP,L,M
c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset.
c
c  *** Attention:  IX, JX, L and M are relative spinor numbers
c  but the IXOOO array needs absolute ones ! 
c  the R spinor index stays relative...
c
c  the rklm values are created by looping through the VOO,O
c  integrals only once. Successively the loop variables
c  L and M are interpreted as L,M/K,M and K,L !
c  
               DO 45 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 CAUX = VOOO2(OFFS)
 
                 XXK = IO(IXREP) + IX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - CAUX
                 ENDIF
                   
                 XXK = IO(LREP)  + L
                 XXL = IO(IXREP) + IX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + CAUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(IXREP) + IX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - CAUX
                 ENDIF

 45            CONTINUE
   
               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 50     CONTINUE
 55   CONTINUE
c|
c|  second treat the delta_j? case
c|
      OFFOVOO = LOVOOT(IXREP) + (IX-1)*NVOOT(IXREP)
      ICNT = 1
      DO 75 LMREP=1,NREP
        DO 70 MREP=1,NREP
          LREP = MULTB(MREP,LMREP+NREP,2)
           IF (LREP.LT.MREP) GOTO 70
           DO M=1,NO(MREP)
             LMIN = 1
             IF (LREP.EQ.MREP) LMIN = M + 1
             DO L = LMIN, NO(LREP)
               RREP=MULTB(IXREP,LMREP+NREP,2)
c
c  from here all the rlm Ireps and rlm indices are known
c  we can assign a destination offset
c
               DO 65 R=1,NV(RREP)

                 OFFS = (OFFOVOO + (ICNT-1) + R)
                 CAUX = VOOO2(OFFS)

                 XXK = IO(JXREP) + JX
                 XXL = IO(LREP)  + L
                 XXM = IO(MREP)  + M
                 XX1 = IXOOO(XXK,XXL,XXM)
                 IF(XX1.NE.0) THEN
                   OFFD = (ORIJK(XX1) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + CAUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(JXREP) + JX
                 XXM = IO(MREP)  + M
                 XX2 = IXOOO(XXK,XXL,XXM)
                 IF(XX2.NE.0) THEN
                   OFFD = (ORIJK(XX2) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) - CAUX
                 ENDIF

                 XXK = IO(LREP)  + L
                 XXL = IO(MREP)  + M
                 XXM = IO(JXREP) + JX
                 XX3 = IXOOO(XXK,XXL,XXM)
                 IF(XX3.NE.0) THEN
                   OFFD = (ORIJK(XX3) + R)
          if(offd.gt.len_rklm) stop 'offd problem'
                   COL(OFFD) = COL(OFFD) + CAUX
                 ENDIF

 65            CONTINUE

               ICNT = ICNT + NV(RREP)

             ENDDO
           ENDDO
 70     CONTINUE
 75   CONTINUE
c|
c|      internal VOOO Loop end
c|  
c|______________________________________________________________
c
c  from here one column of the coupling block matrix is finished
c  and ranges from 1 to LEN_RKLM. Since the c.c. was created according
c  to the formulation of the equations we have
c  to complex conjugate the corresponding vector in the coupling block.
c  the remaining 2h/2h part must not be touched !

          DO JXX = 1,LEN_RKLM
            COL(JXX) = DCONJG(COL(JXX))
          ENDDO

c  Next the corresponding 2h column
c  will be placed in front of it and then written to disk.
c_______________________________________________________________
c   FORM COMPLETE COLUMN AND WRITE TO DISK IN PACKED FORM.     |
c  The packed matrix storage implies column counting starting  |
c  from 1 and the next column is indicated by a row index of   |
c  -1. EOF is indicated by a row index of -255.                |
c                                                              |
          JOFF = NOOT(DESREP)
          GOFFSET = IOOOOTT(DESREP) + (GCOL-1)*JOFF
          DO JXX = LEN_RKLM,1,-1
            COL(JXX+JOFF) = COL(JXX)
          ENDDO
          DO JXX = 1,JOFF
            COL(JXX) = GIJKL(GOFFSET + JXX)
          ENDDO
          JXXM = LEN_RKLM + JOFF
          DO JXX = GCOL, JXXM
            IF(COL(JXX).NE.0.0D0) THEN
              WRITE(ITAPADC+4) JXX,COL(JXX)
            ENDIF
          ENDDO
          JXX = -1
          WRITE(ITAPADC+4) JXX,(0.0d0,0.0d0)
c                                                              |
c______________________________________________________________|
c
c
  80      CONTINUE
  85    CONTINUE
  90  CONTINUE
c
c  close adc matrix file (packed and unpacked)
c  write end of file descriptor.
c
      JXX = -255
      WRITE(ITAPADC+4) JXX,(0.0d0,0.0d0)
      CLOSE(ITAPADC+4)

      RETURN
      END
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_SAT_R(OOOO,VOVO,IROOO,MAXOCCSP,MAXVIRSP,
     &                       SATEV,IRCD,IRCDLEN,LTO,LTV,DESREP,
     &                       ADCTHR,IOI,IOJ,BUF,INTBUF)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/3H1P block of the Double ionization ADC Matrix
C     real number implementation. Attention: row/column matrix storage
C     happens via packed integers ! One can safely assume that the
C     matrix dimension of the ADC matrix never will exceed 2^31-1 !!
C
C---------------Calling variables--------------------------------------
C
      REAL*8 OOOO(0:*),VOVO(0:*)
      INTEGER IROOO(MAXOCCSP,2),MAXOCCSP,MAXVIRSP
      INTEGER IRCD(4,IRCDLEN),IRCDLEN,DESREP
      INTEGER LTO(MAXOCCSP,MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER LTV(MAXVIRSP,MAXOCCSP,MAXVIRSP,MAXOCCSP)
      REAL*8 SATEV(IRCDLEN)
      REAL*8 ADCTHR
c for the packed matrix output file
      INTEGER INTBUF
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      REAL*8 BUF(INTBUF)
c
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"

      COMMON/RDMATND/NBUFS
C
C---------------Local variables--------------------------------------
C
      REAL*8 G
      INTEGER IB,NBUFS,JDUMMY,NTOTELE
      CHARACTER*6 NAME2,NAME3,NAME4
      PARAMETER(NAME2='ADC2PE',NAME3='ADC2PP',NAME4='ADC2PF')
      CHARACTER*9 FNAME2,FNAME3,FNAME4
      INTEGER R,I,J,K,S,L,M,N
      CHARACTER*44 FIELD

      REAL*8 SPARS,XNZEL,XTOTEL
c
c  meaning of the file names:
c  ADC2PC:  control file with the dimensions and file infos
c  ADC2PE:  contains the 3h1p/3h1p diagonal energies
c  ADC2PP:  packed matrix of the 3h1p/2h block. contains also diags
c  ADC2PF:  full packed ADC matrix in the conventional nonzero element
c           storage. This file is created in this routine. Important:
c           In contrast to the one-particle case the diagonal elements
c           are already included in this matrix !
c  ADCDAT:  configuration information for the 2h and 3h1p vector
C
C---------------Executable code--------------------------------------
c
c
c  say hello
c

      CALL PST('Building real SAT block+')
      WRITE(IW,*) 'Threshold for ADC matrix:',ADCTHR

c  create master index table now with *absolute* spinor indices!
c  This index table assigns to a row or column number within the
c  satellite block a quadruple of absolute spinor indices in order
c  to enable the treatment of all the delta conditions.
c  The counting starts with 1 and means the first entry of the
c  3h1p/3h1p block not the start of the total ADC matrix!
c
c  take care with the conjuga routine! This routine *really* assumes a 
c  *complex argument*. It must not be applied to real arrays with the
c  implicit assumption that it will not affect anything.
c

      IXTOT = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            RRP = MULTB(DESREP+NREP,IJKRP,2)
            DO R = 1,NV(RRP)
              IRCD(1,IXTOT) = IV(RRP) + R
              IRCD(2,IXTOT) = I
              IRCD(3,IXTOT) = J
              IRCD(4,IXTOT) = K
              IXTOT = IXTOT + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF( (IXTOT-1) .NE. IRCDLEN) stop 'internal error in BUILD_SAT!'
c
c  read in OOOO and VOVO integrals and cc the OOOO integrals
c
      CALL GETOOOO(OOOO(1))
      CALL GETVOVO(VOVO(1))
c
c  insert the neutral element for addition
c
      OOOO(0)=0.0d0
      VOVO(0)=0.0d0
c
c  create symmetry-specific filenames for 
c  sat-diag eigenvalues, for packed 3h1p/2h matrix 
c  ATT: final packed ADC matrix now only described by generic name!
c  These **huge** files will be deleted immediately after the
c  diagonalization.
c
      IF(DESREP.GT.9) THEN
        WRITE(FNAME2,'(A6,A1,I2)') NAME2,'.',DESREP
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
      ELSE
        WRITE(FNAME2,'(A6,A2,I1)') NAME2,'.0',DESREP
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
      ENDIF

      WRITE(IW,*) '2P-ADC matrix satellite energies in ',FNAME2
      WRITE(IW,*) '2P-ADC 3h1p/2h (packed) data in     ',FNAME3
      WRITE(IW,*) '2P-ADC final packed ADC matrix in   ',NAME4

      ITAP2 = ITAPADC+2
      ITAP3 = ITAPADC+3
      ITAP4 = ITAPADC+4

      OPEN(ITAP2,FILE=FNAME2,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP4,FILE=NAME4, FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c  read packed 3h1p/2h packed matrix created previously and transform
c  it to final format where only the nonzero numbers together with
c  row and column indices are stored.
c
      IB = 0
      NBUFS = 0
      JDUMMY = 0
      NTOTELE = 0
c
c  the column switching is implicitly done. Only the row index was
c  stored. The Column and row indices are absolute, no offset is
c  to be added in this case.
c
c
c  *ATT* *ATT* *ATT*  DSYEV yields eigenvalues sorted from the lowest
c  to the highest. Lanczos is not showing this behavior.
c  It converges fastest to the extremum eigenvalue which would be in our
c  case the highest possible two-hole final state and not the desired
c  lowest one!! Writing out the negative of the ADC matrix leads to the
c  desired behavior for Lanczos!
c
      SATCOL=1
 177  CONTINUE
      READ(ITAP3) SATROW,G
      IF(SATROW.EQ.-1) THEN
        SATCOL = SATCOL+1
        GOTO 177
      ENDIF
      IF(SATROW.EQ.-255) GOTO 180
      IF(G.EQ.0.0D0) STOP 'Read zero element in BUILD_SAT_R'
      IF(SATROW.LT.SATCOL) 
     &     STOP 'ROW/COL inconsistency in BUILD_SAT_R'
c
c  write matrix element only if greater than ADCTHR
c
      IF(ABS(G).GT.ADCTHR) THEN
        NTOTELE = NTOTELE + 1
        IB=IB+1
c
        BUF(IB)=G
        IOI(IB)=SATROW
        IOJ(IB)=SATCOL
        IF(IB.EQ.INTBUF) THEN
          WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                 INTBUF,JDUMMY
          IB = 0
          NBUFS = NBUFS + 1
        ENDIF
      ENDIF
      GOTO 177   !next element
 180  CONTINUE
      CLOSE(ITAP3,STATUS='DELETE')

c
c  from here the output buffer may remain partially filled
c  because more elements are about to come down below.
c

c
c  read in diagonal satellite eigenvalues (SATEV)
c
      READ(ITAP2) (SATEV(IZ),IZ=1,IRCDLEN)
      CLOSE(ITAP2,STATUS='DELETE')
c
c  loop through the columns and check each entry
c  The index conventions are as in the accomp. material
c
      SATOFF = NOOT(DESREP)
      DO 401 SATCOL=1,IRCDLEN

        S   = IRCD(1,SATCOL)
        L   = IRCD(2,SATCOL)
        M   = IRCD(3,SATCOL)
        N   = IRCD(4,SATCOL)

        DO 400 SATROW=SATCOL,IRCDLEN

          R   = IRCD(1,SATROW)
          I   = IRCD(2,SATROW)
          J   = IRCD(3,SATROW)
          K   = IRCD(4,SATROW)
c
c  initialize column entry
c
          G=0.0d0
c
c  start the delta blocks
c
c **** D_rs block
c
          IF(R.EQ.S) THEN

            IF(I.EQ.L) THEN
              G = G + OOOO(LTO(j,k,m,n))
            ENDIF
            IF(I.EQ.M) THEN
              G = G - OOOO(LTO(j,k,l,n))
            ENDIF
            IF(I.EQ.N) THEN
              G = G + OOOO(LTO(j,k,l,m))
            ENDIF

            IF(J.EQ.L) THEN
              G = G - OOOO(LTO(i,k,m,n))
            ENDIF
            IF(J.EQ.M) THEN
              G = G + OOOO(LTO(i,k,l,n))
            ENDIF
            IF(J.EQ.N) THEN
              G = G - OOOO(LTO(i,k,l,m))
            ENDIF

            IF(K.EQ.L) THEN
              G = G + OOOO(LTO(i,j,m,n))
            ENDIF
            IF(K.EQ.M) THEN
              G = G - OOOO(LTO(i,j,l,n))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + OOOO(LTO(i,j,l,m))
            ENDIF

          ENDIF
c
c **** D_occ,occ block, alpha
c
          IF(J.EQ.M) THEN
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, beta
c
          IF(J.EQ.L) THEN
            IF(K.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,i))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, gamma
c
          IF(I.EQ.M) THEN
            IF(J.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,k))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,l,s,j))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, delta
c
          IF(I.EQ.L) THEN
            IF(J.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,k))
            ENDIF
            IF(J.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,k))
            ENDIF
            IF(K.EQ.M) THEN
              G = G + VOVO(LTV(r,n,s,j))
            ENDIF
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,m,s,j))
            ENDIF
          ENDIF
c
c check if diagonal entry has to be added
c
          IF(SATROW.EQ.SATCOL) THEN
            G = G + SATEV(SATROW)
          ENDIF
c
c fill this element in the output buffer and flush, if necessary.
c the buffer counter is set outside the row/column loop
c reminder: satrow and satcol count from 1 to length of sat block,
c but the matrix starts with 2h2h block ==> add offset to these
c variables.
c
          IF(ABS(G).GT.ADCTHR) THEN
            NTOTELE = NTOTELE + 1
            IB=IB+1
c
            BUF(IB)=G
            IOI(IB)=SATROW + SATOFF
            IOJ(IB)=SATCOL + SATOFF
          
            IF(IB.EQ.INTBUF) THEN
              WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                     (IOI(IXX),IXX=1,INTBUF),
     &                     (IOJ(IXX),IXX=1,INTBUF),
     &                     INTBUF,JDUMMY
              IB = 0
              NBUFS = NBUFS + 1
            ENDIF
          ENDIF


 400    CONTINUE
 401  CONTINUE
c
c  write out *remaining* elements of the file buffer
c  and count all nonzero elements in the ADC matrix.
c
      IF(IB.GT.0) THEN
        NBUFS = NBUFS + 1
        WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &               (IOI(IXX),IXX=1,INTBUF),
     &               (IOJ(IXX),IXX=1,INTBUF),
     &               IB,JDUMMY
      ENDIF
      CLOSE(ITAP4)
      WRITE(IW,*)
      WRITE(IW,*) 'Number of elements above threshold: ',NTOTELE
      WRITE(IW,*) 'Number of buffers in ADC matrix file   : ',NBUFS
      WRITE(IW,*) 'Dimension of ADC matrix in this symmetry:',
     6            IRCDLEN+SATOFF
      XNZEL=DBLE(NTOTELE)
      XTOTEL=DBLE(IRCDLEN+SATOFF)
      XTOTEL=XTOTEL*XTOTEL
      SPARS=(XNZEL/XTOTEL)*100.0D0
      SPARS=100.0D0 - SPARS
      WRITE(IW,'(A,F5.2,A)') ' Sparsity of ADC matrix:',SPARS,' %'
      WRITE(IW,*)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUILD_SAT_C(OOOO,VOVO,IROOO,MAXOCCSP,MAXVIRSP,
     &                       SATEV,IRCD,IRCDLEN,LTO,LTV,DESREP,
     &                       ADCTHR,IOI,IOJ,BUF,INTBUF)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes 3H1P/3H1P block of the Double ionization ADC Matrix
C     real number implementation. Attention: row/column matrix storage
C     happens via packed integers ! One can safely assume that the
C     matrix dimension of the ADC matrix never will exceed 2^31-1 !!
C
C---------------Calling variables--------------------------------------
C
      COMPLEX*16 OOOO(0:*),VOVO(0:*)
      INTEGER IROOO(MAXOCCSP,2),MAXOCCSP,MAXVIRSP
      INTEGER IRCD(4,IRCDLEN),IRCDLEN,DESREP
      INTEGER LTO(MAXOCCSP,MAXOCCSP,MAXOCCSP,MAXOCCSP)
      INTEGER LTV(MAXVIRSP,MAXOCCSP,MAXVIRSP,MAXOCCSP)
      REAL*8 SATEV(IRCDLEN)
      REAL*8 ADCTHR
c for the packed matrix output file
      INTEGER INTBUF
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      COMPLEX*16 BUF(INTBUF)
c
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/files.inc"

      COMMON/RDMATND/NBUFS
C
C---------------Local variables--------------------------------------
C
      COMPLEX*16 G
      INTEGER IB,NBUFS,JDUMMY,NTOTELE
      CHARACTER*6 NAME2,NAME3,NAME4
      PARAMETER(NAME2='ADC2PE',NAME3='ADC2PP',NAME4='ADC2PF')
      CHARACTER*9 FNAME2,FNAME3,FNAME4
      INTEGER R,I,J,K,S,L,M,N
      CHARACTER*44 FIELD

      REAL*8 SPARS,XNZEL,XTOTEL
c
c  meaning of the file names:
c  ADC2PE:  contains the 3h1p/3h1p diagonal energies
c  ADC2PP:  packed matrix of the 3h1p/2h block. contains also diags
c  ADC2PF:  full packed ADC matrix in the conventional nonzero element
c           storage. This file is created in this routine. Important:
c           In contrast to the one-particle case the diagonal elements
c           are already included in this matrix !
c  ADCDAT:  configuration information for the 2h and 3h1p vector
C
C---------------Executable code--------------------------------------
c
c
c  say hello
c

      CALL PST('Building complex SAT block+')
      WRITE(IW,*) 'Threshold for ADC matrix:',ADCTHR

c  create master index table now with *absolute* spinor indices!
c  This index table assigns to a row or column number within the
c  satellite block a quadruple of absolute spinor indices in order
c  to enable the treatment of all the delta conditions.
c  The counting starts with 1 and means the first entry of the
c  3h1p/3h1p block not the start of the total ADC matrix!
c
c  take care with the conjuga routine! This routine *really* assumes a 
c  *complex argument*. It must not be applied to real arrays with the
c  implicit assumption that it will not affect anything.
c

      IXTOT = 1
      DO K=1,MAXOCCSP
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            RRP = MULTB(DESREP+NREP,IJKRP,2)
            DO R = 1,NV(RRP)
              IRCD(1,IXTOT) = IV(RRP) + R
              IRCD(2,IXTOT) = I
              IRCD(3,IXTOT) = J
              IRCD(4,IXTOT) = K
              IXTOT = IXTOT + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF( (IXTOT-1) .NE. IRCDLEN) stop 'internal error in BUILD_SAT!'
c
c  read in OOOO and VOVO integrals and cc the OOOO integrals
c  according to the equations on paper the VOVOs must not be
c  changed.
c
      CALL GETOOOO(OOOO(1))
      CALL GETVOVO(VOVO(1))
      DO K=1,IOOOOTT(NREP+1)
        OOOO(K)=DCONJG(OOOO(K))
      ENDDO
c
c  insert the neutral element for addition
c
      OOOO(0)=(0.0d0,0.0d0)
      VOVO(0)=(0.0d0,0.0d0)
c
c  create symmetry-specific filenames for
c  sat-diag eigenvalues, for packed 3h1p/2h matrix
c  ATT: final packed ADC matrix now only described by generic name!
c  These **huge** files will be deleted immediately after the
c  diagonalization.
c
      IF(DESREP.GT.9) THEN
        WRITE(FNAME2,'(A6,A1,I2)') NAME2,'.',DESREP
        WRITE(FNAME3,'(A6,A1,I2)') NAME3,'.',DESREP
      ELSE
        WRITE(FNAME2,'(A6,A2,I1)') NAME2,'.0',DESREP
        WRITE(FNAME3,'(A6,A2,I1)') NAME3,'.0',DESREP
      ENDIF

      WRITE(IW,*) '2P-ADC matrix satellite energies in ',FNAME2
      WRITE(IW,*) '2P-ADC 3h1p/2h (packed) data in     ',FNAME3
      WRITE(IW,*) '2P-ADC final packed ADC matrix in   ',NAME4

      ITAP2 = ITAPADC+2
      ITAP3 = ITAPADC+3
      ITAP4 = ITAPADC+4

      OPEN(ITAP2,FILE=FNAME2,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP3,FILE=FNAME3,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(ITAP4,FILE=NAME4, FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c  read packed 3h1p/2h packed matrix created previously and transform
c  it to final format where only the nonzero numbers together with
c  row and column indices are stored.
c
      IB = 0
      NBUFS = 0
      JDUMMY = 0
      NTOTELE = 0
c
c  the column switching is implicitly done. Only the row index was
c  stored. The Column and row indices are absolute, no offset is
c  to be added in this case.
c
c
c  *ATT* *ATT* *ATT*  DSYEV yields eigenvalues sorted from the lowest
c  to the highest. Lanczos is not showing this behavior. It converges to
c  the highest eigenvalue first. This is not the behavior of Lanczos!
c  It converges fastest to the extremum eigenvalue which would be in our
c  case the highest possible two-hole final state and not the desired
c  lowest one!! Writing out the negative of the ADC matrix leads to the
c  desired behavior for Lanczos!
c
      SATCOL=1
 177  CONTINUE
      READ(ITAP3) SATROW,G
      IF(SATROW.EQ.-1) THEN
        SATCOL = SATCOL+1
        GOTO 177
      ENDIF
      IF(SATROW.EQ.-255) GOTO 180
      IF(G.EQ.(0.0D0,0.0D0)) STOP 'Read zero element in BUILD_SAT_C'
      IF(SATROW.LT.SATCOL) 
     &     STOP 'ROW/COL inconsistency in BUILD_SAT_C!'
      IF(ABS(G).GT.ADCTHR) THEN
        NTOTELE = NTOTELE + 1
        IB=IB+1
        BUF(IB)=G
        IOI(IB)=SATROW
        IOJ(IB)=SATCOL
        IF(IB.EQ.INTBUF) THEN
          WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                 (IOI(IXX),IXX=1,INTBUF),
     &                 (IOJ(IXX),IXX=1,INTBUF),
     &                 INTBUF,JDUMMY
          IB = 0
          NBUFS = NBUFS + 1
        ENDIF
      ENDIF
      GOTO 177   !next element
 180  CONTINUE
      CLOSE(ITAP3,STATUS='DELETE')

c
c  from here the output buffer may remain partially filled
c  because more elements are about to come down below.
c

c
c  read in diagonal satellite eigenvalues (SATEV)
c
      READ(ITAP2) (SATEV(IZ),IZ=1,IRCDLEN)
      CLOSE(ITAP2,STATUS='DELETE')
c
c  loop through the columns and check each entry
c  The index conventions are as in the accomp. material
c
      SATOFF = NOOT(DESREP)
      DO 401 SATCOL=1,IRCDLEN

        S   = IRCD(1,SATCOL)
        L   = IRCD(2,SATCOL)
        M   = IRCD(3,SATCOL)
        N   = IRCD(4,SATCOL)

        DO 400 SATROW=SATCOL,IRCDLEN

          R   = IRCD(1,SATROW)
          I   = IRCD(2,SATROW)
          J   = IRCD(3,SATROW)
          K   = IRCD(4,SATROW)
c
c  initialize column entry
c
          G=(0.0d0,0.0d0)
c
c  start the delta blocks
c
c **** D_rs block
c
          IF(R.EQ.S) THEN

            IF(I.EQ.L) THEN
              G = G + OOOO(LTO(j,k,m,n))
            ENDIF
            IF(I.EQ.M) THEN
              G = G - OOOO(LTO(j,k,l,n))
            ENDIF
            IF(I.EQ.N) THEN
              G = G + OOOO(LTO(j,k,l,m))
            ENDIF

            IF(J.EQ.L) THEN
              G = G - OOOO(LTO(i,k,m,n))
            ENDIF
            IF(J.EQ.M) THEN
              G = G + OOOO(LTO(i,k,l,n))
            ENDIF
            IF(J.EQ.N) THEN
              G = G - OOOO(LTO(i,k,l,m))
            ENDIF

            IF(K.EQ.L) THEN
              G = G + OOOO(LTO(i,j,m,n))
            ENDIF
            IF(K.EQ.M) THEN
              G = G - OOOO(LTO(i,j,l,n))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + OOOO(LTO(i,j,l,m))
            ENDIF

          ENDIF
c
c **** D_occ,occ block, alpha
c
          IF(J.EQ.M) THEN
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, beta
c
          IF(J.EQ.L) THEN
            IF(K.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,i))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,i))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, gamma
c
          IF(I.EQ.M) THEN
            IF(J.EQ.N) THEN
              G = G - VOVO(LTV(r,l,s,k))
            ENDIF
            IF(K.EQ.N) THEN
              G = G + VOVO(LTV(r,l,s,j))
            ENDIF
          ENDIF
c
c **** D_occ,occ block, delta
c
          IF(I.EQ.L) THEN
            IF(J.EQ.M) THEN
              G = G - VOVO(LTV(r,n,s,k))
            ENDIF
            IF(J.EQ.N) THEN
              G = G + VOVO(LTV(r,m,s,k))
            ENDIF
            IF(K.EQ.M) THEN
              G = G + VOVO(LTV(r,n,s,j))
            ENDIF
            IF(K.EQ.N) THEN
              G = G - VOVO(LTV(r,m,s,j))
            ENDIF
          ENDIF
c
c according to the formulation on paper the c.c. satellite
c block entry is created here. We therefore have to c.c. the 
c entry.
c
          G=DCONJG(G)
c
c check if diagonal entry has to be added
c
          IF(SATROW.EQ.SATCOL) THEN
            G = G + DCMPLX(SATEV(SATROW),0.0D0)
          ENDIF
c
c fill this element in the output buffer and flush, if necessary.
c the buffer counter is set outside the row/column loop
c reminder: satrow and satcol count from 1 to length of sat block,
c but the matrix starts with 2h2h block ==> add offset to these
c variables.
c
          IF(ABS(G).GT.ADCTHR) THEN
            NTOTELE = NTOTELE + 1
            IB=IB+1
c
c   !!!
c
            BUF(IB)=G
            IOI(IB)=SATROW + SATOFF
            IOJ(IB)=SATCOL + SATOFF
          
            IF(IB.EQ.INTBUF) THEN
              WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &                     (IOI(IXX),IXX=1,INTBUF),
     &                     (IOJ(IXX),IXX=1,INTBUF),
     &                     INTBUF,JDUMMY
              IB = 0
              NBUFS = NBUFS + 1
            ENDIF
          ENDIF


 400    CONTINUE
 401  CONTINUE
c
c  write out *remaining* elements of the file buffer
c  and count all nonzero elements in the ADC matrix.
c
      IF(IB.GT.0) THEN
        NBUFS = NBUFS + 1
        WRITE(ITAP4) (BUF(IXX),IXX=1,INTBUF),
     &               (IOI(IXX),IXX=1,INTBUF),
     &               (IOJ(IXX),IXX=1,INTBUF),
     &               IB,JDUMMY
      ENDIF
      CLOSE(ITAP4)
      WRITE(IW,*)
      WRITE(IW,*) 'Total of nonzero elements in ADC matrix: ',NTOTELE
      WRITE(IW,*) 'Number of buffers in ADC matrix file   : ',NBUFS
      WRITE(IW,*) 'Dimension of ADC matrix in this symmetry:',
     6            IRCDLEN+SATOFF
      XNZEL=DBLE(NTOTELE)
      XTOTEL=DBLE(IRCDLEN+SATOFF)
      XTOTEL=XTOTEL*XTOTEL
      SPARS=(XNZEL/XTOTEL)*100.0D0
      SPARS=100.0D0 - SPARS
      WRITE(IW,'(A,F5.2,A)') ' Sparsity of ADC matrix:',SPARS,' %'
      WRITE(IW,*)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_EIJKL (EPS,E)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Create energy array 0.5*(ei + ej + ek + el) in the same order as a
C     <OO||OO> integral stream is stored !
C     this energy array is always real, no complexation considerd.
C     ATT: one factor of 0.5 is absorbed due to summation r,s -->  r>s !
C
C    **** important remark ***
C
C  is is very error-prone to change all the indices from IJ to
C  KL and from AB to IJ because this is the denomvvoo source!
C  We therefore are lazy but safe and *interpret* just for this case
C  the AB variables as IJ and the IJ variable as KL. Therefore AB are
C  now *occupieds*  !!! We just adjust the energy offsets and signs (all
C  energies are *positive*), everything else stays identical
C  due to the twofold triangular storage.
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 E(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC,FAC1,FAC2,FAC3
C
C---------------Executable code--------------------------------------
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         FAC1 = + EPS(JJ)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            II = IOFF + I
            FAC2 = + EPS(II) + FAC1
            DO 20 BRP = 1, NREP
            BB = IO(BRP)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IO(ARP)
            DO B = 1, NO(BRP)
               BB = BB + 1
               FAC3 = FAC2 + EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NO(ARP)
                  AA = AOFF + A
                  FAC = FAC3 + EPS(AA)
                  ABIJ = ABIJ + 1

                  E(ABIJ) = 0.5D0 * FAC

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      IF(ABIJ.NE.IOOOOTT(NREP+1)) CALL QUIT('Error in MAKE_EIJKL')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DENOMY (EPS,T2,CT2)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Multiply  VVOO integrals by -(e_a + e_b) where the energies
C     comprise the virtuals only !
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      REAL*8 T2(*)
      COMPLEX*16 CT2(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 FAC3,FAC
      COMPLEX*16 CFAC
C
C---------------Executable code--------------------------------------
C
      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      JJ = IO(JRP)
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      IOFF = IO(IRP)
      DO J = 1, NO(JRP)
         JJ = JJ + 1
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            BB = IV(BRP) + IO(NREP+1)
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            AOFF = IV(ARP) + IO(NREP+1)
            DO B = 1, NV(BRP)
               BB = BB + 1
               FAC3 = EPS(BB)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NV(ARP)
                  AA = AOFF + A
                  FAC = -(FAC3 + EPS(AA))
                  ABIJ = ABIJ + 1
                  IF (CARITH) THEN
                     CFAC = DCMPLX(FAC,0.0D0)
                     CT2(ABIJ) = CT2(ABIJ)*CFAC
                  ELSE
                     T2(ABIJ) = T2(ABIJ)*FAC
                  ENDIF
               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE AWEAVERR(EPS,G,A)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C Scan the totally antisymmetric and triagularly stored Gamma_ijkl
C matrix and put at the proper locations the A_ik terms. A brief check
C reveals that G_ijkl is also totally antisymmetric with respect to the
C A_ik terms. At the same time we can "weave into" 
c the diagonal elements e_i + e_j to the proper positions.   
C
C remark: also derived from DENOMVVOO. therefore AB = IJ and IJ = KL !!
C
C          real version
C
C remark: it is easier to write the real/complex version separately.
C even if there is more code, it is for the sake of clarity avoiding a
C lot more ifs!
C Mention: all four indices refer to *occupied* orbitals!
C
C---------------Calling variables--------------------------------------
C
      REAL*8     EPS(*)
      REAL*8     G(*)
      REAL*8     A(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
       REAL*8 ESUM
       INTEGER IAOFF
C
C---------------Executable code--------------------------------------
C

      IJKL = 0
      DO KLRP = 1, NREP
      DO 10 LRP = 1, NREP
      KRP = MULTB(LRP,KLRP+NREP,2)
      IF (KRP.LT.LRP) GOTO 10
      DO L = 1, NO(LRP)
         KMIN = 1
         IF (KRP.EQ.LRP) KMIN = L + 1
         DO K = KMIN, NO(KRP)
            DO 20 JRP = 1, NREP
            IRP = MULTB(JRP,KLRP+NREP,2)
            IF (IRP.LT.JRP) GOTO 20
            DO J = 1, NO(JRP)
               IMIN = 1
               IF (IRP.EQ.JRP) IMIN = J + 1
               DO I = IMIN, NO(IRP)

                  IJKL = IJKL + 1
c
c from here we know lrep,krep,jrep,irep  AND  I,J,K,L
c we collect the A_ik terms and add them properly with *CORRECT SIGNS*
c IMPORTANT FOR DEBUGGING: From Abelian group theory one can derive
c that if e.g. JRP.eq.LRP then IRP must also be KRP!
c
c  debugging: this check is included and does not eat up a lot of time
c  because the 2h/2h block is small.
c
                  IF(JRP.eq.LRP) THEN   !positive sign
                    IF(J.eq.L) THEN
                      if(irp.ne.krp) call quit('weaver_error')
                      IAOFF = IOO(KRP) + (I-1)*NO(KRP) + K
                      G(IJKL) = G(IJKL) + A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  

                  IF(IRP.eq.LRP) THEN   !negative sign
                    IF(I.eq.L) THEN
                      if(jrp.ne.krp) call quit('weaver_error')
                      IAOFF = IOO(KRP) + (J-1)*NO(KRP) + K
                      G(IJKL) = G(IJKL) - A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  

                  IF(JRP.eq.KRP) THEN   !negative sign
                    IF(J.eq.K) THEN
                      if(irp.ne.lrp) call quit('weaver_error')
                      IAOFF = IOO(LRP) + (I-1)*NO(LRP) + L
                      G(IJKL) = G(IJKL) - A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  

                  IF(IRP.eq.KRP) THEN   !positive sign
                    IF(I.eq.K) THEN
                      if(jrp.ne.lrp) call quit('weaver_error')
                      IAOFF = IOO(LRP) + (J-1)*NO(LRP) + L
                      G(IJKL) = G(IJKL) + A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  
c
C  treat the diagonal element
c
             IF(IRP.eq.KRP.and.JRP.eq.LRP.and.i.eq.k.and.j.eq.l) THEN
                IOFF = IO(IRP) + I
                JOFF = IO(JRP) + J
                ESUM = EPS(IOFF) + EPS(JOFF)
                G(IJKL) = G(IJKL) - ESUM 
             ENDIF

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE AWEAVERC(EPS,G,A)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C  see description in routine AWEAVERR, but here we have the
C
C          complex version
C
C---------------Calling variables--------------------------------------
C
      REAL*8         EPS(*)
      COMPLEX*16     G(*)
      COMPLEX*16     A(*)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
       REAL*8 ESUM
       INTEGER IAOFF
C
C---------------Executable code--------------------------------------
C

      IJKL = 0
      DO KLRP = 1, NREP
      DO 10 LRP = 1, NREP
      KRP = MULTB(LRP,KLRP+NREP,2)
      IF (KRP.LT.LRP) GOTO 10
      DO L = 1, NO(LRP)
         KMIN = 1
         IF (KRP.EQ.LRP) KMIN = L + 1
         DO K = KMIN, NO(KRP)
            DO 20 JRP = 1, NREP
            IRP = MULTB(JRP,KLRP+NREP,2)
            IF (IRP.LT.JRP) GOTO 20
            DO J = 1, NO(JRP)
               IMIN = 1
               IF (IRP.EQ.JRP) IMIN = J + 1
               DO I = IMIN, NO(IRP)

                  IJKL = IJKL + 1
c
c from here we know lrep,krep,jrep,irep  AND  I,J,K,L
c we collect the A_ik terms and add them properly with *CORRECT SIGNS*
c IMPORTANT FOR DEBUGGING: From Abelian group theory one can derive
c that if e.g. JRP.eq.LRP then IRP must also be KRP!
c
c
                  IF(JRP.eq.LRP) THEN   !positive sign
                    IF(J.eq.L) THEN
                      if(irp.ne.krp) call quit('weaver_error')
                      IAOFF = IOO(KRP) + (I-1)*NO(KRP) + K
                      G(IJKL) = G(IJKL) + A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  

                  IF(IRP.eq.LRP) THEN   !negative sign
                    IF(I.eq.L) THEN
                      if(jrp.ne.krp) call quit('weaver_error')
                      IAOFF = IOO(KRP) + (J-1)*NO(KRP) + K
                      G(IJKL) = G(IJKL) - A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  

                  IF(JRP.eq.KRP) THEN   !negative sign
                    IF(J.eq.K) THEN
                      if(irp.ne.lrp) call quit('weaver_error')
                      IAOFF = IOO(LRP) + (I-1)*NO(LRP) + L
                      G(IJKL) = G(IJKL) - A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  

                  IF(IRP.eq.KRP) THEN   !positive sign
                    IF(I.eq.K) THEN
                      if(jrp.ne.lrp) call quit('weaver_error')
                      IAOFF = IOO(LRP) + (J-1)*NO(LRP) + L
                      G(IJKL) = G(IJKL) + A(IAOFF)
c            write(*,*) 'weaver active:',IRP,JRP,KRP,LRP,I,J,K,L
                    ENDIF
                  ENDIF                  
c
c this is the check for the diagonal element
c
             IF(IRP.eq.KRP.and.JRP.eq.LRP.and.i.eq.k.and.j.eq.l) THEN
                IOFF = IO(IRP) + I
                JOFF = IO(JRP) + J
                ESUM = EPS(IOFF) + EPS(JOFF)
                G(IJKL) = G(IJKL) - DCMPLX(ESUM,0.0D0) 
             ENDIF

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EWARRMULR (RA,RE,NUM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Multiplies the real array RA elementwise with the real array RE
C     and stores the result in RE
C
C---------------Calling parameters ------------------------------------
C
      INTEGER NUM
      REAL*8 RA(NUM),RE(NUM)

      DO I=1,NUM
        RE(I) = RE(I)*RA(I)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE EWARRMULC (RA,RE,NUM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Multiplies the complex vector RE with the real vector RA
c     and stores the result in RE
C
C---------------Calling parameters ------------------------------------
C
      INTEGER NUM
      REAL*8 RA(NUM)
      COMPLEX*16 RE(NUM)

      DO I=1,NUM
        RE(I) = RE(I)*DCMPLX(RA(I),0.0D0)
      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKE_AIK(AIK14,VVOO1,VVOO2,EPS,EABLKKS,KREP)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Computes the A_ik elements required for term 1-4 in the formulas
C
C---------------Calling variables--------------------------------------
C
      REAL*8 AIK14(*),VVOO1(*),VVOO2(*),EPS(*),EABLKKS(*)
      INTEGER KREP
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/complex.inc"
#include  "../relccsd/param.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
C  Create the v_rskm integrals (i.e. the VVOO integrals divided
C  by e_r + e_s - e_k - e_m). v_rskm and v*_rsim are afterwards
C  available as IJK,L:LREP in VVOO1 and VVOO2

      CALL GETVVOO(VVOO1)
      CALL DENOMVVOO(EPS,VVOO1,VVOO1)
      CALL SRT1T3 (NREP,MULTB,LFA,NVVT,NO,NO,NVVOT,KVVOOT,KKVVOT,
     &             VVOO1,VVOO2)
      NX=KVVOOT(NREP+1)
      CALL XCOPY(NX,VVOO2,1,VVOO1,1)
      IF(CARITH) THEN
        CALL CONJUGA(NX,VVOO2,1)
      ENDIF
C
C  Generate the epsilon array of VVOT,O,O structure in the symmetry KREP
C
      CALL EPSARR(EPS,EABLKKS,KREP,ISIZE)
C
C  contract the arrays v_rsm,k and v*_rsm,i and E_rsm,ik in symmetry
C  KREP over the common index range rsm.
C  The resulting A_ik block has the BKC row: |>  col: <| =   |><|
C
      CALL CNTRCT3(VVOO1,VVOO2,EABLKKS,AIK14,
     &             NVVOT(KREP),NO(KREP),KVVOOT(KREP))
c
c  negate the result because CNTRCT3 produces a negative result
c  ATT the AIK(*) is declared as REAL! We cannot use an operation
c  which assumes AIK to be complex! Special operation routine!
c
      NX = NO(KREP)*NO(KREP)
      CALL XSCAL(NX,-A1,AIK14,1)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKELOOKUPO(LT,LENLT)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C    **** important remark ***
C
C  is is very error-prone to change all the indices from IJ to
C  KL and from AB to IJ because this is the denomvvoo source!
C  We therefore are lazy but safe and *interpret* just for this case
C  the AB variables as IJ and the IJ variable as KL. Therefore AB are
C  now *occupieds*  !!! 
C
C---------------Calling variables--------------------------------------
C
      INTEGER LT(LENLT,LENLT,LENLT,LENLT)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER AABSL,BABSL,IABSL,JABSL
C
C---------------Executable code--------------------------------------
C
c
c  initialize lookup table (extremely important due to null offsets)
c  loop structure follows storage mode of the OOOO integrals.
c
      do i=1,lenlt
      do j=1,lenlt
      do k=1,lenlt
      do l=1,lenlt
        LT(l,k,j,i) = 0
      enddo
      enddo
      enddo
      enddo

c     WRITE(IW,*) 'Preparing OOOO Lookup table...'

      ABIJ = 0
      DO IJRP = 1, NREP
      DO 10 JRP = 1, NREP
      IRP = MULTB(JRP,IJRP+NREP,2)
      IF (IRP.LT.JRP) GOTO 10
      DO J = 1, NO(JRP)
         IMIN = 1
         IF (IRP.EQ.JRP) IMIN = J + 1
         DO I = IMIN, NO(IRP)
            DO 20 BRP = 1, NREP
            ARP = MULTB(BRP,IJRP+NREP,2)
            IF (ARP.LT.BRP) GOTO 20
            DO B = 1, NO(BRP)
               AMIN = 1
               IF (ARP.EQ.BRP) AMIN = B + 1
               DO A = AMIN, NO(ARP)
                  ABIJ = ABIJ + 1

                  AABSL = IO(ARP) + A
                  BABSL = IO(BRP) + B
                  IABSL = IO(IRP) + I
                  JABSL = IO(JRP) + J

                  LT(AABSL,BABSL,IABSL,JABSL)=ABIJ

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO
C
      IF(ABIJ.NE.IOOOOTT(NREP+1)) CALL QUIT('Error in MAKELOOKUPO')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKELOOKUPV(LT,LENLTO,LENLTV)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C    **** important remark ***
C
C  is is very error-prone to change all the indices from IJ to
C  KL and from AB to IJ because this is the denomvvoo source!
C  We therefore are lazy but safe and *interpret* just for this case
C  the AB variables as IJ and the IJ variable as KL. Therefore AB are
C  now *occupieds*  !!! 
C
C---------------Calling variables--------------------------------------
C
      INTEGER LT(LENLTV,LENLTO,LENLTV,LENLTO)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER AABSL,BABSL,IABSL,JABSL
C
C---------------Executable code--------------------------------------
C
c
c  initialize lookup table (extremely important due to null offsets)
c  loop structure follows storage mode of the OOOO integrals.
c
      do i=1,lenlto
      do j=1,lenltv
      do k=1,lenlto
      do l=1,lenltv
        LT(l,k,j,i) = 0
      enddo
      enddo
      enddo
      enddo

c     WRITE(IW,*) 'Preparing VOVO Lookup table...'

      AIBJ = 0
      DO BJRP = 1, NREP
      DO 10 JRP = 1, NREP
      BRP = MULTB(JRP,BJRP+NREP,2)
      DO J = 1, NO(JRP)
         DO B = 1, NV(BRP)
            DO 20 IRP = 1, NREP
            ARP = MULTB(IRP,BJRP+NREP,2)
            DO I = 1, NO(IRP)
               DO A = 1, NV(ARP)
                  AIBJ = AIBJ + 1

                  AABSL = IV(ARP) + A
                  IABSL = IO(IRP) + I
                  BABSL = IV(BRP) + B
                  JABSL = IO(JRP) + J

                  LT(AABSL,IABSL,BABSL,JABSL)=AIBJ

               ENDDO
            ENDDO
 20         CONTINUE
         ENDDO
      ENDDO
 10   CONTINUE
      ENDDO

      IF(AIBJ.NE.IVOVO(NREP+1)) CALL QUIT('Error in MAKELOOKUPV')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MAKESPTRA(IROOO,ILENO)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Creates the auxiliary array for conversion of relative/absolute
C     spinor numbers for occupied spinors (only these are needed)
C
C---------------Calling variables--------------------------------------
C
      INTEGER IROOO(ILENO,2)
C
C---------------Common Blocks--------------------------------------
C
#include  "../relccsd/symm.inc"
#include  "../relccsd/files.inc"
C
C---------------Executable code--------------------------------------
c
      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
      If((IX-1).ne.ileno) STOP 'Consistency violation in MAKESPTRA!'

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



