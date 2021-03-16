module extra

   use interface_ao
   use interface_mo

   implicit none

   public gtdoav
   public timtxt
   public reacmo
   public qtrans90

   private

contains

      SUBROUTINE TIMTXT(TEXT,TIMUSD,LUPRIN)
!
! TIMTXT based on TIMER by TUH //900709-hjaaj
!
      CHARACTER*(*) TEXT
      CHARACTER AHOUR*6, ASEC*8, AMIN*8
      real(8) :: timusd
      integer :: luprin, IHOURS, MINUTE, ISECND
!
      ISECND = NINT(TIMUSD)
      IF (ISECND .GE. 60) THEN
         MINUTE = ISECND/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = ISECND - 3600*IHOURS - 60*MINUTE
         IF (IHOURS .EQ. 1) THEN
            AHOUR = ' hour '
         ELSE
            AHOUR = ' hours'
         END IF
         IF (MINUTE .EQ. 1) THEN
            AMIN = ' minute '
         ELSE
            AMIN = ' minutes'
         END IF
         IF (ISECND .EQ. 1) THEN
            ASEC = ' second '
         ELSE
            ASEC = ' seconds'
         END IF
         IF (IHOURS .GT. 0) THEN
            WRITE(LUPRIN,100) &
     &            TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
         ELSE
            WRITE(LUPRIN,200) TEXT, MINUTE, AMIN, ISECND, ASEC
         END IF
      ELSE
         WRITE(LUPRIN,300) TEXT,TIMUSD
      END IF
  100 FORMAT(1X,A,I4,A,I3,A,I3,A)
  200 FORMAT(1X,A,     I3,A,I3,A)
  300 FORMAT(1X,A,F7.2,' seconds')
      end subroutine

         SUBROUTINE DZERO(DX,LENGTH)
!
! Last revision 5-May-1984 by Hans Jorgen Aa. Jensen
!
!   Subroutine DZERO sets a real array of length *LENGTH*
!   to zero.
!...................................................................
      real(8) :: DX(LENGTH)
      integer :: i, length
!
      IF (LENGTH.LE.0) RETURN
!
      DO 100 I = 1,LENGTH
  100    DX(I) = 0.0D00
!
      end subroutine

      SUBROUTINE REACMO(IUNIT,CMOFIL,CMO,EIG,IBEIG,TOTERG,IOPT)
!***********************************************************************
!
!     Read DHF-coefficients and eigenvalues from unformatted file
!     Read option is provided from bit-packed IOPT:
!       0001 - give restart info
!       0010 - read coefficients
!       0100 - read eigenvalues
!       1000 - read boson irrep identification
!
!     IDIM(1,IFRP) - number of positronic solutions
!     IDIM(2,IFRP) - number of electronic solutions
!     IDIM(3,IFRP) - number of AO-basis functions
!
!     Written by T.Saue Sept 1 1995
!
!***********************************************************************

!
      LOGICAL INFO, TEST
      logical, external :: mylbit
      CHARACTER TEXT*74,CMOFIL*6,FMT*6
      real(8) :: CMO((nr_mo_gerade*nr_ao_gerade + nr_mo_ungerade*nr_ao_ungerade)*nr_quaternion_blocks)
      real(8) :: EIG(nr_mo)
      integer :: IBEIG(nr_mo)
      integer :: IDIM(3,2)
      integer :: iunit
      real(8), allocatable :: kc(:)
      real(8), allocatable :: ko(:)
      real(8) :: toterg
      integer, allocatable :: ki(:)
      integer :: iopt
      integer :: ifrp, jc, i, j, jo, nc
      integer :: nsym, no, ne, np, ioff, ni
      integer :: nr_ao_in_fermion_corep(2), icmoq(2)
!
      REWIND IUNIT
!
!     Read title line
!     ===============
!
      READ (IUNIT,END=10,ERR=20) TEXT,NSYM,                             &
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG

!
!     File info
!     =========
!
      IF(mylbit(IOPT,1)) THEN
        WRITE(*,*)                                                      &
     &   '* REACMO: Coefficients read from file ',CMOFIL,               &
     &   ' - Total energy: ',TOTERG
        WRITE(*,'(2A)') '* Heading :',TEXT
      ENDIF
!
      IF(IOPT.LT.2) RETURN
!
!     Compatibility tests for number of orbitals
!     ==========================================
!
      nr_ao_in_fermion_corep(1) = nr_ao_gerade
      nr_ao_in_fermion_corep(2) = nr_ao_ungerade

      TEST = .FALSE.
      NC = 0
      NO = 0
      DO IFRP = 1,nr_fermion_coreps
             if (ifrp == 1) then
                np = nr_mo_gerade_negative_secondary
                ni = nr_mo_gerade_positive_inactive
                no = nr_mo_gerade
                ioff = 0
             else
                np = nr_mo_ungerade_negative_secondary
                ni = nr_mo_ungerade_positive_inactive
                no = nr_mo_ungerade
                ioff = nr_mo_gerade
             end if
        IDIM(3,IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
        IF(IDIM(3,IFRP).NE.NO) TEST=.TRUE.
        NC = NC + IDIM(3,IFRP)*nr_ao_in_fermion_corep(ifrp)*nr_quaternion_blocks
        NO = NO + IDIM(3,IFRP)
      ENDDO
!
!     Read coefficients
!     =================
!
      IF(mylbit(IOPT,2)) THEN
        IF(TEST) THEN
          allocate(kc(nc))
          CALL READT(IUNIT,NC, kc)
          JC = 1
          DO IFRP = 1,nr_fermion_coreps
             if (ifrp == 1) then
                np = nr_mo_gerade_negative_secondary
                ni = nr_mo_gerade_positive_inactive
                no = nr_mo_gerade
                ioff = 0
             else
                np = nr_mo_ungerade_negative_secondary
                ni = nr_mo_ungerade_positive_inactive
                no = nr_mo_ungerade
                ioff = nr_mo_gerade
             end if
             ne = no - np
             icmoq = 0
             if (nr_fermion_coreps > 1) then
                icmoq(2) = nr_mo_gerade*nr_ao_gerade
             end if
            CALL ADACMO(nr_ao_in_fermion_corep(ifrp),nr_quaternion_blocks,                               &
     &         CMO(ICMOQ(IFRP)+1),NO,NP,NE,     &
     &         kc(JC),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            JC = JC + IDIM(3,IFRP)*nr_ao_in_fermion_corep(ifrp)*nr_quaternion_blocks
          ENDDO
          deallocate(kc)
        TEST = .FALSE.
        ELSE
          READ (IUNIT,END=30,ERR=40) CMO
        ENDIF
        IF(IOPT.LT.4) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
!
!     Read eigenvalues
!     ================
!
      IF(mylbit(IOPT,3)) THEN
        IF(TEST) THEN
          allocate(ko(no))
          CALL READT(IUNIT,NO,ko)
          JO = 1
          DO IFRP = 1,nr_fermion_coreps
             if (ifrp == 1) then
                np = nr_mo_gerade_negative_secondary
                ni = nr_mo_gerade_positive_inactive
                no = nr_mo_gerade
                ioff = 0
             else
                np = nr_mo_ungerade_negative_secondary
                ni = nr_mo_ungerade_positive_inactive
                no = nr_mo_ungerade
                ioff = nr_mo_gerade
             end if
             ne = no - np
            CALL ADAEIG(                                                &
     &         EIG(Ioff+1),NO,NP,NE,      &
     &         ko(JO),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            JO = JO + IDIM(1,IFRP)+IDIM(2,IFRP)
          ENDDO
          deallocate(ko)
        ELSE
          READ (IUNIT,END=50,ERR=60) EIG
        ENDIF
        IF(IOPT.LT.8) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
      return
!
!     Read single group irrep (spinfree and Levy-Leblond calcs)
!     =========================================================
!
      IF(mylbit(IOPT,4)) THEN
        IF(TEST) THEN
          allocate(ki(no))
          read(iunit) ki
          CALL ADAIBP(IBEIG,ki,IDIM)
          deallocate(ki)
        ELSE
          READ (IUNIT,END=70,ERR=80) IBEIG
        ENDIF
        IF(IOPT.LT.16) RETURN
      ELSE
        READ(IUNIT)
      ENDIF
!
      RETURN
!
!     Error messages
!
 10   CONTINUE
      write(*,'(//A/A/A)') &
     &   'REACMO: END OF FILE reading first record on DFCOEF.', &
     &   '        Possible reason: DFCOEF was generated with'// &
     &      ' 32-bit integer Dirac version.', &
     &   '        This Dirac is compiled with 64 bit integers.'
      call my_quit('REACMO: END OF FILE reading first record of DFCOEF')
 20   CONTINUE
      write(*,'(//A/A/A)') &
     &   'REACMO: ERROR reading first record on DFCOEF.', &
     &   '        Possible reason: DFCOEF was generated with'// &
     &      ' 32-bit integer Dirac version.', &
     &   '        This Dirac is compiled with 64 bit integers.'
      call my_quit('REACMO: ERROR reading first record of DFCOEF')
 30   CONTINUE
      call my_quit('REACMO: END OF FILE reading coefficients')
 40   CONTINUE
      call my_quit('REACMO: ERROR reading coefficients')
 50   CONTINUE
      call my_quit('REACMO: END OF FILE reading eigenvalues')
 60   CONTINUE
!     call my_quit('REACMO: ERROR reading eigenvalues')
      WRITE(*,'(/A)') 'REACMO: WARNING !!! No eigenvalues read'
      RETURN
 70   CONTINUE
      call my_quit('REACMO: END OF FILE reading irrep identification')
 80   CONTINUE
      call my_quit('REACMO: ERROR reading irrep identification')
      end subroutine

      SUBROUTINE READT (IT,N,X)
      CHARACTER*30 STRING
      real(8) :: X(N)
      integer :: it, n
      READ (IT,END=10,ERR=20) X
      RETURN
 10   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: END reading unit  ',IT
      call my_quit(STRING)
 20   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: Error reading unit',IT
      call my_quit(STRING)
      end subroutine

      SUBROUTINE QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,              &
     &                  FAO,LRAO,LCAO,NZAO,IQAO,                        &
     &                  FMO,LRMO,LCMO,NZMO,IQMO,                        &
     &                  TM1,LR1,LC1,NZTM1,IQTM1,                        &
     &                  TM2,LR2,LC2,NZTM2,IQTM2,                        &
     &                  BUF,LBUF,IPRINT)
!***********************************************************************
!
!     This routine performs unitary transformations indicated by TYP:
!
!       TYP = AOMO : AO-to-MO-transformation  FMO = (C+)FAOC + FADD*FMO
!       TYP = MOAO : MO-to-AO-transformation  FAO = CFMO(C+) + FADD*FAO
!
!     FADD gives the possibility of adding the results
!
!     TREV indicates symmetry of F-matrix under timereversal
!
!       TREV = 'S' - symmetric     (t =  1)
!       TREV = 'A' - anti-symmtric (t = -1)
!
!     The unitary C-matrix has a time-symmetric structure.
!
!     The AO-to-MO-transformation can be set up quaternionically as
!
!       FMO = [(Ca+)-t(CbT)j][(FAOa)+(FAOb)j][(Ca)+(Cb)j]
!
!       F(NBAS,NBAS,NZ) --->  F(NORB,NORB,NZ)
!
!     and the MO-to-Ao-transformation can be set up as:
!
!       FAO = [(Ca)+t(Cb)j][(FMOa)+(FMOb)j][(Ca+)-(Cb)j]
!
!       F(NORB,NORB,NZ) --->  F(NBAS,NBAS,NZ)
!
!     Written by T.Saue,January 1995
!     Last revision: 29 January 1997 : Luuk Visscher
!**********************************************************************
!
      integer :: NRAO,NCAO,NRMO,NCMO
      integer :: LRAO,LCAO,NZAO
      integer :: LRMO,LCMO,NZMO
      integer :: LR1,LC1,NZTM1
      integer :: LR2,LC2,NZTM2
      integer :: iprint
      CHARACTER TREV*1,TRA*1,TYP*4
      real(8) :: FAO(LRAO,LCAO,*),                                      &
     &          FMO(LRMO,LCMO,*),                                       &
     &          TM1(LR1,LC1,*),                                         &
     &          TM2(LR2,LC2,*),                                         &
     &          BUF(LBUF)
      real(8) :: fadd
      integer :: lbuf
      integer :: IQAO(*),                                               &
     &          IQMO(*),                                                &
     &          IQTM1(*),                                               &
     &          IQTM2(*)
      integer :: IQW(4)
      integer :: nbuf, nzw

!
!
!     Add the results on to the contents of the MO matrix ?
!
!
!     Time reversal symmetry of F-matrix
!
      IF    (TREV.EQ.'S') THEN
        TRA = 'N'
      ELSEIF(TREV.EQ.'A') THEN
        TRA = 'I'
      ELSE
        call my_quit('QTRANS: Unknown keyword TREV '//TREV)
      ENDIF
!
!     *******************************************************
!     *****  AO - to - MO -transformation : F' = (C+)FC *****
!     *******************************************************
!
      IF(TYP.EQ.'AOMO') THEN
!
!       Look for the most efficient transformation
!       ------------------------------------------
!
        IF (NCAO.GE.NRAO) THEN
!
!         Perform the two steps:
!           1. W  = FC
!           2. F' = (C+)W
!
          CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (*,'(/A,I10/A,I10)')                                 &
     &   ' >>> QTRANS error, need work space    ',NBUF,                 &
     &   '                   current work space ',LBUF
          ENDIF

!
!         First part of transformation: W = FC
!         ------------------------------------
!
          CALL QGEMM(NRAO,NCMO,NCAO,1.0d0,                              &
     &              'N','N',IQAO,FAO,LRAO,LCAO,NZAO,                    &
     &              'N','N',IQTM2,TM2,LR2,LC2,NZTM2,                    &
     &               0.0d0,IQW,BUF,NRAO,NCMO,NZW)

!
!         Second part of transformation: H = (C+)W
!         ----------------------------------------
!
          CALL QGEMM(NRMO,NCMO,NRAO,1.0d0,                              &
     &               'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,                   &
     &               'N','N',IQW,BUF,NRAO,NCMO,NZW,                     &
     &               FADD,IQMO,FMO,LRMO,LCMO,NZMO)


        ELSE
!
!         Perform the two steps:
!           1. W  = (C+)F
!           2. F' = WC
!
          CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (*,'(/A,I10/A,I10)')                                 &
     &   ' >>> QTRANS error, need work space    ',NBUF,                 &
     &   '                   current work space ',LBUF
          ENDIF

!
!         First part of transformation: W = (C+)F
!         ----------------------------------------
!
          CALL QGEMM(NRMO,NCAO,NRAO,1.0d0,                              &
     &       'H',TRA,IQTM1,TM1,LR1,LC1,NZTM1,                           &
     &       'N','N',IQAO,FAO,LRAO,LCAO,NZAO,                           &
     &        0.0d0,IQW,BUF,NRMO,NCAO,NZW)
!
!         Second part of transformation: H = WC
!         -------------------------------------
!
          CALL QGEMM(NRMO,NCMO,NCAO,1.0d0,                              &
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,                            &
     &       'N','N',IQTM2,TM2,LR2,LC2,NZTM2,                           &
     &       FADD,IQMO,FMO ,LRMO,LCMO,NZMO)
!
        ENDIF
!
!     *******************************************************
!     *****  MO - to - AO -transformation : F' = CF(C+) *****
!     *******************************************************
!
      ELSEIF(TYP.EQ.'MOAO') THEN
!
!       Look for the most efficient transformation
!       ------------------------------------------
!
        IF(NCMO.GE.NRMO) THEN
          CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
          NBUF = NRMO*NCAO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (*,'(/A,I10/A,I10)')                                 &
     &   ' >>> QTRANS error, need work space    ',NBUF,                 &
     &   '                   current work space ',LBUF
          ENDIF
!
!         First part of transformation: W = F(C+)
!
          CALL QGEMM(NRMO,NCAO,NCMO,1.0d0,                              &
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,                           &
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,                           &
     &       0.0d0,IQW,BUF,NRMO,NCAO,NZW)
!
!         Second part of transformation: H = CW
!
          CALL QGEMM(NRAO,NCAO,NRMO,1.0d0,                              &
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,                          &
     &       'N','N',IQW ,BUF,NRMO,NCAO,NZW,                            &
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ELSE
          CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
          NBUF = NRAO*NCMO*NZW
          IF(NBUF.GT.LBUF) THEN
             WRITE (*,'(/A,I10/A,I10)')                                 &
     &   ' >>> QTRANS error, need work space    ',NBUF,                 &
     &   '                   current work space ',LBUF
          ENDIF
!
!         First part of transformation: W = CF
!
          CALL QGEMM(NRAO,NCMO,NRMO,1.0d0,                              &
     &       'N',TRA,IQTM1,TM1 ,LR1,LC1,NZTM1,                          &
     &       'N','N',IQMO,FMO,LRMO,LCMO,NZMO,                           &
     &       0.0d0,IQW,BUF,NRAO,NCMO,NZW)
!
!         Second part of transformation: H = CW
!
          CALL QGEMM(NRAO,NCAO,NCMO,1.0d0,                              &
     &       'N','N',IQW ,BUF,NRAO,NCMO,NZW,                            &
     &       'H','N',IQTM2,TM2,LR2,LC2,NZTM2,                           &
     &       FADD,IQAO,FAO,LRAO,LCAO,NZAO)
        ENDIF
      ELSE
        call my_quit('FTRANS: Unknown keyword TYP '//TYP)
      ENDIF
!
      end subroutine

      SUBROUTINE QTRANS90(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,            &
     &                  FAO,LRAO,LCAO,NZAO,IQAO,                        &
     &                  FMO,LRMO,LCMO,NZMO,IQMO,                        &
     &                  TM1,LR1,LC1,NZTM1,IQTM1,                        &
     &                  TM2,LR2,LC2,NZTM2,IQTM2,IPRINT)
!     Like QTRANS, but uses fortran 90 memory allocation internally
      real(8), allocatable :: work(:)
      integer :: NRAO,NCAO,NRMO,NCMO
      integer :: LRAO,LCAO,NZAO
      integer :: LRMO,LCMO,NZMO
      integer :: LR1,LC1,NZTM1
      integer :: LR2,LC2,NZTM2
      integer :: iprint
      CHARACTER TREV*1,TRA*1,TYP*4
      real(8) :: FAO(LRAO,LCAO,*),                                      &
     &          FMO(LRMO,LCMO,*),                                       &
     &          TM1(LR1,LC1,*),                                         &
     &          TM2(LR2,LC2,*)
      real(8) :: fadd
      integer :: lbuf
      integer :: IQAO(*),                                               &
     &          IQMO(*),                                                &
     &          IQTM1(*),                                               &
     &          IQTM2(*)
      integer :: IQW(4)
      integer :: nbuf, nzw
!     find out how much buffer is needed in qtrans
      IF(TYP.EQ.'AOMO') THEN
         IF (NCAO.GE.NRAO) THEN
            CALL IQPACK(IQAO,NZAO,IQTM2,NZTM2,IQW,NZW)
            NBUF = NRAO*NCMO*NZW
         else
            CALL IQPACK(IQAO,NZAO,IQTM1,NZTM1,IQW,NZW)
            NBUF = NRMO*NCAO*NZW
         endif
      else
         IF(NCMO.GE.NRMO) THEN
            CALL IQPACK(IQMO,NZMO,IQTM2,NZTM2,IQW,NZW)
            NBUF = NRMO*NCAO*NZW
         else
            CALL IQPACK(IQMO,NZMO,IQTM1,NZTM1,IQW,NZW)
            NBUF = NRAO*NCMO*NZW
         endif
      endif

      if(nbuf.le.0)then
        call my_quit(' *** error in QTRANS90: buffer space allocation not  &
     & possible for length <= 0.***')
      end if

      allocate(work(nbuf))
      work = 0
      call QTRANS(TYP,TREV,FADD,NRAO,NCAO,NRMO,NCMO,                    &
     &                  FAO,LRAO,LCAO,NZAO,IQAO,                        &
     &                  FMO,LRMO,LCMO,NZMO,IQMO,                        &
     &                  TM1,LR1,LC1,NZTM1,IQTM1,                        &
     &                  TM2,LR2,LC2,NZTM2,IQTM2,                        &
     &                  work,nbuf,IPRINT)
      deallocate(work)

      end subroutine
      SUBROUTINE QGEMM(M,N,K,ALPHA,FA,TA,IQPA,A,LRA,LCA,NZA,            &
     &                             FB,TB,IQPB,B,LRB,LCB,NZB,            &
     &                              BETA,IQPC,C,LRC,LCC,NZC)
!*****************************************************************************
!
!     Performs matrix-matrix operations on general matrices:
!
!     C := alpha * op(A) * op(B) + beta * C
!
!     where alpha and beta are real constants. A,B and C
!     are matrices that may be real(NZ = 1), complex(NZ = 2)
!     or quaternionic (NZ = 4). Calls BLAS routine DGEMM.
!
!     LRX and LCX are leading rows and columns of matrix X.
!     op(A) is an M by K matrix, op(B) is a K by N matrix
!     and C an M by N matrix.
!
!     op(X) is determined by the corresponding character variable
!     FX:
!       FX = 'N'        Normal matrix
!       FX = 'T'        Transpose..
!       FX = 'C'        Complex conjugate
!       FX = 'H'        Hermitian conjugate
!     and TX
!       TX = 'N'        No transform    a+bj
!       TX = 'I'        i-transform   -i(a+bj)i = a - bj
!       TX = 'J'        j-transform   -j(a+bj)j = a*+b*j
!       TX = 'K'        k-transform   -k(a+bj)k = a*-b*j
!
!     AFUL and BFUL indicates on input what components of
!     matrices A and B are non-zero, whereas CFUL on output
!     indicates what components of C are non-zero.
!
!     Written by T.Saue, University of Oslo, Nov. 1994
!
!*****************************************************************************
      use quaternion_algebra
!
!     Global variables
!
      real(8) ::A(LRA*LCA,NZA),B(LRB*LCB,NZB),C(LRC*LCC,NZC)
      integer ::IQPA(NZA),IQPB(NZB),IQPC(NZC), ifa, ifb
      CHARACTER  FA*1,FB*1,TA*1,TB*1,OA*1,OB*1
      real(8) :: faca, facb, beta
      real(8) :: alpha
      integer :: iza, izb, nza, nzb, nzc, izc, ita, itb, iqc, iqa, iqb
      integer :: m, n, k
      integer :: LRA,LCA,                                               &
     &           LRB,LCB,                                               &
     &           LRC,LCC
!
!     Local variables
!
      integer :: IPB(4)
      LOGICAL   LBUF(4)
!
!
!     Initialize LBUF
!     ===============
!
      lbuf = .false.
!
!     Make pointer IPB
!     ================
!
      IPB = 0
      DO IZB = 1,NZB
        IPB(IQPB(IZB)) = IZB
      ENDDO
!
!     MATRIX A
!

!     Determine FX:
!     =============
!     Normal matrix...
      IF    (FA.EQ.'N') THEN
        OA  = 'N'
        IFA = 1
!     Transpose...
      ELSEIF(FA.EQ.'T') THEN
        OA  = 'T'
        IFA = 1
!     Complex conjugate....
      ELSEIF(FA.EQ.'C') THEN
        OA  = 'N'
        IFA = 2
!     Hermitian conjugate...
      ELSEIF(FA.EQ.'H') THEN
        OA  = 'T'
        IFA = 2
      ELSE
        call my_quit('QGEMM:Unknown FA '//FA//' of matrix A!')
      ENDIF
!     Determine TX:
!     =============
!     No transform
      IF    (TA.EQ.'N') THEN
        ITA = 1
!     i-transform
      ELSEIF(TA.EQ.'I') THEN
        ITA = 2
!     j-transform
      ELSEIF(TA.EQ.'J') THEN
        ITA = 3
!     k-transform
      ELSEIF(TA.EQ.'K') THEN
        ITA = 4
      ELSE
        call my_quit('QGEMM:Unknown TA '//TA//' of matrix A!')
      ENDIF
!
!     MATRIX B
!
!     Determine FX:
!     =============
!     Normal matrix...
      IF    (FB.EQ.'N') THEN
        OB  = 'N'
        IFB = 1
!     Transpose...
      ELSEIF(FB.EQ.'T') THEN
        OB  = 'T'
        IFB = 1
!     Complex conjugate....
      ELSEIF(FB.EQ.'C') THEN
        OB  = 'N'
        IFB = 2
!     Hermitian conjugate...
      ELSEIF(FB.EQ.'H') THEN
        OB  = 'T'
        IFB = 2
      ELSE
        call my_quit('QGEMM:Unknown format '//FA//' of matrix B!')
      ENDIF
!     Determine TX:
!     =============
!     No transform
      IF    (TB.EQ.'N') THEN
        ITB = 1
!     i-transform
      ELSEIF(TB.EQ.'I') THEN
        ITB = 2
!     j-transform
      ELSEIF(TB.EQ.'J') THEN
        ITB = 3
!     k-transform
      ELSEIF(TB.EQ.'K') THEN
        ITB = 4
      ELSE
        call my_quit('QGEMM:Unknown TB '//TB//' of matrix B!')
      ENDIF
!
      DO 10 IZC = 1,NZC
        IQC = IQPC(IZC)
        DO 20 IZA = 1,NZA
          IQA = IQPA(IZA)
          IQB = IQMULT(IQA,IQC,1)
          IZB = IPB(IQB)
          IF(IZB.EQ.0) GOTO 20
          FACA = ALPHA*(IQSIGN(IQA,IFA,ITA)*IQSIGN(IQB,IFB,ITB)         &
     &                *IQPHASE(IQA,IQB,1))
          IF(LBUF(IZC)) THEN
            FACB = 1.0d0
          ELSE
            FACB = BETA
          ENDIF
          CALL DGEMM(OA,OB,M,N,K,FACA,A(1,IZA),LRA,B(1,IZB),LRB,        &
     &               FACB,C(1,IZC),LRC)
          LBUF(IZC) = .TRUE.
   20   CONTINUE
   10 CONTINUE
!
!     Check that matrix C is filled
!
      DO IZC = 1,NZC
      IF(.NOT.LBUF(IZC))                                                &
     &    WRITE(*,'(A,I5)') '* WARNING(QGEMM): Not filled: ',IZC
      ENDDO
      end subroutine

      SUBROUTINE IQPACK(IQPA,NZA,IQPB,NZB,IQPC,NZC)
!***********************************************************************
!
!     IQPA and IQPB are pointers to quaternion units 1,i,j and k
!     IQPACK considers what quaternion units span the quaternion product
!     A * B and makes packs these in pointer IQPC as generated.
!
!     Written by T.Saue - June 27 1996
!     Last revision: June 27 1996 - tsaue
!
!***********************************************************************
      use quaternion_algebra
!
      LOGICAL LBUF(4)
      integer :: IQPA(NZA),IQPB(NZB),IQPC(*)
      integer :: iqa, iqb, iqc
      integer :: iza, izb, izc
      integer :: NZA, nzb, nzc
!
!     Initialize LBUF
!     ===============
!
      LBUF(1) = .FALSE.
      LBUF(2) = .FALSE.
      LBUF(3) = .FALSE.
      LBUF(4) = .FALSE.
!
!     Determine quaternion units that contribute to A*B
!     =================================================
!
      NZC = 0
      DO IZA = 1,NZA
        IQA = IQPA(IZA)
        DO IZB = 1,NZB
          IQB = IQPB(IZB)
          IQC = IQMULT(IQA,IQB,1)
          IF(.NOT.LBUF(IQC)) THEN
            NZC = NZC + 1
            IQPC(NZC) = IQC
            LBUF(IQC) = .TRUE.
          ENDIF
        ENDDO
      ENDDO
!
      end subroutine

      SUBROUTINE ADACMO(NB,NZ,CF1,NO1,NP1,NE1,CF2,NO2,NP2,NE2)
!***********************************************************************
!
!     Adapt coefficients CF2 to current dimensions CF1
!     Written by T.Saue
!
!***********************************************************************
      real(8) :: CF1(NB,NO1,NZ),CF2(NB,NO2,NZ)
      integer :: NB,NZ, NO1,NP1, NE1, NO2,NP2,NE2
      integer :: iv1, iv2, iz, nr, nn, nc
!     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        NN = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)
      ENDDO
!     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)
      ENDDO
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        NN  = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
      ENDIF
!
      end subroutine

      SUBROUTINE ADAEIG(EIG1,NO1,NP1,NE1,EIG2,NO2,NP2,NE2)
!***********************************************************************
!
!     Adapt eigenvalues EIG2 to current dimensions EIG1
!     Written by T.Saue
!
!***********************************************************************
      real(8) :: EIG1(NO1),EIG2(NO2)
      integer :: NO1,NP1,NE1, NO2,NP2,NE2
      integer :: iv1, iv2, nr, nc
!     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        CALL DZERO(EIG1(IV1),NR)
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)
!     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        CALL DZERO(EIG1(IV1),NR)
      ENDIF
!
      end subroutine
      SUBROUTINE ADAIBP(IBO1,IBO2,IDIM)
!***********************************************************************
!
!     Pre-routine to adaibo
!     Written by T.Saue Feb 1 2001
!
!***********************************************************************
!
      integer :: IBO1(*),IBO2(*),IDIM(3,2), i2, ifrp
      integer :: ioff, np, ni, no, ne
      I2 = 1
      DO IFRP = 1,nr_fermion_coreps
             if (ifrp == 1) then
                np = nr_mo_gerade_negative_secondary
                ni = nr_mo_gerade_positive_inactive
                no = nr_mo_gerade
                ioff = 0
             else
                np = nr_mo_ungerade_negative_secondary
                ni = nr_mo_ungerade_positive_inactive
                no = nr_mo_ungerade
                ioff = nr_mo_gerade
             end if
             ne = no - np
        CALL ADAIBO(                                                    &
     &         IBO1(ioff      +1),NO,NP,NE,     &
     &         IBO2(I2),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
          I2 = I2 + IDIM(2,IFRP)+IDIM(1,IFRP)
      ENDDO
!
      end subroutine

      SUBROUTINE ADAIBO(IBO1,NO1,NP1,NE1,IBO2,NO2,NP2,NE2)
!***********************************************************************
!
!     Adapt boson irrep array IBO2 to current dimensions IBO1
!     Written by T.Saue
!
!***********************************************************************
      integer :: IBO1(NO1),IBO2(NO2)
      integer :: NO1,NP1,NE1, NO2,NP2,NE2, iv1, iv2, nc, nr
!     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        ibo1(iv1:iv1+nr) = 0
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      ibo1(iv1:iv1+nc) = ibo2(iv2:iv2+nc)
!     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      ibo1(iv1:iv1+nc) = ibo2(iv2:iv2+nc)
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        ibo1(iv1:iv1+nr) = 0
      ENDIF
!
      end subroutine

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GTDOAV(mat_dim,                                        &
     &                  nz,                                             &
     &                  nr_dmat,                                        &
     &                  dmat,                                           &
     &                  ifac)
!***********************************************************************
!
!     Get average density from open shells; only real part !
!     Written by T.Saue March 2004
!
!***********************************************************************
      integer :: mat_dim, nz, nr_dmat
      integer :: ifac, ishell
      real(8) :: DMAT(mat_dim, mat_dim, nz, nr_dmat)
      real(8) :: fac
      if (nr_open_shells > 0) then
        DO ISHELL = 1, nr_open_shells
          FAC = IFAC*df_open(ISHELL)
          CALL DAXPY(mat_dim*mat_dim,FAC,                               &
     &        DMAT(1,1,1,ISHELL+1),1,DMAT,1)
        ENDDO
      end if
!
      end subroutine

end module
      SUBROUTINE my_Q2BPHASE(TYP,IM,JFAC,AMAT)
!***********************************************************************
!
!     Quaternion phase insertion in square matrix
!       TYP = F:  (e_K1)(e_IM)(e_K2*)   as Fock matrix
!       TYP = D:  (e_K1*)(e_IM)(e_K2)   as density matrix
!
!     Written by T.Saue, January 1995
!     Last revision: June 25 1996 - tsaue
!
!***********************************************************************
      use interface_ao
      use quaternion_algebra
      implicit none
      CHARACTER TYP*1
      real(8) :: AMAT(*)
      integer :: ic1, ic2, irp1, irp2, isy1, isy2
      integer :: k1, k2, kfac
      integer :: nb1, nb2, ioff, icol
      integer :: im, jfac

      integer, parameter :: iqph(4, 2) = reshape((/1, 1, 1, 1, 1, -1, -1, -1/), (/4, 2/))

      DO 10 IRP1 = 0,nr_boson_ireps-1
        ISY1 = IRP1 + 1
        DO 20 IC1 = 1,2
          NB1 = ao_nr(IRP1,IC1)
          IF(NB1.EQ.0) GOTO 20
          K1  = ao_q(IRP1,IC1)
          DO 30 IRP2 = 0,nr_boson_ireps-1
            ISY2 = IRP2 + 1
            DO 40 IC2 = 1,2
              NB2 = ao_nr(IRP2,IC2)
              IF(NB2.EQ.0) GOTO 40
              K2   = ao_q(IRP2,IC2)
              IF    (TYP.EQ.'F') THEN
                KFAC = JFAC*IQPHASE(K1,K2,IM)*IQPH(K1,1)*IQPH(K2,2)
              ELSEIF(TYP.EQ.'D') THEN
                KFAC = JFAC*IQPHASE(K1,K2,IQMULT(K1,K2,IM))*IQPH(K1,2)*IQPH(K2,1)
              ENDIF
              IF(KFAC.EQ.-1) THEN

                IOFF = ao_off(irp2, ic2)*nr_ao + ao_off(irp1, ic1)
                DO 50 ICOL = 1,NB2
                  CALL DSCAL(NB1,-1.0d0,AMAT(IOFF+1),1)
                  IOFF = IOFF + nr_ao
   50           CONTINUE
              ENDIF
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
!
      end subroutine

      SUBROUTINE my_QUIT(TEXT)
      implicit none
      CHARACTER TEXT*(*)
      stop
      end subroutine
      LOGICAL FUNCTION mylbit(I,N)
!       Based on an analogous DISCO routine(Jan Almloef)
      integer :: i, n
      mylbit = .FALSE.
      IF (iand(ishft(I,1-N),1).EQ.1) mylbit = .TRUE.
      end function


      subroutine my_fortwrt(str, slength)
      implicit none
#ifdef INT_STAR8
      integer(8) :: slength
#else
      integer    :: slength
#endif /* ifdef INT_STAR8 */
      character  :: str(slength)

      write(*, *) str

      end subroutine
