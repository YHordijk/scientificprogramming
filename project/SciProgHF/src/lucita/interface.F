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

***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
**********************************************************************
*                                                                    *
*                                                                    *
*          MOLCAS - LUCIA interface                                  *
*          DIRAC  - LUCIA interface                                  *
*                                                                    *
*          written by Timo Fleig, Feb.-Apr. 1999                     *
*          written by Timo Fleig, Dec       2001                     *
*                                                                    *
*                                                                    *
*                                                                    *
**********************************************************************
*
      subroutine dirluc
*
      implicit real*8 (A-H,O-Z)
*
#include "clunit.inc"
*
*  Say hello:
      call hello_dirluc
*
*  Get it going:
      call dirluct
*
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*********************************************************************
*  Note concerning EXPERT mode:
*   If you want to use it, the first 2 lines of your input for
*   MOLUC-LUCIA (LUCITA) should be
*       &LUCITA &END
*      EXPERT
*      * &LUCIA &END
*       ...    (first line of regular input for LUCIA)
*       ...    (regular input for LUCIA)
*       ...
*       END OF INPUT
*   And ONE blank line after 'END OF INPUT'. This is essential,
*   otherwise the program will crash.
*
*  Feb-Mar 1999,   Timo Fleig
*  Feb     2000,   Timo Fleig
*
*  Modified for DIRAC - LUCITA information transfer.
*  Dec-Mar 2001/02,   Timo Fleig
*********************************************************************
*
      subroutine dirluct
*
      implicit real*8 (A-H,O-Z)
*
#include "luctinp.h"
*
#include "clunit.inc"
#include "mxpdim.inc"
#include "luci.inc"
#include "units.inc"
#include "infpar.h"
#include "parluci.h"
*
      parameter (NCMD = 21)
      parameter (MXPLNC = 72)
C     dimension IMOKW(NTABLE)
      character*1 CRDHEL
      character*3 SCR3,YESSTATEMENT
      character*10 PERFORMCI
      character*4 COMMAND,CMD(NCMD)
      character*6 SCR6,MOLUCENV
      character*72 LINE,CARD(MXPNGAS),CARD2(MXPNGAS,MXPICI),CARD3(5)
*
      common/MLCENV/MOLUCENV
      real(8) :: energy_convergence
*
*
* Set logical units
      call setunits
*
*----------------------------------------------------------------------*
*  Insert defaults and print (error) messages where appropriate
*----------------------------------------------------------------------*
*
*  Code for keyword status vector:
*
*     IMOKW(N) = 0           keyword not provided in input,
*                            default or error
*     IMOKW(N) = 1           keyword given in input
*     IMOKW(N) = 2           keyword has been set to default value
*     IMOKW(N) = 3           keyword ignored, default is used
*
*
*
*  1: Default Title
*
      if (IMOKW(1).eq.0) then
         NTIT = 1
         TITLE(1) = ' Running LUCIA under DIRAC. No title supplied. '
         IMOKW(1) = 2
      else if (IMOKW(1).eq.1) then
         NTIT = 1
         TITLE(1) = TITLUC
      end if
*
*  2: Type of initial wave function (no default)
*
      if (IMOKW(2).eq.0) then
         write(6,*)
     & ' Keyword for type of initial wave function missing. '
         write(6,*) ' This keyword is mandatory. '
         Call Abend2('Input error')
      else
         call lftpos(WAFFCD,MXPLNC)
         call uppcas(WAFFCD,MXPLNC)
         if (WAFFCD.eq.'DHFSCF') WAFFCT = 'HF_SCF'
         if (WAFFCD.eq.'OSHSCF') WAFFCT = 'RASSCF'
         if (WAFFCT.ne.'HF_SCF'.and.WAFFCT.ne.'RASSCF') then
            write(6,*)
     &  ' Type of initial wave function not implemented. '
            write(6,*) ' You have chosen: ',WAFFCD
            write(6,*) ' Allowed types: '
            write(6,*) ' DHFSCF, OSHSCF '
            Call Abend2('Input error')
         end if
      end if
*
*  3: CI type for LUCIA (no default)
*
      if (IMOKW(3).eq.0) then
         write(6,*) ' Keyword for type of CI calculation missing. '
         write(6,*) ' This keyword is mandatory. '
         Call Abend2('Input error')
      else
         call lftpos(CALCTP,MXPLNC)
         call uppcas(CALCTP,MXPLNC)
         TYPE = CALCTP
         if (TYPE.ne.'FCI   '.and.TYPE.ne.'SDCI  '.and.
     &       TYPE.ne.'GASCI '.and.TYPE.ne.'RASCI '.and.
     &       TYPE.ne.'SDTQ  ') then
            write(6,*) ' Type of CI calculation not specified.'
            write(6,*) ' You have chosen: ',TYPE
            write(6,*) ' Allowed types: '
            write(6,*) ' FCI, SDCI, GASCI, RASCI, SDTQ '
            Call Abend2('Input error')
         end if
      end if
*
*  4: Number of roots to be obtained
*
      if (IMOKW(4).eq.0) then
         NROOT = 1
         IMOKW(4) = 2
      else if (IMOKW(4).eq.1) then
         NROOT = NROOTD
      end if
*
*  5: State symmetry in point group
*
      if (IMOKW(5).eq.0) then
         IRREP = 1
         IMOKW(5) = 2
      else if (IMOKW(5).eq.1) then
         IRREP = ISSYMD
      end if
*
*  6: Active electrons
*
      if (IMOKW(6).eq.0) then
         if (TYPE.eq.'GASCI '.or.TYPE.eq.'RASCI ') then
            write(6,*) 'Number of active electrons NACTEL'
            write(6,*) 'has to be specified in this type of'
            write(6,*) 'calculation. Quitting.'
            Call Abend()
         else
            NACTEL = -1
            IMOKW(6) = 2
         end if
      else if (IMOKW(6).eq.1) then
         NACTEL = NACTED
      end if
*
*  7: Spin multiplicity (no default)
*
      if (IMOKW(7).eq.1) then
         MULTIP = IMULTD
         if (IMOKW(6).eq.1) then
            call evenodd(IEONAC,NACTEL)
            call evenodd(IEOMUL,MULTIP)
            if (IEONAC.eq.2.and.IEOMUL.eq.2) then
               write(6,*) 'Illegal spin multiplicity given.'
               write(6,*) 'Read a book about fermions.'
               Call Abend2('quitting')
            else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
               write(6,*) 'Illegal spin multiplicity given.'
               write(6,*) 'Read a book about fermions.'
               Call Abend2('quitting')
            else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
               if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
                  write(6,*) 'Illegal spin multiplicity given.'
                  write(6,*) 'Compare with number of active'
                  write(6,*) 'electrons.'
                  Call Abend2('quitting')
               end if
            else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
               if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
                  write(6,*) 'Illegal spin multiplicity given.'
                  write(6,*) 'Compare with number of active'
                  write(6,*) 'electrons.'
                  Call Abend2('quitting')
               end if
            end if
         end if
      else if (IMOKW(7).eq.0) then
         write(6,*) 'Spin multiplicity is a MANDATORY keyword.'
         write(6,*) 'Specify and restart.'
         Call Abend2('Quitting.')
      end if
*
*  8: Global print parameter
*
      if (IMOKW(8).eq.0) then
         PRILUC = 'NON'
         IMOKW(8) = 2
      else if (IMOKW(8).eq.1) then
         if (IPRNGD.eq.0) then
            PRILUC = 'NON'
         else if (IPRNGD.eq.1) then
            PRILUC = 'LOW'
         else if (IPRNGD.eq.2) then
            PRILUC = 'MED'
         else if (IPRNGD.eq.3) then
            PRILUC = 'HIG'
         else if (IPRNGD.eq.4) then
            PRILUC = 'VHI'
         else
            write(6,*) 'Invalid print flag. Check PRINTG.'
            write(6,*) 'IPRNGD = ',IPRNGD
            Call Abend2('quitting')
         end if
      end if
*
*  9: Local print parameter
*
      if (IMOKW(9).eq.0) then
         IPLOCAL = 0
         IMOKW(9) = 2
      else if (IMOKW(9).eq.1) then
         IPLOCAL = IPRNLD
      end if
*
* 10: Approximate size of calculation
*
      if (IMOKW(10).eq.0) then
         CALSIZ = 'NOR'
         IMOKW(10) = 2
      else
         call lftpos(SZCALD,MXPLNC)
         call uppcas(SZCALD,MXPLNC)
         CALSIZ = SZCALD(1:3)
         if (CALSIZ.ne.'NOR'.and.CALSIZ.ne.'LAR'.and.
     &       CALSIZ.ne.'HUG') then
           write(6,'(A,A3)') 'SIZE specified to ',CALSIZ
           write(6,*) 'This is an unknown type.'
           Call Abend2('Quitting.')
         end if
      end if
*
* 11: Orbital distribution in inactive space (no defaults if RASCI)
*
      if (IMOKW(11).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'Number of inactive orbitals per sym'
            write(6,*) 'has to be specified. This is mandatory'
            write(6,*) 'in a RASCI calculation.'
            Call Abend2('Quitting.')
         else
            INACOB = 0
            IMOKW(11) = 2
         end if
      else if (IMOKW(11).eq.1) then
         if (TYPE.eq.'GASCI ') then
            write(6,*) 'INACtive orbitals not allowed in GASCI.'
            write(6,*) 'Use GAS 1 for these orbitals and restart.'
            Call Abend2('Quitting.')
         end if
         CARD3(2) = CRDINA
         INACOB = 1
      end if
*
* 12: Orbital distribution in GAS spaces (no defaults if GASCI)
*
      if (IMOKW(12).eq.0) then
         if (TYPE.eq.'GASCI ') then
            write(6,*) 'GASCI type requires GASS to be specified.'
            write(6,*) 'Else, I do not know what to do.'
            Call Abend2('Quitting.')
         end if
         IMOKW(12) = 2
      end if
      if (IMOKW(12).eq.1) then
         if (TYPE.ne.'GASCI ') then
           write(6,*) 'GASS specified and not doing GASCI.'
           write(6,*) 'This input will be ignored.'
           write(6,*) 'Program will continue with ',TYPE,'.'
           IMOKW(12) = 3
         else
           NGAS = INGASD
           do IGS = 1,NGAS,1
             CARD(IGS) = CRDGAS(IGS)
           end do
         end if
      end if
*
* 13: Cumulative min. and max. numbers of electrons in GAS spaces
*
      if (IMOKW(13).eq.0) then
         if (TYPE.eq.'GASCI ') then
            write(6,*) 'GSSP has to be specified in GASCI calculation.'
            write(6,*) 'No defaults in this type of CI.'
            Call Abend2('Quitting.')
         end if
         IMOKW(13) = 2
      end if
      if (IMOKW(13).eq.1) then
         if (TYPE.ne.'GASCI ') then
           write(6,*) 'GSSP specified and not doing GASCI.'
           write(6,*) 'This input will be ignored.'
           write(6,*) 'Program will continue with ',TYPE,'.'
           write(6,*) 'Defaults will be set.'
           IMOKW(13) = 3
         else
           NCISPC = NSEQCD
           do ISPC = 1,NCISPC,1
             do IGS = 1,NGAS,1
               CARD2(IGS,ISPC) = CRDGOC(IGS)
             end do
           end do
         end if
      end if
*
* 14: FRMO (frozen orbitals from DIRAC) specification; no default
*     can be used with any CI or WF type
*
      if (IMOKW(14).eq.0) then
         CRDHEL(1:1) = 'N'
         CARD3(1) = CRDHEL
         IMOKW(14) = 2
      else if (IMOKW(14).eq.1) then
         CARD3(1) = CRDFRO
      end if
*
* 15: RAS1 specification; no default
*
      if (IMOKW(15).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'RAS1 keyword missing. Add and restart.'
            Call Abend2('Quitting.')
         end if
         IMOKW(15) = 2
      else if (IMOKW(15).eq.1) then
         if (TYPE.ne.'RASCI ') then
            write(6,*) 'Keyword RAS1 not compatible with TYPE ',TYPE
            write(6,*) 'I will ignore your this input.'
            write(6,*) 'Program will continue.'
            IMOKW(15) = 3
         end if
         CARD3(3) = CRDRS1
         MXHOL1 = MXHL1D
      end if
*
* 16: RAS2 specification; no default
*
      if (IMOKW(16).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'RAS2 keyword missing. Add and restart.'
            Call Abend2('Quitting.')
         end if
         IMOKW(16) = 2
      else if (IMOKW(16).eq.1) then
         if (TYPE.ne.'RASCI ') then
            write(6,*) 'Keyword RAS2 not compatible with TYPE ',TYPE
            write(6,*) 'I will ignore your this input.'
            write(6,*) 'Program will continue.'
            IMOKW(16) = 3
         end if
         CARD3(4) = CRDRS2
      end if
*
* 17: RAS3 specification; no default
*
      if (IMOKW(17).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'RAS3 keyword missing. Add and restart.'
            Call Abend2('Quitting.')
         end if
         IMOKW(17) = 2
      else if (IMOKW(17).eq.1) then
         if (TYPE.ne.'RASCI ') then
            write(6,*) 'Keyword RAS3 not compatible with TYPE ',TYPE
            write(6,*) 'I will ignore your this input.'
            write(6,*) 'Program will continue.'
            IMOKW(17) = 3
         end if
         CARD3(5) = CRDRS3
         MXELR3 = MXEL3D
      end if
*
* 18: Default setting for DENSity matrices
*
      if (IMOKW(18).eq.0) then
         IDEMOL = 0
         IMOKW(18) = 2
      else if (IMOKW(18).eq.1) then
         IDEMOL = IDENSD
      end if
*
* 19: Restart option from CI vector(s)
*
      if (IMOKW(19).eq.0) then
         IRSTRT = 0
         IMOKW(19) = 2
      else if (IMOKW(19).eq.1) then
         IRSTRT = IRSTLT
      end if
*
* 21: Type of parallel calculation
*
      if (IMOKW(21).eq.0) then
#if defined (VAR_MPI2)
         IODENSEPAR = 2
         IRUNPA     = 1
         IPARMODEL  = 2
#else
         IODENSEPAR = 0
         IRUNPA     = 0
         IPARMODEL  = 0
#endif
         IMOKW(21) = 2
      end if
*
* 22: Type of I/O model in parallel calculation
*
      if (IMOKW(22).eq.0) then
#if defined (VAR_MPI2)
         IIOMOD = 1
         YESSTATEMENT = 'Yes'
#else
         IIOMOD = 0
         YESSTATEMENT = 'No'
#endif
         IMOKW(22) = 2
      end if
*
* 23: Number of maximum CI Iterations
*
      if (IMOKW(23).eq.0) then
         IMAXCIITER = 100
         IMOKW(23) = 2
      else if (IMOKW(23).eq.1) then
         if(ICIMAXITER.lt.0) then
           IMAXCIITER = -1
         else
           IMAXCIITER = ICIMAXITER
         end if
      end if
*
* 24: Largest allowed batchsize (if larger than biggest block!) 
*
      if (IMOKW(24).eq.0) then
         IMAXBLKSIZE = 100 000 000
         IMAXLBLKSZ = IMAXBLKSIZE
         IMOKW(24) = 2
      else if (IMOKW(24).eq.1) then
         IMAXLBLKSZ = MAX(100 000,IMAXBLKSIZE)
C        .. we do not accept user input of less than 100 000 !
      end if

!     25: analyze wave function     
      if (IMOKW(25).eq.1) PRILUC = 'LOW'
*
* 26: Distribution routine to use (only relevant in parallel case) 
*
      IDISTROUTE = 0
      if (IMOKW(26).eq.0) then
        IDISTROUTE = 2
        IMOKW(26)  = 2
      else if (IMOKW(26).eq.1) then
        if( I_USE_DIST_ROUTE .eq. 1 ) then
          IDISTROUTE = I_USE_DIST_ROUTE
        else if ( I_USE_DIST_ROUTE .eq. 2 ) then
          IDISTROUTE = I_USE_DIST_ROUTE
        else
          write(6,*) 'Value for keyword DISTRT incorrect ',
     &               I_USE_DIST_ROUTE
          write(6,*) 'I will ignore your input.'
          write(6,*) 'Program will continue.'
          IDISTROUTE = 2   
          IMOKW(26) = 2
        endif
      end if
      if (IMOKW(27).eq.1) then
        energy_convergence = my_convergence
      else
        energy_convergence = 1.0D-8
      end if
*
*
* 28: Truncate residual vectors before creating new trial vector?
*     (14-jun-07, hjaaj)
*
      if (IMOKW(28).eq.1) then
        trunc_fac = ctrunc_fac
      else
      ! default: remove only very small numbers, less than round-off
      ! factor
        trunc_fac = 1.0D-10
      end if
*
*
* 29: memory multiplier ...
*    
*
      if (IMOKW(29).eq.1) then
         ISMEMFAC  = in_memfac 
      else
         ISMEMFAC  = 3
      end if

      if(irunpa.eq.1) then
        performci = 'parallel'
      else if(irunpa.eq.0) then
        performci = 'sequential'
      endif
*
*
*-----------------------------------------------------------------*
*  Print section
*-----------------------------------------------------------------*
*
*  1: Title
*
      write(LUOUT,*)
      write(LUOUT,'(1X,80A1)') ('*',I=1,80)
      write(LUOUT,'(1X,80A1)') '*',(' ',I=1,78),'*'
      write(LUOUT,'(1X,37A1,A6,37A1)')
     &        '*',(' ',I=1,36),'Title:',(' ',I=1,36),'*'
      do I=1,NTIT
         write(LUOUT,'(1X,4A1,A72,4A1)')
     &        '*',(' ',J=1,3),Title(I),(' ',J=1,3),'*'
      end do
      write(LUOUT,'(1X,80A1)') '*',(' ',I=1,78),'*'
      write(LUOUT,'(1X,80A1)') ('*',I=1,80)
      write(LUOUT,*)
*
*  2: Initial wave function
*
      write(LUOUT,'(A42,A6)')
     &  ' Orbitals as initial wave function .... ',WAFFCT
      write(LUOUT,*)
*
*  3: Type of CI calculation
*
      write(LUOUT,'(A42,A6)')
     &  ' Type of calculation .................. ',TYPE
      write(LUOUT,*)
*
*  4: Number of roots to be treated
*
      write(LUOUT,'(A42,I3)')
     &  ' Number of roots to be obtained ....... ',NROOT
      write(LUOUT,*)
*
*  5: State symmetry
*
      write(LUOUT,'(A42,I3)')
     &  ' Calculation carried out in irrep ..... ',IRREP
      write(LUOUT,*)
*
*  6: Number of active electrons
*
      if (IMOKW(6).eq.1) then
        write(LUOUT,'(A42,I3)')
     &  ' Number of active electrons ........... ',NACTEL
        write(LUOUT,*)
      end if
*
*  7: Spin multiplicity
*
      write(LUOUT,'(A42,I3)')
     &  ' Spin multiplicity .................... ',MULTIP
      write(LUOUT,*)
*
*  8: LUCIA global print parameter
*
      write(LUOUT,'(A42,A3)')
     &  ' Global print level is ................ ',PRILUC
      write(LUOUT,*)
*
*  9: MOLUC local print parameter
*
      write(LUOUT,'(A42,1I3)')
     &  ' Local print level is ................. ',IPLOCAL
      write(LUOUT,*)
*
* 10: SIZE of CI calculation
*
      write(LUOUT,'(A42,A3)')
     &  ' Approximate size of CI calculation ... ',CALSIZ
      write(LUOUT,*)
*
* 10b: runtype of CI calculation
*
      write(LUOUT,'(A42,A8)')
     &  ' Running the CI calculation ........... ',performci
      write(LUOUT,*)
*
*
* 10c: runtype of parallel calculation
*
      if(irunpa.gt.0)then
        write(LUOUT,'(A42,I3)')
     &  ' Applying parallel scheme ............. ',iparmodel
        write(LUOUT,*)
      endif
*
*
* 10d: file version to use in parallel calculation
*
      write(LUOUT,'(A42,A3)')
     &  ' Using MPI-FILE I/O ................... ',YESSTATEMENT
      write(LUOUT,*)
*
*
* 10e: file system in parallel calculation
*
      if( irunpa .gt. 0 ) then
        write(LUOUT,'(A42,I1)')
     &  ' Parallel distribution routine ........ ',IDISTROUTE
      end if
*
* 10f: truncation factor
*
      write(LUOUT,'(/A42,1P,D10.2)')
     &  ' Truncation Factor ..................... ',ctrunc_fac
      write(LUOUT,*)
                   
* 11: INAC specification of inactive orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 12: GASS specification of orbitals per GAS per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 13: GSSP specification of electrons in GAS Spaces
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 14: FRMO specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 15: RAS1 specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 16: RAS2 specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 17: RAS3 specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 18: DENSity matrix level
*
      if (IDEMOL.ge.1) then
        write(LUOUT,'(A42,1I3)')
     &    ' Calculating density matrices at level  ',IDEMOL
        write(LUOUT,*)
      end if
*
* 19: Restart option activated
*
      if (IRSTRT.eq.1) then
        write(LUOUT,'(A42,A3)')
     &    ' Restarted calculation from file ...... ','YES'
        write(LUOUT,*)
      end if
*
*
*
*-----------------------------------------------------------------*
*  Call conversion routine for creating LUCIA input
*-----------------------------------------------------------------*
      call mol2luc(NTIT,NGAS,NCISPC,IPLOCAL,CARD,CARD2,CARD3,
     &             energy_convergence)
*
6666  continue
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Distribute externally transformed (MOLTRA) dirac one-electron
*  integrals to double group lists
*
*----------------------------------------------
*  Timo Fleig, Oct. 2000
*
*   modified for scalar relativistic runs with LUCITA
*     Nov. 2000
*   extended for use of groups C2, Cs, C2h, D2h
*     Feb. 2001
*   generalized to all kinds of group generators and
*     arbitrary order.
*     28.Feb.2001
*----------------------------------------------
*
      subroutine dist_1eints_dirac2(HELPMAT,ONEINT,
     &                             OTRINT,IREO,
     &                             NTRELM,ITROFF,
     &                             ITOTDIM2,ITOTDIM,
     &                             IDIM,NIRREP,IRR,
     &                             NTOORB,NTEST)
*
      implicit real*8 (A-H,O-Z)
*
      dimension ONEINT(ITOTDIM2,ITOTDIM2),HELPMAT(IDIM,IDIM)
      dimension OTRINT(*)
      dimension IREO(IDIM),NTOORB(NIRREP)
*
      NTESTL = 00
      NTEST = max(NTESTL,NTEST)
*
*  Zero out help matrix
      ZERO = 0.D0
      do I=1,IDIM,1
        call setvec(HELPMAT(1,I),ZERO,IDIM)
      end do
*
      do IR=1,IDIM,1
        do IL=1,IDIM,1
          HELPMAT(IL,IR) = ONEINT(IREO(IL),IREO(IR))
        end do
      end do
*
*  Check picked diagonal partition
      if (NTEST.ge.3) then
        write(6,'(A,I4)') 'Checking input matrix for fold, Sym',IRR
        call wrtmat(HELPMAT,
     &              NTOORB(IRR),NTOORB(IRR),
     &              NTOORB(IRR),NTOORB(IRR))
      end if
*
*  Fold to lower triangle list
      call fold(OTRINT(1+ITROFF),HELPMAT,NTOORB(IRR),1,NTRELM)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Create symmetry reduced 2-el. integral lists from
* read full Dirac integral lists
* Distribute integrals to double group lists in the
* same fashion as scalar-relativistic integrals (spinfree
* formalism)
*
*  Input:  RKL
*  Output: TWELINT
*
      subroutine dist_2eints_dirac(INVREO,ITOTDIM,
     &                             RKL,KR,IOFF12,NUMINT,
     &                             INDK,INDL,TWELINT,IGLSUM,
     &                             SYDI,NIRR_PN,IPRINT)
*
      implicit real*8 (A-H,O-Z)
*
#include "multd2h.inc"
#include "parluci.h"
*
*--------------------------------------
*  Timo Fleig, October 2000
*
*   modified for scalar relativistic runs with LUCITA
*   Nov. 2000
*
*   New DIRAC integrals have full 8-fold permutation
*   symmetry. No introduction of quaternion factors in
*   spinfree run (Luuk Visscher, Feb. 2001)
*
*   implemented quaternion groups Ci and C1
*   May 2001
*--------------------------------------
*
      integer SYDI
*
      dimension INVREO(ITOTDIM/2),KR(-ITOTDIM/2:ITOTDIM/2),
     &          IOFF12(ITOTDIM,ITOTDIM),
     &          NUMINT(ITOTDIM,ITOTDIM),
     &          INDK(ITOTDIM**4),INDL(ITOTDIM**4),SYDI(NIRR_PN)
      dimension RKL(ITOTDIM**4),TWELINT(IGLSUM)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
*  Create reduced integral list as input for GETINCN
*  according to:
*   1) Kramers symmetry
*   2) Particle symmetry
*   3) Real/imaginary symmetry
*
!     initialize TWELINT array
      call dzero(twelint,iglsum)
!
      ZERO = 0.D0
*  Proceed according to spinfree formalism, i.e. only (IJ|KL) list
      ICT = 1
      IIOFF = 1
*  loop over non-redundant output (boson) symmetry blocks
      do ISM=1,NIRR_PN,1
        NORBI = SYDI(ISM)
        ISYM = ISM
        IJOFF = 1
        do JSM=1,ISM,1
          NORBJ = SYDI(JSM)
          JSYM = JSM
          IKOFF = 1
          do KSM=1,ISM,1
            NORBK = SYDI(KSM)
            KSYM = KSM
            LLSM = KSM
            if (KSM.eq.ISM) LLSM = JSM
            ILOFF = 1
            do LSM=1,LLSM,1
              NORBL = SYDI(LSM)
              LSYM = LSM
              IJSYM = MULTD2H(ISYM,JSYM)
              KLSYM = MULTD2H(KSYM,LSYM)
              IJKLSYM = MULTD2H(IJSYM,KLSYM)
              if (IJKLSYM.eq.1) then
                if (NTEST.ge.5) then
                  write(6,*)
                  write(6,*) 'checksym'
                  write(6,'(A,4I4)') 'ISYM,JSYM,KSYM,LSYM ',
     &                                ISYM,JSYM,KSYM,LSYM
                  write(6,*) 'check current offsets'
                  write(6,'(A,4I4)') 'IIOFF,IJOFF,IKOFF,ILOFF ',
     &                                IIOFF,IJOFF,IKOFF,ILOFF
                end if
                CALL FLSHFO(6)
*
*  Index restrictions
*   Real/imaginary symmetry
                if (ISYM.eq.JSYM) then
                  IJ_RES = 1
                else
                  IJ_RES = 0
                end if
                if (KSYM.eq.LSYM) then
                  KL_RES = 1
                else
                KL_RES = 0
                end if
*
*   Particle symmetry
                if (ISYM.eq.KSYM.and.JSYM.eq.LSYM) then
                  IJKL_RES = 1
                else
                  IJKL_RES = 0
                end if
*
* K
*  L
*   I
*    J
                if (NTEST.ge.5) 
     &            write(6,*) ' loop over non-redundant output indices'
*  loop over non-redundant output indices including symmetry offsets
                do K=IKOFF,IKOFF+NORBK-1,1
                  if (KL_RES.eq.1) then
                    L_LAST = K
                  else
                    L_LAST = ILOFF + NORBL -1
                  end if
                  do L=ILOFF,L_LAST,1
                    if (IJKL_RES.eq.1) then
                      I_FIRST = K
                    else
                      I_FIRST = IIOFF
                    end if
                    do I=I_FIRST,IIOFF+NORBI-1,1
                      if (IJ_RES.eq.1) then
                        J_LAST = I
                      else
                        J_LAST = IJOFF + NORBJ -1
                      end if
                      if (IJKL_RES.eq.1.and.I.eq.K) then
                        J_FIRST = L
                      else
                        J_FIRST = IJOFF
                      end if
                      do J=J_FIRST,J_LAST,1
                        if (NTEST.ge.5) then
                          write(6,'(/2x,a)') 
     &        '-------------------------------------------------------'
                          write(6,'(2x,a,4i4)') 'next indices pairs:'//
     &                      ' I, J, K, L = ',I,J,K,L
                          write(6,'(2x,a/)') 
     &        '-------------------------------------------------------'
                        end if
*
*  Translate to input indices from Dirac lists
                        IIN = INVREO(I)
                        JIN = INVREO(J)
                        KIN = INVREO(K)
                        LIN = INVREO(L)
                        if (NTEST.ge.5) then
                          write(6,*) 'Reordered indices:'
                          write(6,'(A,4I6/)')
     &                      'xIN',IIN,JIN,KIN,LIN
                        end if
*
                        if (NTEST.ge.5) then
                          write(6,*) 'Kramers indices:'
                          write(6,'(A,4I6/)')
     &                      'KR(xIN) ',KR(IIN),KR(JIN),KR(KIN),KR(LIN)
                        end if
*
*  Find non-vanishing Kramers block
                        IBLLEN = NUMINT(KR(IIN),KR(JIN))
                        if (JIN.gt.IIN) then
*  Use complex conjugation symmetry. Real integrals.
                          IIN_SAVE = IIN
                          JIN_SAVE = JIN
                          IIN = JIN_SAVE
                          JIN = IIN_SAVE
                          KIN_SAVE = KIN
                          LIN_SAVE = LIN
                          KIN = LIN_SAVE
                          LIN = KIN_SAVE
                          IBLLEN = NUMINT(KR(IIN),KR(JIN))
                        end if
*
                        if (IBLLEN.eq.0) then
                          write(6,*) 'Dist_2eints_dirac: '
                          write(6,*) 'Seeking indices ',IIN,JIN
                          write(6,'(A,I4,A,I4,A,2I4)')
     &                    'KR(',IIN,'), KR(',JIN,') = ',KR(IIN),KR(JIN)
                          write(6,*) 'Vanishing Kramers block.'
                          call abend2('Stop.')
                        end if
*
                        if (NTEST.ge.5) then
                          write(6,*)
     &                    'Kramers indices after symm. reordering:'
                          write(6,'(A,4I6)')
     &                    'KR(xIN) ',KR(IIN),KR(JIN),KR(KIN),KR(LIN)
                        end if
*
*  Fetch integral
                        IBLKOFF = IOFF12(KR(IIN),KR(JIN))
                        do INRUN=1,IBLLEN,1
                          if (KR(INDK(INRUN+IBLKOFF)).eq.KR(KIN).and.
     &                        KR(INDL(INRUN+IBLKOFF)).eq.KR(LIN))
     &                    then
                            TWELINT(ICT) = RKL(INRUN+IBLKOFF)
                            if (NTEST.ge.10) then
                              write(6,'(2x,a/,4i6,1F18.12/)')
     &                          ' integral fetched: integral number,'//
     &                          ' (IBLKOFF + INRUN) = offset on input'//
     &                          ' block, integral value:',ICT,
     &                            IBLKOFF,INRUN,IBLKOFF+INRUN,
     &                            TWELINT(ICT)
                            end if
                            ICT = ICT + 1
                            goto 10
                          end if
                        end do
*
*  Integral has not been found on existing symmetry (IJ| block
*  of input integrals. It has been deleted in DIRAC due to small
*  value. Insert ZERO in its place.
                        TWELINT(ICT) = ZERO
                        if (NTEST.ge.10) then
                          write(6,'(10X,A,I9,A,1F18.12)')
     &                        'Integral ',ICT,' set to ',ZERO
                        end if
                        ICT = ICT + 1
   10                   continue
                        CALL FLSHFO(6)
*                       ^ Generate next output index quadruple
                      end do
                    end do
                  end do
                end do
*               ^ End loops over 4 indices
              end if
*             ^ End if non-vanishing symmetry block
              ILOFF = ILOFF + NORBL
            end do
            IKOFF = IKOFF + NORBK
          end do
          IJOFF = IJOFF + NORBJ
        end do
        IIOFF = IIOFF + NORBI
      end do
*     ^ End loops over symmetries of 4 indices
*
      if (NTEST.ge.3) then
        write(6,'(//a/,a,i8/)') ' (dist_2eints_dirac): symmetry'//
     &                          ' reduced 2-el. integrals:',
     &                          ' total # =',IGLSUM
        write(6,'(2x,4F18.12)') (TWELINT(I), I=1,IGLSUM)
      end if
*
      if (myproc.eq.master) then
        write(LUWRT,'(/A)') 
     &       ' Complete real list of 2-el. ints processed. '
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Driver routine for distribution of DIRAC integrals to
* double group GAS lists
*  Creates some generally required information
*  Loops over spinor lists and real/imaginary parts
*
*----------------------------------------------
*  Timo Fleig, Oct. 2000
*
*   modified for scalar relativistic runs with LUCITA
*     Nov. 2000
*   extended for use of groups C2, Cs, C2h, D2h
*     Feb. 2001
*   generalized version independent of the choice and
*     order of symmetry generators
*     28.Feb. 2001
*----------------------------------------------
*
      subroutine dist_ints_dirac2(HELPMAT,ONEINTR,
     &                           OTRINT,
     &                           RKLR,KR,IOFF12,
     &                           NUMINT,INDK,INDL,
     &                           TWELINT,
     &                           IBOSYM,IREO,
     &                           INVREO,NTOOBS,NSTR,
     &                           ITOTDIM,ITOTDIM2,
     &                           INVERSM,
     &                           IHAM12,
     &                           ISPINFREE,IMXDIM,NIRREP,IPRINT)
*
      implicit real*8 (A-H,O-Z)
*
#include"dcbbas.h"
#include"dcborb.h"
*
      dimension OTRINT(*),ONEINTR(ITOTDIM2,ITOTDIM2)
      dimension IREO(IMXDIM,NIRREP),INVREO(ITOTDIM),IBOSYM(ITOTDIM2)
      dimension NTOOBS(NIRREP),NSTR(INVERSM+1)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      ICT = 1
      JCT = 1
*
* NSTR contains only active set of Kramers pairs.
*
!     initialize
      call izero(ireo,  imxdim*nirrep)
      call izero(invreo,itotdim)
!
* Determine reordering arrays for 1- and 2-electron integrals
      do IRR=1,NIRREP,1
        IOFFBS = 0
        do ISP=1,ITOTDIM,1
          if (IBOSYM(ISP+IOFFBS)+1.eq.IRR) then
            INVREO(ICT) = ISP
            ICT = ICT + 1
            IREO(JCT,IRR) = ISP + IOFFBS
            JCT = JCT + 1
            if (NTEST.ge.5) then
         write(6,'(A,I4,A,I4)')  'ISP,INVREO(',ICT-1,') = ',ISP
         write(6,'(A,2I4,A,I4)') 'ISP+IOFFBS,IREO(',JCT-1,IRR,') = ',
     &                            ISP+IOFFBS
            end if
          end if
          if (INVERSM.eq.1.and.ISP.eq.NSTR(1)) then
            IOFFBS = NSTR(1)
          end if
        end do
        JCT = 1
      end do
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) 'Inverse reordering vector for 2-el. ints:'
        call iwrtma(INVREO,1,ITOTDIM,1,ITOTDIM)
        write(6,*) 'Reordering array for 1-el. ints (index/irrep):'
        call iwrtma(IREO,IMXDIM,NIRREP,IMXDIM,NIRREP)
      end if
*
*  ... and distribute
      call numints_red(IGLSUM,NSMBLK,NTOOBS,NIRREP,NTEST)
      call dist_2eints_dirac(INVREO,ITOTDIM2,
     &                       RKLR,KR,IOFF12,NUMINT,
     &                       INDK,INDL,TWELINT,IGLSUM,
     &                       NTOOBS,NIRREP,NTEST)
*
      ITROFF = 0
      NELMTOT = 0
      do IRR=1,NIRREP,1
        if (IRR.ge.2) then
          ITROFF = ITROFF + (NTOOBS(IRR-1)**2 + NTOOBS(IRR-1)) / 2
        end if
        NTRELM = (NTOOBS(IRR)**2 + NTOOBS(IRR)) / 2
        NELMTOT = NELMTOT + NTRELM
        if (NTEST.ge.3) then
          write(6,*) 'Calling 1ints distribution with:'
          write(6,'(A,I4)') 'Irrep:           ',IRR
          write(6,'(A,I4)') 'Triangle offset: ',ITROFF
          write(6,'(A,I4)') 'Num. of el.:     ',NTRELM
          write(6,'(A,I4)') 'Total number:    ',NELMTOT
          write(6,*) 'Partition of reordering array:'
          call iwrtma(IREO(1,IRR),NTOOBS(IRR),1,IMXDIM,NIRREP)
          if (NTEST.ge.10) then
            write(6,*) 'Full matrix of 1-electron integrals:'
            call wrtmat(ONEINTR,ITOTDIM2,ITOTDIM2,ITOTDIM2,ITOTDIM2)
          end if
        end if
        call dist_1eints_dirac2(HELPMAT,ONEINTR,
     &                         OTRINT,IREO(1,IRR),
     &                         NTRELM,ITROFF,
     &                         ITOTDIM2,ITOTDIM,
     &                         NTOOBS(IRR),NIRREP,IRR,
     &                         NTOOBS,NTEST)
      end do
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) 'Triangle list of 1-el. integrals:'
        call wrtmat(OTRINT,NELMTOT,1,NELMTOT,1)
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DMPINT(LUINT)
*
* Dump integrals in WORK(KINT1),WORK(KINT2) on file LUINT
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "glbbas.inc"
#include "wrkspc.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX
*
      CALL REWINE(LUINT,-1)
*.1 : One-electron integrals
      WRITE(LUINT,'(E22.15)')
     &     (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
      WRITE(LUINT,'(E22.15)')
     &     (WORK(KINT2-1+INT2),INT2=1,NINT2)
*.3 : Core energy
      WRITE(LUINT,'(E22.15)') ECORE_ORIG
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_CMOAO(CMO)
*
* Obtain AO-MO transformation matrix
*
* Jeppe Olsen, November 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "clunit.inc"
#include "orbinp.inc"
*. Output
      DIMENSION CMO(*)

      IF(ENVIRO(1:6).EQ.'DALTON') THEN
        CALL GET_CMOAO_DALTON(CMO,NMOS_ENV(1),NAOS_ENV(1),NSMOB)
      ELSE IF(ENVIRO(1:6).EQ.'MOLCAS') THEN
*. Read in from LUMORB file
        CALL GETMOAO_MOLCAS(CMO,LUMOIN)
      ELSE IF(ENVIRO(1:5).EQ.'LUCIA' ) THEN
*. Read in from LUCIA 1e file : unit 91
        LU91 = 91
        CALL GET_CMOAO_LUCIA(CMO,NMOS_ENV,NAOS_ENV,LU91)
      ELSE IF(ENVIRO(1:4).EQ.'NONE') THEN
        WRITE(6,*) ' GET_CMOAO, Warning : Called with ENVIRO = NONE'
        WRITE(6,*) ' No coefficients read in '
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_CMOAO_DALTON(CMO,NBAS,NMO,NSM)
*
* Obtain MO-AO expansion matrix from SIRIUS/DALTON file SIRGEOM
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER NBAS(*), NMO(*)
*. Output
      DIMENSION CMO(*)
*

      ITAP30 = 16
      OPEN(ITAP30,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND ITAP30
      CALL MOLLAB('TRCCINT ',ITAP30,6)
*. Skip record containing dimensions of orbitals
      READ(ITAP30)
*. And skip record containing eigenvalues etc
      READ(ITAP30)
C     READ (ITAP30) NSYMHF,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYMHF),
C    *              (NLAMDA(I),I=1,NSYMHF),(NORB(I),I=1,NSYMHF),
C    *              POTNUC,EMCSCF
C
C
C     READ (ITAP30) (WRK(KEIGVL+I-1),I=1,NORBT),
C    *              (IWRK(KEIGSY+I-1),I=1,NORBT)
*. And then the MO-AO expansion matrix
      NCOEF = 0
      DO ISM = 1, NSM
        NCOEF = NCOEF + NMO(ISM)*NBAS(ISM)
      END DO
      READ (ITAP30) (CMO(I),I=1,NCOEF)
      CLOSE(ITAP30,STATUS='KEEP')
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  MO - AO expansion matrix '
        WRITE(6,*) '============================='
        WRITE(6,*)
        CALL APRBLM2(CMO,NBAS,NMO,NSM,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_CMOAO_LUCIA(CMO,NMOS,NAOS,LUH)
*
* Obtain CMOAO expansion matrix from LUCIA formatted file LUH
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
*. Input
      INTEGER NMOS(*),NAOS(*)
*. Output
      DIMENSION CMO(*)
*
* Structure of file
* 1 : Number of syms
* 2 : NMO's per sym
* 3 : NAO's per SYM
* 4 : Number of elements in CMOAO
* Note : CMOAO and property integrals written in form
*     given by ONEEL_MAT_DISC
*
* Jeppe Olsen, Feb. 98
*
      WRITE(6,*)  ' GET_CMOAO_LUCIA, LUH = ', LUH
      CALL REWINE(LUH,-1)
*. skip Number of orbital symmetries
      READ(LUH,*)
*. skip Number of MO's per symmetry
      READ(LUH,*)
*. skip Number of AO's per symmetry
      READ(LUH,*)
*. skip read Length of CMO-AO expansion
      READ(LUH,*)
*. read CMO-AO expansion matrix
      CALL ONEEL_MAT_DISC(CMO,1,NSMOB,NAOS,NMOS,LUH,1)
C          ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' MO-AO transformation read in '
        CALL PRHONE(CMO,NMOS,1,NSMOB,0)
C            PRHONE(C,NFUNC,M,NSM,IPACK)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_CMOMO(CMOMO)
*
* Obtain MO-MO transformation matrix CMOMO for transforming to
* final set of orbitals
*
* Output matrix CMOMO is returned in symmetry packed form
*
*. Density matrix is assumed in place
*
* Type of final orbitals is provided by the keyword
* keywords ITRACI_CR, ITRACI_CN
*
* ITRACI_CR : COMP => Rotate all orbitals
*             REST => Rotalte only inside orbital subspaces
*
* ITRACI_CN : NATU => Transform to natural orbitals
* ITRACI_CR : CANO => Transform to canonical orbitals
*
* Jeppe Olsen, February 1998 ( from FINMO)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cgas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Output
      DIMENSION CMOMO(*)
*
      NTEST = 000
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' GET_CMOMO in action'
        WRITE(6,*) ' ===================='
        WRITE(6,*)
      END IF

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETCMO')
      CALL MEMMAN(KMAT1,NTOOB**2,'ADDL  ',2,'MAT1  ')
      CALL MEMMAN(KMAT2,NTOOB**2,'ADDL  ',2,'MAT2  ')
      CALL MEMMAN(KMAT3,NTOOB**2,'ADDL  ',2,'MAT3  ')
      CALL MEMMAN(KMAT4,2*NTOOB**2,'ADDL  ',2,'MAT4  ')
*
*. Matrix defining final orbitals
*
      IF(ITRACI_CN(1:4).EQ.'CANO' ) THEN
*. Construct FI+FA in WORK(KMAT1)
        CALL COPVEC(WORK(KINT1O),WORK(KMAT1),NINT1)
        CALL FIFAM(WORK(KMAT1))
      ELSE IF(ITRACI_CN(1:4).EQ.'NATU' ) THEN
*. Symmetry order density matrix
        CALL TYPE_TO_SYM_REO_MAT(WORK(KRHO1),WORK(KMAT2))
*. Pack to triangular form
        CALL TRIPAK_BLKM(WORK(KMAT2),WORK(KMAT1),1,NTOOBS,NSMOB)
*. multiply by minus one to get natural orbitals
*. with largest occupations first
        ONEM = -1.0D0
        LDIM = 0
        DO ISM = 1, NSMOB
          LDIM = LDIM + NTOOBS(ISM)*(NTOOBS(ISM)+1)/2
        END DO
        CALL SCALVE(WORK(KMAT1),ONEM,LDIM)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Packed density matrix ( times - 1 )'
          CALL APRBLM2(WORK(KMAT1),NACOBS,NACOBS,NSMOB,1)
        END IF
      END IF
*
* Diagonalize
*
      IF(ITRACI_CR(1:4).EQ.'REST') THEN
*. Diagonalize symmetry-type blocks
        CALL DIAG_BLKS(WORK(KMAT1),CMOMO,NGSOB,NTOOBS,MXPOBS,
     &                 NSMOB,NGAS,WORK(KMAT3),WORK(KMAT4))
*. Reorder to assure max diag dominance
        IREO = 1
        IF(IREO.NE.0) THEN
          WRITE(6,*) ' CMOMO reordered to assure max. diag. dom.'
          DO ISM = 1, NSMOB
            IF(ISM.EQ.1) THEN
              IOFF = 1
            ELSE
              IOFF = IOFF + NTOOBS(ISM-1)**2
            END IF
            L  = NTOOBS(ISM)
            CALL GET_DIAG_DOM(CMOMO(IOFF),WORK(KMAT1),L,WORK(KMAT2))
            CALL COPVEC(WORK(KMAT1),CMOMO(IOFF),L*L)
          END DO
        END IF
      ELSE IF (ITRACI_CR(1:4).EQ.'COMP') THEN
*. Diagonalize symmetry blocks
        CALL DIAG_BLKS(WORK(KMAT1),CMOMO,NACOBS,NTOOBS,MXPOBS,
     &                 NSMOB,1,WORK(KMAT3),WORK(KMAT4))
*. Reorder to assure max diag dominance
        IREO = 1
        IF(IREO.NE.0) THEN
          WRITE(6,*) ' CMOMO reordered to assure max. diag. dom.'
          DO ISM = 1, NSMOB
            IF(ISM.EQ.1) THEN
              IOFF = 1
            ELSE
              IOFF = IOFF + NTOOBS(ISM-1)**2
            END IF
            L  = NTOOBS(ISM)
            CALL GET_DIAG_DOM(CMOMO(IOFF),WORK(KMAT1),L,WORK(KMAT2))
            CALL COPVEC(WORK(KMAT1),CMOMO(IOFF),L*L)
          END DO
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Output set of MO''s '
         CALL APRBLM2(CMOMO,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETCMO')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
* Obtain property integrals with LABEL LABEL from LU91,
* LUCIA format
*
* Jeppe Olsen, Feb.98

      SUBROUTINE GET_H1AO(LABEL,H1AO,IHSM,NBAS)
*
* Obtain 1 electron integrals with label LABEL
*
* Jeppe Olsen, Feb.98
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "crun.inc"
#include "orbinp.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
*
      CHARACTER*8 LABEL
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GT_H1A')
*
      IF(ENVIRO(1:6).EQ.'DALTON') THEN
        LSCR = NTOOB**2
        CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'GTH1SC')
        CALL GET_H1AO_DALTON(LABEL,H1AO,IHSM,WORK(KLSCR),NBAS,NSMOB)
C            GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS,NSM)
      ELSE IF (ENVIRO(1:5).EQ.'LUCIA') THEN
        LU91 = 91
        CALL GET_H1AO_LUCIA(LABEL,H1AO,LU91)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GT_H1A')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS,NSM)
*
*. Obtain one-electron integrals in ao basis from dalton
*
* Label of integrals LABEL from FILE AORPROPER
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      CHARACTER*8 LABEL
      DIMENSION NBAS(*)
#include "multd2h.inc"
*. output
      DIMENSION H1AO(*)
*. Scratch
      DIMENSION SCR(*)
*
      LOGICAL FNDLAB
*
      NTEST =   02
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Fetching one-electron integrals with Label ',
     &  LABEL
        WRITE(6,*) ' IHSM NSM', IHSM,NSM
      END IF
*
*. Number of elements in : Complete lower half array
*                          Symmetry restricted complete matrix
*                          Symmetry restricted lower half matrix
*-- I am not completely sure about the input format of the integrals
      NBAST = 0
      DO ISM = 1, NSM
       NBAST = NBAST + NBAS(ISM)
      END DO
      NINT01 = NBAST*(NBAST+1)/2
C     write(6,*) ' IHSM = ', IHSM
*
      NINT10 = 0
      DO IRSM = 1, NSM
       ICSM = MULTD2H(IHSM,IRSM)
       NINT10 = NINT10 + NBAS(IRSM)*NBAS(ICSM)
      END DO
*
      NINT11 = 0
      DO IRSM = 1, NSM
       ICSM = MULTD2H(IHSM,IRSM)
       IF(IRSM.GT.ICSM) THEN
        NINT11 = NINT11 + NBAS(IRSM)*NBAS(ICSM)
       ELSE IF(IRSM.EQ.ICSM) THEN
        NINT11 = NINT11 + NBAS(IRSM)*(NBAS(IRSM)+1)/2
       END IF
      END DO
*
*. Read in integrals, assumed in complete lower half format
*
         LUPRP = 15
         OPEN (LUPRP,STATUS='OLD',FORM='UNFORMATTED',FILE='AOPROPER')
         REWIND (LUPRP)
         IF (FNDLAB(LABEL,LUPRP)) THEN
C           write(6,*) ' Label obtained'
            READ(LUPRP) (SCR(I),I=1,NINT01)
C           write(6,*) 'integrals read in'
C           call prsym(scr,NBAST)
C           CALL READT(LUPRP,NBAST*(NBAST+1)/2,WRK(KSCR2))
         ELSE
            WRITE(6,*) 'Property lable: ',LABEL ,'  not found on file'
            Call Abend2( 'Wrong input or integrals not generated' )
         ENDIF
        CLOSE(LUPRP,STATUS='KEEP')
*
C        WRITE(6,*) ' Number of symmetry apdapted integrals',NINT10
*
*. Transfer integrals to symmetry adapted form, complete form
*
         IBINT = 1
*. Loop over symmetry blocks
         DO IRSM = 1, NSM
           ICSM = MULTD2H(IHSM,IRSM)
           NR = NBAS(IRSM)
           NC = NBAS(ICSM)
*. Offsets
           IBR = 1
           DO ISM = 1, IRSM - 1
             IBR = IBR + NBAS(ISM)
           END DO
           IBC = 1
           DO ISM = 1, ICSM - 1
             IBC = IBC + NBAS(ISM)
           END DO
*. Complete block, stored in usual column wise fashion
           DO ICORB = 1, NC
             DO IRORB = 1, NR
               ICABS = IBC + ICORB -1
               IRABS = IBR + IRORB -1
               ICRMX = MAX(ICABS,IRABS)
               ICRMN = MIN(ICABS,IRABS)
               H1AO(IBINT-1 + (ICORB-1)*NR+IRORB) =
     &         SCR(ICRMX*(ICRMX-1)/2+ICRMN)
             END DO
           END DO
           IBINT = IBINT + NR*NC
         END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' One-electron integrals obtained from AOPROPER'
        CALL PRSYM(SCR,NBAST)
*
        WRITE(6,*) ' One-electron integrals in packed form'
        CALL PRHONE(H1AO,NBAS,IHSM,NSM,0)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_H1AO_LUCIA(LABEL,H1,LUH)
*
*
* Obtain property integrals with LABEL LABEL from LU91,
* LUCIA format
*
* Jeppe Olsen, Feb.98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
#include "wrkspc.inc"
C     CHARACTER*1 XYZ(3)
C     DATA XYZ/'X','Y','Z'/
      CHARACTER*8 LABEL, LABEL2, LABELX
*. Output
      DIMENSION H1(*)
*
* Structure of file
* 1 : Number of syms
* 2 : NMO's per sym
* 3 : NAO's per SYM
* 4 : Number of elements in CMOAO
* 4 : CMOAO-expansion matrix (in symmetry packed form)
* 5 : Number of property AO lists
*     Loop over number of properties
*     Label, offset and length of each proprty list
*
*     Property integrals for prop1,prop2 ...
*
* Note : CMOAO and property integrals written in form
*     given by ONEEL_MAT_DISC
*
* Jeppe Olsen, Feb. 98
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETH1A')
*
*. DIPOLE => DIPLEN
      IF(LABEL(1:6).EQ.'DIPOLE') THEN
        LABELX = 'DIPLEN  '
      ELSE
        LABELX = LABEL
      END IF
*
      CALL REWINE(LUH,-1)
*. Skip Number of orbital symmetries
      READ (LUH,*) NSMOB
*. Skip Number of MO's per symmetry
      READ (LUH,*) (NMOS_ENV(ISM),ISM=1,NSMOB)
*. Skip Number of AO's per symmetry
      READ (LUH,*) (NAOS_ENV(ISM),ISM=1,NSMOB)
*. Length of CMO-AO expansion
      READ(LUH,*) LENGTH
*. And skip
      DO IJ = 1, LENGTH
        READ(LUH,'(E22.15)')
      END DO
*. Total number of properties ( 3 for each rank1, 6 for each rank 2)
      READ(LUH,*) NPROP_COMP
      IFOUND = 0
      WRITE(6,*) ' NPROP_COMP = ', NPROP_COMP
      DO IPROP_COM = 1, NPROP_COMP
        READ(LUH,'(A,I6,I6)') LABEL2,IOFF,LENGTH
        IF(LABEL2.EQ.LABELX) THEN
          IOFFA = IOFF
          LENGTHA = LENGTH
          IFOUND = 1
        END IF
      END DO
      IF(IFOUND.EQ.0) THEN
        WRITE(6,*) ' Label not found on file 91'
        WRITE(6,'(A,A)' ) ' Label = ', LABELX
        Call Abend2( ' Label not found on file 91' )
      END IF
*. Skip to start of integrals
      WRITE(6,*) ' IOFFA, LENGTHA ', IOFFA,LENGTHA
      DO IJ = 1, IOFFA - 1
        READ(LUH,*)
      END DO
*. and read
      CALL SYM_FOR_OP(LABEL,IXYZSYM,IOPSM)
      CALL ONEEL_MAT_DISC(H1,IOPSM,NSMOB,
     &                    NAOS_ENV,NAOS_ENV,LUH,1)
C          ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Property integrals read in '
        CALL PRHONE(H1,NAOS_ENV,IOPSM,NSMOB,0)
C            PRHONE(H,NFUNC,IHSM,NSM,IPACK)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETH1A')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_ORB_DIM_ENV(ECORE_ENV)
*
* Obtain number of orbitals and basis functions from the
* programming environment.
* results stored in NAOS_ENV, NMOS_ENV
*
* Obtain environments CORE energy, ECORE_ENV
*
* Jeppe Olsen, December 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "crun.inc"
#include "lucinp.inc"
#include "parluci.h"
*
      IF(ENVIRO(1:6).EQ.'MOLCAS') THEN
C       CALL GETOBS2(ECORE_ENV,NAOS_ENV,NMOS_ENV)
        stop 'No MOLCAS ENVIRONMENT!'
      ELSE IF(ENVIRO(1:6).EQ.'DALTON' ) THEN
        CALL GETOBS_DALTON(ECORE_ENV,NAOS_ENV,NMOS_ENV)
      ELSE IF(ENVIRO(1:5).EQ.'LUCIA') THEN
*. Lucia : core energy is obtained from 2-e file
        CALL GETOBS_LUCIA(NAOS_ENV,NMOS_ENV)
      else if (ENVIRO(1:5).eq.'DIRAC') then
        if (myproc.eq.master) then
*        write(6,*)
        write(6,*) 'Integral import from DIRAC.'
        write(6,*) 'No checking of environment dims.'
        end if
      ELSE IF(ENVIRO(1:4).EQ.'NONE') THEN
*. No environment,
        WRITE(6,*) 'GET_ORB_DIM_ENV  in problems '
        WRITE(6,*) 'No ENVIRO parameter defined '
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' From GET_ORB_FROM_ENV : '
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' NAOS_ENV'
        CALL IWRTMA(NAOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' NMOS_ENV'
        CALL IWRTMA(NMOS_ENV,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' ECORE_ENV=', ECORE_ENV
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE GET_PROPINT(H,IHSM,LABEL,SCR,NMO,NBAS,NSM,ILOW)
*
*. Obtain Property integrals in MO basis for operator with
*  label LABEL.
*
* If ILOW = 1, only the elements below the diagonal are
* obtained.
*
* Jeppe Olsen, June 1997
*              September 97 : ILOW added
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NMO(*),NBAS(*)
#include "multd2h.inc"
*. Output
      DIMENSION H(*)
*. Scratch
      DIMENSION SCR(*)
*. Scratch should atleaest be of length  **
*
      NTEST = 000
*. Integrals in AO basis, neglect symmetry
      NBAST = 0
      NMOT = 0
      DO ISM = 1, NSM
        NBAST = NBAST + NBAS(ISM)
        NMOT  = NMOT  + NMO(ISM)
      END DO
C?    WRITE(6,*) ' Total number of basis functions ',NBAST
      LINTMX = NBAST*NBAST
*
      KLH1AO = 1
      KLFREE = KLH1AO + LINTMX
*
      KLC = KLFREE
      KLFREE = KLC + LINTMX
*
      KLSCR = KLFREE
*. Currently only DALTON route is working
      IDALTON = 1
      IF(IDALTON.EQ.1) THEN
C?      WRITE(6,*) ' Dalton route in action'
*. Obtain AO property integrals
C            GET_H1AO_DALTON(LABEL,H1AO,IHSM,SCR,NBAS)
C            GET_H1AO(LABEL,H1AO,IHSM,NBAS)
        CALL GET_H1AO(LABEL,SCR(KLH1AO),IHSM,NBAS)
C       CALL GET_H1AO_DALTON(LABEL,SCR(KLH1AO),IHSM,
C    &       SCR(KLSCR),NBAS,NSM)
*. Obtain MO-AO transformation matrix
        CALL GET_CMOAO(SCR(KLC))
*. Transform from AO to MO basis
C            TRAH1(NBAS,NORB,NSYM,HAO,C,HMO,IHSM,SCR)
        CALL TRAH1(NBAS,NMO,NSM,SCR(KLH1AO),SCR(KLC),H,IHSM,
     &             SCR(KLSCR))
      END IF
*
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) 'electron integrals in MO basis, full format '
        CALL PRHONE(H,NMO,IHSM,NSM,0)
      END IF
      IF(ILOW.EQ.1) THEN
*. Complete to lower half form
        IOFF_IN = 1
        IOFF_OUT = 1
        DO ISM = 1, NSM
          JSM = MULTD2H(ISM,IHSM)
          IF(ISM.EQ.JSM) THEN
*. Copy lower half
            LDIM = NMO(ISM)
            NELMNT_IN = LDIM * LDIM
            NELMNT_OUT = LDIM * (LDIM + 1)/2
            CALL COPVEC(H(IOFF_IN),SCR(KLSCR),NELMNT_IN)
            SIGN = 1.0D0
            CALL TRIPAK(SCR(KLSCR),H(IOFF_OUT),1,LDIM,LDIM,SIGN)
            IOFF_IN = IOFF_IN + NELMNT_IN
            IOFF_OUT = IOFF_OUT + NELMNT_OUT
          ELSE IF(ISM.LT.JSM) THEN
*. Just skip block in input matrix
            LIDIM = NMO(ISM)
            LJDIM = NMO(JSM)
            IOFF_IN = IOFF_IN + LIDIM*LJDIM
          ELSE IF(ISM.GT.JSM) THEN
*. Copy block to block
            LIDIM = NMO(ISM)
            LJDIM = NMO(JSM)
            NELMNT = LIDIM*LJDIM
C           CALL TRPMAT(H(IOFF_IN),LIDIM,LJDIM,H(IOFF_OUT))
            CALL COPVEC(H(IOFF_IN),H(IOFF_OUT),NELMNT)
            IOFF_IN = IOFF_IN + NELMNT
            IOFF_OUT = IOFF_OUT + NELMNT
          END IF
        END DO
      END IF
*. The one-electron integrals reside in a NMOT X NMOT matrix.
*. Zero trivial integrals
      IF(ILOW.EQ.1) THEN
        NELMNT = IOFF_OUT-1
      ELSE
        LENGTH = 0
        DO ISM = 1, NSM
          JSM = MULTD2H(ISM,IHSM)
          NELMNT = NELMNT + NMO(ISM)*NMO(JSM)
        END DO
        IFREE = NELMNT + 1
      END IF
C?    WRITE(6,*) ' GET_PROP : NELMNT= ', NELMNT
      ZERO = 0.0D0
      NZERO = NMOT*NMOT - NELMNT
      IFREE = NELMNT + 1
      CALL SETVEC(H(IFREE),ZERO,NZERO)

      IF(NTEST .GE. 50 ) THEN
        WRITE(6,*) 'electron integrals in MO basis '
        CALL PRHONE(H,NMO,IHSM,NSM,ILOW)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      Subroutine GetH0(H)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Load one electron integrals                                      *
*     File assumed opened by GETOBS                                    *
*                                                                      *
*     Calling parameters:                                              *
*     H   : core Hamiltonian matrix                                    *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)
*
#include "SysDef.inc"
*
      Parameter( LuOne = 14)
*
      COMMON/MOLOBS/
     & IOList(64),iToc(64),nBas(8),nOrb(8),nFro(8),nDel(8),
     & nSym
*
      Dimension H(*)
      Call qEnter('GetH0')
*
*----------------------------------------------------------------------*
*     Determine the number of integrals (symmetry blocked)             *
*----------------------------------------------------------------------*
      NorbTT=0
      Do iSym=1,nSym
        NorbTT=NorbTT+(Norb(iSym)*Norb(iSym)+Norb(iSym))/2
      End Do
*----------------------------------------------------------------------*
*     Load the core Hamiltonian matrix                                 *
*----------------------------------------------------------------------*
      iDisk=iToc(2)
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
      Call qExit('GetH0')
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      Subroutine GetH0S(H,NTORB)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Obtain one electron integrals                                    *
*     SIRIUS interface                                                 *
*                                                                      *
*     Calling parameters:                                              *
*     H   : core Hamiltonian matrix                                    *
*                                                                      *
*****  Author : Unknown                                *****************
*
      Implicit Real*8 (A-H,O-Z)
*
      Parameter (LUONE = 19)
*
      Dimension H(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GETH0S')
*. Square form of H0
      CALL MEMMAN(KLH0,NTORB**2,'ADDL  ',2,'KLH0  ')
*. Allocate scratch memory
      LSCR = 100000
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'H0SCR   ')
*. Get one body matrix in MO basis
      CALL INFSIR(WORK(KLSCR),WORK(KLSCR),LSCR)
C
C----------------------------------------
C     Read in the one electron integrals.
C----------------------------------------
C
      OPEN (LUONE,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *      FILE='MOONEINT')
      READ(LUONE) NCMOT,(H(I),I=1,NCMOT)
      CLOSE (LUONE,STATUS='DELETE')
*
      NTEST = 0
      IF( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' One electron integrals '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        CALL WRTMAT(H,1,NCMOT,1,NCMOT)
      END IF
C
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GETH0S')
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETH1(H,ISM,ITP,JSM,JTP)
*
* One-electron integrals over orbitals belonging to
* given OS class
*
*
* The orbital symmetries  are used to obtain the total
* symmetry of the one-electron integrals.
* It is therefore assumed that ISM, JSM represents 
*   a correct symmetry block
* of the integrals
*
* Jeppe Olsen, Version of fall 97
*              Summer of 98 : CC options added
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
*.Global pointers
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cc_exc.inc"
*.Output
      DIMENSION H(*)
*
      NI = NOBPTS(ITP,ISM)
      NJ = NOBPTS(JTP,JSM)
*
      IF(ICC_EXC.EQ.0) THEN
*
* Normal one-electron integrals
*
        IJ = 0
        DO J = 1, NJ
          DO I = 1, NI
            IJ = IJ+1
            H(IJ) = GETH1E(I,ITP,ISM,J,JTP,JSM)
          END DO
        END DO
      ELSE
*
* Single excitation coefficients dressed up as integrals
* taken from KCC
C           GET_SX_BLK(HBLK,H,IGAS,ISM,JGAS,JSM)
*. Note : WORK(KCC1) not perfect choice
       CALL GET_SX_BLK(H,WORK(KCC1),ITP,ISM,JTP,JSM)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
        CALL WRTMAT(H,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION GETH1E(IORB,ITP,ISM,JORB,JTP,JSM)
*
* One-electron integral for active
* orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
*
* The orbital symmetries are used to obtain the
* total symmetry of the operator
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "multd2h.inc"
#include "intform.inc"
*
      IJSM = MULTD2H(ISM,JSM)
      IF(IH1FORM.EQ.1) THEN
*. Normal integrals, lower triangular packed
        IF(IJSM.EQ.1) THEN
          GETH1E =
     &    GTH1ES(IREOTS,WORK(KPINT1),WORK(KINT1),IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        ELSE
          GETH1E =
     &    GTH1ES(IREOTS,WORK(KPGINT1(IJSM)),WORK(KINT1),IBSO,MXPNGAS,
     &              IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
        END IF
      ELSE
*. Integrals are in full blocked form
        GETH1E =
     &  GTH1ES(IREOTS,WORK(KPGINT1A(IJSM)),WORK(KINT1),IBSO,MXPNGAS,
     &         IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION GETH1I(IORB,JORB)
*
* Obtain one -electron integral H(IORB,JOB)
*
* Interface from EXPHAM to LUCIA
      IMPLICIT REAL*8 (A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
*
      ISM = ISMFTO(IORB)
      ITP = ITPFSO(IREOTS(IORB))
      IREL = IORB - IOBPTS(ITP,ISM) + 1
*
      JSM = ISMFTO(JORB)
      JTP = ITPFSO(IREOTS(JORB))
      JREL = JORB - IOBPTS(JTP,JSM) + 1
*
      GETH1I = GETH1E(IREL,ITP,ISM,JREL,JTP,JSM)
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' GETH1I : IORB JORB ', IORB, JORB
        WRITE(6,*) ' ISM ITP IREL ', ISM,ITP,IREL
        WRITE(6,*) ' JSM JTP JREL ', JSM,JTP,JREL
        WRITE(6,*) ' GETH1I = ', GETH1I
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETINCN2(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,I2INDX,
     &                  ICOUL)
*
* Obtain integrals
*
*     ICOUL = 0 :
*                  XINT(IK,JL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 1 :
*                  XINT(IJ,KL) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
*     ICOUL = 2 :  XINT(IL,JK) = (IJ!KL)         for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
*
* Storing for ICOUL = 1 not working if IKSM or JLSM .ne. 0
*
*
* Version for integrals stored in INTLST
*
* If type equals zero, all integrals of given type are fetched
* ( added aug8, 98)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
*. Pair of orbital indeces ( symmetry ordered ) => address in symmetry packed
*. matrix
      Dimension I2INDX(*)
*.Output
      DIMENSION XINT(*)
*. Local scratch
      DIMENSION IJARR(MXPORB)
*
      IF(ITP.GE.1) THEN
        iOrb=NOBPTS(ITP,ISM)
      ELSE
        IORB = NTOOBS(ISM)
      END IF
*
      IF(JTP.GE.1) THEN
        jOrb=NOBPTS(JTP,JSM)
      ELSE
        JORB = NTOOBS(JSM)
      END IF
*
      IF(KTP.GE.1) THEN
        kOrb=NOBPTS(KTP,KSM)
      ELSE
        KORB = NTOOBS(KSM)
      END IF
*
      IF(LTP.GE.1) THEN
        lOrb=NOBPTS(LTP,LSM)
      ELSE
        LORB = NTOOBS(LSM)
      END IF
*
*. Offsets relative to start of all orbitals, symmetry ordered
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO

*
*     Collect Coulomb terms
*
      ijblk = max(ism,jsm)*(max(ism,jsm)-1)/2 + min(ism,jsm)
      klblk = max(ksm,lsm)*(max(ksm,lsm)-1)/2 + min(ksm,lsm)
*
      IF(IJBLK.GT.KLBLK) THEN
       IJRELKL = 1
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.EQ.KLBLK) THEN
       IJRELKL = 0
       IBLOFF=IJKLOF(MAX(ISM,JSM),MIN(ISM,JSM),MAX(KSM,LSM))
      ELSE IF (IJBLK.LT.KLBLK) THEN
       IJRELKL = -1
       IBLOFF = IJKLOF(MAX(KSM,LSM),MIN(KSM,LSM),MAX(ISM,JSM))
      END IF
*
      itOrb=NTOOBS(iSm)
      jtOrb=NTOOBS(jSm)
      ktOrb=NTOOBS(kSm)
      ltOrb=NTOOBS(lSm)
*
      If(ISM.EQ.JSM) THEN
       IJPAIRS = ITORB*(ITORB+1)/2
      ELSE
       IJPAIRS = ITORB*JTORB
      END IF
*
      IF(KSM.EQ.LSM) THEN
        KLPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KLPAIRS = KTORB*LTORB
      END IF
*
      iInt=0
      Do lJeppe=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lJeppe
        Do jJeppe=jMin,jOff+jOrb-1
*
*
*. Set up array IJ*(IJ-1)/2
          IF(IJRELKL.EQ.0) THEN
            DO II = IOFF,IOFF+IORB-1
              IJ = I2INDX((JJEPPE-1)*NTOOB+II)
              IJARR(II) = IJ*(IJ-1)/2
            END DO
          END IF
*
          Do kJeppe=kOff,kOff+kOrb-1
            iMin = iOff
            kl = I2INDX(KJEPPE+(LJEPPE-1)*NTOOB)
            If(IKSM.ne.0) iMin = kJeppe
            IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &               + (KJEPPE-KOFF)*Jorb*Iorb
     &               + (JJEPPE-JOFF)*Iorb
            ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)
                IINT = (KJEPPE-KOFF)*JORB*LORB*IORB
     &               + (JJEPPE-JOFF)     *LORB*IORB
     &               + (LJEPPE-LOFF)          *IORB
            END IF
*
            IF(IJRELKL.EQ.1) THEN
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) > (ksm,lsm)
              IJKL0 = IBLOFF-1+(kl-1)*ijPairs
              IJ0 = (JJEPPE-1)*NTOOB
              Do iJeppe=iMin,iOff+iOrb-1
                  ijkl = ijkl0 + I2INDX(IJEPPE+IJ0)
                  iInt=iInt+1
                  Xint(iInt) = Intlst(ijkl)
              End Do
            END IF
*
*. block (ISM JSM !ISM JSM)
            IF(IJRELKL.EQ.0) THEN
              IJ0 = (JJEPPE-1)*NTOOB
              KLOFF = KL*(KL-1)/2
              IJKL0 = (KL-1)*IJPAIRS-KLOFF
              Do iJeppe=iMin,iOff+iOrb-1
                ij = I2INDX(IJEPPE+IJ0   )
                If ( ij.ge.kl ) Then
C                 ijkl=ij+(kl-1)*ijPairs-klOff
                  IJKL = IJKL0 + IJ
                Else
                  IJOFF = IJARR(IJEPPE)
                  ijkl=kl+(ij-1)*klPairs-ijOff
                End If
                iInt=iInt+1
                Xint(iInt) = Intlst(iblOff-1+ijkl)
              End Do
            END IF
*
*. Block (ISM JSM ! KSM LSM ) with (Ism,jsm) < (ksm,lsm)
            IF(IJRELKL.EQ.-1) THEN
              ijkl0 = IBLOFF-1+KL - KLPAIRS
              IJ0 = (JJEPPE-1)*NTOOB
              Do iJeppe=iMin,iOff+iOrb-1
                IJKL = IJKL0 + I2INDX(IJEPPE + IJ0)*KLPAIRS
                iInt=iInt+1
                Xint(iInt) = Intlst(ijkl)
              End Do
            END IF
*
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
*
      IF(ISM.EQ.LSM) THEN
       ILPAIRS = ITORB*(ITORB+1)/2
      ELSE
       ILPAIRS = ITORB*LTORB
      END IF
*
      IF(KSM.EQ.JSM) THEN
        KJPAIRS = KTORB*(KTORB+1)/2
      ELSE
        KJPAIRS = KTORB*JTORB
      END IF
*
        ilblk = max(ism,lsm)*(max(ism,lsm)-1)/2 + min(ism,lsm)
        kjblk = max(ksm,jsm)*(max(ksm,jsm)-1)/2 + min(ksm,jsm)
        IF(ILBLK.GT.KJBLK) THEN
          ILRELKJ = 1
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.EQ.KJBLK) THEN
          ILRELKJ = 0
          IBLOFF = IJKLOF(MAX(ISM,LSM),MIN(ISM,LSM),MAX(KSM,JSM))
        ELSE IF(ILBLK.LT.KJBLK) THEN
          ILRELKJ = -1
          IBLOFF = IJKLOF(MAX(KSM,JSM),MIN(KSM,JSM),MAX(ISM,LSM))
        END IF
*
        iInt=0
        Do lJeppe=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lJeppe
*
          IF(ILRELKJ.EQ.0) THEN
           DO II = IOFF,IOFF+IORB-1
             IL = I2INDX(II+(LJEPPE-1)*NTOOB)
             IJARR(II) = IL*(IL-1)/2
           END DO
          END IF
*
          Do jJeppe=jMin,jOff+jOrb-1
            Do kJeppe=kOff,kOff+kOrb-1
              KJ = I2INDX(KJEPPE+(JJEPPE-1)*NTOOB)
              KJOFF = KJ*(KJ-1)/2
              iMin = iOff
*
              IF(ICOUL.EQ.1)  THEN
*. Address before integral (1,j!k,l)
                  IINT = (LJEPPE-LOFF)*Jorb*Korb*Iorb
     &                  + (KJEPPE-KOFF)*Jorb*Iorb
     &                  + (JJEPPE-JOFF)*Iorb
              ELSE IF (ICOUL.EQ.2) THEN
*  Address before (1L,JK)
                IINT = (KJEPPE-KOFF)*JORB*LORB*IORB
     &               + (JJEPPE-JOFF)     *LORB*IORB
     &               + (LJEPPE-LOFF)          *IORB
              END IF
*
              If(IKSM.ne.0) iMin = kJeppe
*
              IF(ILRELKJ.EQ.1) THEN
                ILKJ0 = IBLOFF-1+( kj-1)*ilpairs
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE + IL0)
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.0) THEN
                IL0 = (LJEPPE-1)*NTOOB
                ILKJ0 = (kj-1)*ilPairs-kjOff
                Do iJeppe=iMin,iOff+iOrb-1
                  IL = I2INDX(IJEPPE + IL0 )
                  If ( il.ge.kj ) Then
C                     ilkj=il+(kj-1)*ilPairs-kjOff
                      ILKJ = IL + ILKJ0
                    Else
                      ILOFF = IJARR(IJEPPE)
                      ilkj=kj+(il-1)*kjPairs-ilOff
                    End If
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(iBLoff-1+ilkj)
                End Do
              END IF
*
              IF(ILRELKJ.EQ.-1) THEN
                ILKJ0 = IBLOFF-1+KJ-KJPAIRS
                IL0 = (LJEPPE-1)*NTOOB
                Do iJeppe=iMin,iOff+iOrb-1
                  ILKJ = ILKJ0 + I2INDX(IJEPPE+ IL0)*KJPAIRS
                  iInt=iInt+1
                  XInt(iInt)=XInt(iInt)-Intlst(ilkj)
                End Do
              END IF
*
            End Do
          End Do
        End Do
      End If
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE LGETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,ICOUL)

*
* Outer routine for accessing integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cc_exc.inc"
#include "crun.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*
      CALL QENTER('GETIN')
      NTEST = 00
*
      IF(NTEST.GE.5)
     &WRITE(6,*) ' GETINT : ICC_EXC and ICOUL = ', ICC_EXC, ICOUL
      IF(ICC_EXC.EQ.0) THEN
*
* =======================
* Usual/Normal  integrals
* =======================
*
*. Integrals in core in internal LUCIA format
        IF(ICOUL.NE.2) THEN
          CALL GETINCN2(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,WORK(KINT2),
     &                  WORK(KPINT2),NSMOB,WORK(KINH1),ICOUL)
        ELSE
          CALL GETINCN2(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,WORK(KINT2),
     &                  WORK(KPINT2),NSMOB,WORK(KINH1),ICOUL)
        END IF
      ELSE IF (ICC_EXC.EQ.1) THEN
*
* ============================
* Coupled Cluster coefficients
* ============================
*
        IF(ICOUL.EQ.1) THEN
          IKLJ = 0
          IJ_TRNSP = 1
        ELSE
          IKLJ = 1
          IJ_TRNSP = 0
        END IF
*. IJ_TRNSP : RSBB2BN requires blocks for e(ijkl) in the form C(ji,kl)
*. Amplitudes fetched from KCC1, KCC2 used as scratch
        CALL GET_DX_BLK(ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,WORK(KCC1+NSXE),
     &                  XINT,1,IXCHNG,IKLJ,IKSM,JLSM,WORK(KCC2),
     &                  IJ_TRNSP )
C            GET_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
C    &                  C,CBLK,IEXP,IXCHNG,IKLJ,IKSM,JLSM,SCR)
      END IF
*
      IF(NTEST.NE.0) THEN
        IF(ITP.EQ.0) THEN
          NI = NTOOBS(ISM)
        ELSE
          NI = NOBPTS(ITP,ISM)
        END IF
*
        IF(KTP.EQ.0) THEN
          NK = NTOOBS(KSM)
        ELSE
          NK = NOBPTS(KTP,KSM)
        END IF
*
        IF(IKSM.EQ.0) THEN
          NIK = NI * NK
        ELSE
          NIK = NI*(NI+1)/2
        END IF
*
        IF(JTP.EQ.0) THEN
          NJ = NTOOBS(JSM)
        ELSE
          NJ = NOBPTS(JTP,JSM)
        END IF
*
        IF(LTP.EQ.0) THEN
          NL = NTOOBS(LSM)
        ELSE
          NL = NOBPTS(LTP,LSM)
        END IF
*
        IF(JLSM.EQ.0) THEN
          NJL = NJ * NL
        ELSE
          NJL = NJ*(NJ+1)/2
        END IF
        WRITE(6,*) ' 2 electron integral block for TS blocks '
        WRITE(6,*) ' Ixchng :', IXCHNG
        WRITE(6,*) ' After GETINC '
        WRITE(6,'(1X,4(A,I2,A,I2,A))')
     &  '(',ITP,',',ISM,')','(',JTP,',',JSM,')',
     &  '(',KTP,',',KSM,')','(',LTP,',',LSM,')'
        CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
      END IF
*
      CALL QEXIT('GETIN')
C     Call Abend2( ' Jeppe forced me to stop in GETINT ' )
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETMOAO(CMOAO)
*
* Obtain MOAO matrix and save in CMOAO
*
* A sunny day in April 96
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "crun.inc"
*
      IF(INTIMP.EQ.1) THEN
*. MOLCAS environment
        WRITE(6,*) ' GETMOAO : MOLCAS environment'
        CALL GETMOAO_MOLCAS(CMOAO,LUMOIN)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETMOAO_MOLCAS(CMOAO,LU)
*
* THE MO-AO file is assumed to be a NBAS X NBAS file in LUMORB format
* as delivered by SCF or RASREAD
*
* Obtain MOAO transformation matrix from
* MOLCAS file
*
* GETOBS assumed called to define /MOLOBS/
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 TITLEMO
      COMMON/MOLOBS/
     & IOList(64),iToc(64),nBas(8),nOrb(8),nFro(8),nDel(8),
     & Nsym
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Information available in GETMOAO_MOLCAS'
        WRITE(6,*) ' ======================================='
        WRITE(6,*)
        WRITE(6,*) ' NSYM = ', NSYM
        WRITE(6,*) ' NBAS : '
        CALL IWRTMA(NBAS,1,NSYM,1,NSYM)
        WRITE(6,*) 'NORB : '
        CALL IWRTMA(NORB,1,NSYM,1,NSYM)
        WRITE(6,*) 'NFRO : '
        CALL IWRTMA(NFRO,1,NSYM,1,NSYM)
      END IF
*
      LOCC = 0
*. Full NBAS X NBAS matrix assumed, truncation only in in int transformation
      CALL RDVEC('INPROB',LU,NSYM,NBAS,NBAS,CMOAO,OCC,LOCC,TITLEMO)
      WRITE(6,*) ' Header from MOAO file (LUMOIN)'
      WRITE(6,'(80A)') TITLEMO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input MOAO transformation matrix '
        CALL APRBLM2(CMOAO,NBAS,NBAS,NSYM,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      subroutine getmolinf(WAFFCT,IPRT,NACRAS,CARD,ITEM,MXITEM,
     &                     IRETISH,IDBGRP)
*
      implicit real*8 (A-H,O-Z)
*
      character*6 WAFFCT
      character*8 FNCOM,FNIPH
      character*72 CARD
      dimension ITEM(MXITEM)
*
#include "units.inc"
*
      if (WAFFCT.eq.'HF_SCF') then
*  1) HF_SCF, input orbitals from MOLCAS-SCF calculation
*       get info from SCF
*
C        FNCOM ='COMFILE'
C        call rdcomfile(FNCOM,LUCOM,IPRT,CARD,ITEM)
C        write(6,*)
C        write(6,'(A42,A7)')
C    & ' Orbital info from file ............... ','COMFILE'
*
*  1) DHFSCF wave function, closed shell
         call info_dirac(CARD,WAFFCT,NACRAS,IDBGRP,IRETISH,IPRT)
         if (IPRT.ge.1) write(6,*)
*
      else if (WAFFCT.eq.'RASSCF') then
*  2) RASSCF, input orbitals from MOLCAS-RASSCF calculation
*       get info from JOBIPH
*
C        FNIPH = 'JOBIPH'
C        call rdiph(FNIPH,LUIPH,NACRAS,IPRT,CARD,ITEM,IRETISH)
C        write(6,*)
C        write(6,'(A42,A6)')
C    & ' Orbital info from file ............... ','JOBIPH'
*
*  2) DHFSCF wave function, open shell
         call info_dirac(CARD,WAFFCT,NACRAS,IDBGRP,IRETISH,IPRT)
         if (IPRT.ge.1) write(6,*)
      end if
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GETOBS_DALTON(ECORE_ENV,NAOS_ENV,NMOS_ENV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Scratch
      DIMENSION TITLE(24),NBAS(8), NOCC(8), NLAMDA(8), NORB(8)
*. Output
      DIMENSION NAOS_ENV(*), NMOS_ENV(*)

*
* AO info
*
C     Read information on file AONEINT from HERMIT.
      ITAP34 = 66
      OPEN (ITAP34,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
      REWIND ITAP34
      READ (ITAP34) TITLE,NST,(NBAS(I),I=1,NST),ENUC
      CLOSE(ITAP34,STATUS='KEEP')
      ECORE_ENV = ENUC
C     WRITE(6,'(//A,2(/12A6)/)')
C    *   ' Dalton   title from basis set input :',(TITLE(I),I=1,24)
*
C     WRITE(6,*) ' Number of basis functions per sym '
C     CALL IWRTMA(NBAS,NST,1,NST,1)
*
      CALL ICOPVE(NBAS,NAOS_ENV,NST)
C
C     Read information on file SIRIFC written from SIRIUS.
C
*
* MO info
*
*. By trial and error - EKD + JO, NLAMDA was identified as
*. the array holding number of MO's
*
      ITAP30 = 16
      OPEN(ITAP30,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND ITAP30
      CALL MOLLAB('TRCCINT ',ITAP30,6)
      READ (ITAP30) NSYMHF,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYMHF),
     *              (NLAMDA(I),I=1,NSYMHF),(NORB(I),I=1,NSYMHF),
     *              POTNUC,EMCSCF
      CALL ICOPVE(NLAMDA,NMOS_ENV,NST)
C?    WRITE(6,*) ' Norb as delivered from environment '
C?    CALL IWRTMA(NORB,1,8,1,8)
*
C?    WRITE(6,*) ' NOCC NLAMDA  as delivered from DALTON'
C?    CALL IWRTMA(NOCC,1,8,1,8)
C?    CALL IWRTMA(NLAMDA,1,8,1,8)
*.
C?    WRITE(6,*) ' NORBT, NCMOT = ', NORBT,NCMOT
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
       SUBROUTINE GETOBS_LUCIA(NAOS_ENV,NMOS_ENV)
*
* Obtain info on orbital dimensions from LU91 - LUCIA format
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      INTEGER NMOS_ENV(*),NAOS_ENV(*)
*
      LUH = 91
      CALL REWINE(LUH,-1)
*.
      READ(LUH,*) NSMOB
*.
      READ(LUH,*) (NMOS_ENV(ISM),ISM=1, NSMOB)
*
      READ(LUH,*) (NAOS_ENV(ISM),ISM=1, NSMOB)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION GMIJKL(IORB,JORB,KORB,LORB,INTLST,IJKLOF)
*
* Obtain integral (IORB JORB ! KORB LORB) MOLCAS version
* Integrals assumed in core
*
* Version for integrals stored in INTLST
*
      IMPLICIT REAL*8(A-H,O-Z)
*.ORBINP
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),
C    &              NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),
C    &              ITOOBS(MXPOBS)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*. Integral list
      Real * 8 Intlst(*)
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
*. The orbital list corresponds to type ordered indeces, reform to
*. symmetry ordering
*
      IABS = IREOTS(IORB)
      ISM = ISMFTO(IORB)
      IOFF = IBSO(ISM)
*
      JABS = IREOTS(JORB)
      JSM = ISMFTO(JORB)
      JOFF = IBSO(JSM)
*
      KABS = IREOTS(KORB)
      KSM = ISMFTO(KORB)
      KOFF = IBSO(KSM)
*
      LABS = IREOTS(LORB)
      LSM = ISMFTO(LORB)
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(6,*) ' GMIJKL at your service '
        WRITE(6,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(6,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(6,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(6,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
      END IF
*
      If ( jSm.gt.iSm .or. ( iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JABS - JOFF + 1
        J = IABS - IOFF + 1
      Else
        iSym=iSm
        jSym=jSm
        I = IABS - IOFF + 1
        J = JABS - JOFF + 1
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( lSm.gt.kSm  .or. ( kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LABS -LOFF + 1
        L = KABS - KOFF + 1
      Else
        kSym=kSm
        lSym=lSm
        K = KABS - KOFF + 1
        L = LABS -LOFF + 1
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(6,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      GMIJKL = Intlst(iblOff-1+ijkl)
      If( ntest .ge. 100 )
     & write(6,*) ' GMIJKL ', GMIJKL
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine gslist
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION GTIJKL(I,J,K,L)
*
* Obtain  integral (I J ! K L )
* where I,J,K and l refers to active orbitals in
* Type ordering
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*.GLobal pointers
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "crun.inc"
      IF(INTIMP .EQ. 2 ) THEN
*. LUCAS ordering
        I12S = 0
        I34S = 0
        I1234S = 1
CINA    GTIJKL = GIJKLL(IREOTS(1+NINOB),WORK(KPINT2),WORK(KLSM2),
CINA &           WORK(KINT2),
CINA &           ISMFTO,IBSO,NACOB,NSMOB,NOCOBS,I,J,K,L)
        GTIJKL = GIJKLL(IREOTS(1),WORK(KPINT2),WORK(KLSM2),
     &           WORK(KINT2),
     &           ISMFTO,IBSO,NACOB,NSMOB,NOCOBS,I,J,K,L)
       ELSE IF (INTIMP.EQ.1.OR.INTIMP.EQ.5.or.INTIMP.eq.6) THEN
*. MOLCAS OR SIRIUS IMPORT ( I hope integrals are in core !! )
          GTIJKL = GMIJKL(I,J,K,L,WORK(KINT2),WORK(KPINT2))
       END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GTJK(RJ,RK,NTOOB,SCR,IREOTS)
*
* Interface routine for obtaining Coulomb (RJ) and
* Exchange integrals (RK)
*
* Ordering of intgrals is in the internal order
      IMPLICIT REAL*8(A-H,O-Z)
*
*.CRUN
C     COMMON/CRUN/MAXIT,IRESTR,INTIMP,NP1,NP2,NQ,INCORE,MXCIV,ICISTR,
C    &            NOCSF,IDIAG
#include "mxpdim.inc"
#include "crun.inc"
*.Input
      DIMENSION IREOTS(*)
*.Output
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*.Scratch
      DIMENSION SCR(2*NTOOB ** 2)
* Commented out
C     IF(INTIMP.EQ.1) THEN
      IF(INTIMP.EQ.11) THEN
*. Interface to MOTRA integrals
        CALL GTJKM(SCR(1),SCR(1+NTOOB**2))
*. Reorder to LUCIA order
        DO 100 J = 1, NTOOB
          JJ = IREOTS(J)
          DO 50 I = 1, NTOOB
            II = IREOTS(I)
C?          WRITE(6,*) ' I II J JJ ',I,II,J, JJ
            RJ(I,J) = SCR((JJ-1)*NTOOB+II)
            RK(I,J) = SCR((JJ-1)*NTOOB+II+ NTOOB **2)
   50     CONTINUE
  100  CONTINUE
      ELSE IF(INTIMP.EQ.1.OR.INTIMP.EQ.5.or.INTIMP.eq.6) THEN
*. Interface to SIRIUS
        CALL GTJKS(RJ,RK,NTOOB)
      ELSE
*. Interface to LUCAS integrals
        CALL GTJKL(RJ,RK,NTOOB)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' RJ and RK from GTJK '
        CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
        CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GTJKL(RJ,RK,NTOOB)
*
* Obtain Coulomb  integrals (II!JJ)
*        exchange integrals (IJ!JI)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*
      DO 100 IORB = 1, NTOOB
        DO 50 JORB = 1, NTOOB
          RJ(IORB,JORB) = GTIJKL(IORB,IORB,JORB,JORB)
          RK(IORB,JORB) = GTIJKL(IORB,JORB,JORB,IORB)
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' RJ and RK from GTJK '
        CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
        CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      Subroutine GTJKM(RJ,RK)
*
*     Gather all integrals RJ(I,J) = (II!JJ)
*     Gather all integrals RK(I,J) = (IJ!IJ)
*
      IMPLICIT REAL*8(A-H,O-Z)

#include "SysDef.inc"

#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIO / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*.Output
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
*.Local
      Parameter ( lBuf    = 9600 )
      Dimension Scr(lBuf)
*
      Do iSym=1,nSmOb
        itOrb=NTOOBS(iSym)
        itOff=ITOOBS(iSym)
        iiBlk=iSym*(iSym+1)/2
        iiPairs=itOrb*(itOrb+1)/2
        Do jSym=1,iSym
          jtOrb=NTOOBS(jSym)
          jtOff=ITOOBS(jSym)
          jjBlk=jSym*(jSym+1)/2
          jjPairs=jtOrb*(jtOrb+1)/2
          ijPairs=itOrb*jtOrb
          If ( iSym.eq.jSym ) ijPairs=jtOrb+itOrb*(itOrb-1)/2
          ijBlk=jSym+iSym*(iSym-1)/2
*
*     collect all RJ(iOrb,jOrb)=(II,JJ)
*
          iRecOld=-1
          iSyBlk=jjBlk+iiBlk*(iiBlk-1)/2
          iBatch=nBatch(iSyBlk)
          iDisk=iDAdr(iBatch)
          nInts=iiPairs*jjPairs
          Do i=1,itOrb
            ii=i*(i+1)/2
*JOS
            MaxJ = jtOrb
            If(Isym.eq.Jsym) MaxJ = i
            Do j=1,MaxJ
*JOE
              jj=j*(j+1)/2
              iijj=ii+(jj-1)*iiPairs
*JOS
              If ( Isym.eq.Jsym ) Then
                jjOff=jj+(jj-1)*(jj-2)/2-1
                iijj =  iijj - jjOff
              End If
*JOE
              iRec=(iijj-1)/lTr2Rec
              If ( iRec.eq.iRecOld ) then
                iijj=iijj-iRec*lTr2Rec
              Else
                iDisk=iDAdr(iBatch)
                Do iSkip=1,iRec
                End Do
                iijj=iijj-iRec*lTr2Rec
                iRecOld=iRec
              End If
              RJ(i+itOff-1,j+jtOff-1)=Scr(iijj)
              RJ(j+jtOff-1,i+itOff-1)=Scr(iijj)
            End Do
          End Do
*
*     collect all RK(iOrb,jOrb)=(IJ,IJ)
*
          iRecOld=-1
          iSyBlk=ijBlk*(ijBlk+1)/2
          iBatch=nBatch(iSyBlk)
          iDisk=iDAdr(iBatch)
          nInts=ijPairs*(ijPairs+1)/2
          ij=0
          Do i=1,itOrb
            jMax=jtOrb
            If ( iSym.eq.jSym ) jMax=i
            Do j=1,jMax
              ij=ij+1
              ijOff=ij+(ij-1)*(ij-2)/2-1
              ijij=ij+(ij-1)*ijPairs-ijOff
              iRec=(ijij-1)/lTr2Rec
              If ( iRec.eq.iRecOld ) then
                ijij=ijij-iRec*lTr2Rec
              Else
                iDisk=iDAdr(iBatch)
                Do iSkip=1,iRec
                End Do
                ijij=ijij-iRec*lTr2Rec
                iRecOld=iRec
              End If
              RK(i+itOff-1,j+jtOff-1)=Scr(ijij)
              RK(j+jtOff-1,i+itOff-1)=Scr(ijij)
            End Do
          End Do
*
        End Do
      End Do
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Working on EXPHAM
* some known problems :
*     1 : if CSF are used diagonal is not delivered to H0mat
      SUBROUTINE GTJKS(J,K,NORB)
*
* Obtain Coulomb and Exchange integrals
* from complete integral list stored in core
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL*8           J(NORB,NORB),K(NORB,NORB)
*
      DO 200 IORB = 1, NORB
	DO 100 JORB = 1, NORB
	  J(IORB,JORB) = GTIJKL(IORB,IORB,JORB,JORB)
	  K(IORB,JORB) = GTIJKL(IORB,JORB,JORB,IORB)
  100   CONTINUE
  200 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine hello
************************************************************************
*                                                                      *
*     Print the program banner, date and time of execution             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Character*8   Fmt
      Character*120  Line,BlLine,StLine
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      Call qEnter('Hello')
      lPaper=132
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      left=(lPaper-lLine)/2
      Write(Fmt,'(A,I3.3,A)') '(',left,'X,A)'
*----------------------------------------------------------------------*
*     Print the program header                                         *
*----------------------------------------------------------------------*
      nLine=18
      Do i=1,nLine
        Line=BlLine
        If ( i.eq.1 .or. i.eq.nLine )
     &  Line=StLine
c       If ( i.eq.3 )
c    &  Line=_MOLCAS_VERSION_
        If ( i.eq.4 )
     &  Line='L U C I A'
        If ( i.eq.6 )
     &  Line='A direct CI program'
        If ( i.eq.8 )
     &  Line='author:'
        If ( i.eq.9 )
     &  Line='J. Olsen'
        If ( i.eq.10 )
     &  Line='Theoretical Chemistry, '//
     &       'University of Aarhus (Denmark)'
        If ( i.eq.12 )
     &  Line='(C) copyright, all rights reserved.'
        If ( i.eq.13 )
     &  Line='Permission is hereby granted to use but not to '
        If ( i.eq.14 )
     &  Line='reproduce or distribute any part of this program'
c       If ( i.eq.16 )
c    &  Line='module created at : '//
c    &       _BUILD_DATE_
c       Write(*,Fmt) '*'//Line//'*'
      End Do
*----------------------------------------------------------------------*
*     Print the date and time of execution                             *
*----------------------------------------------------------------------*
      Write(*,*)
      Line=BlLine
      Line(1:31)='Job run on:                 at '
      Line(13:22)=Line(41:50)
      Line(23:27)=Line(60:64)
      Line(32:39)=Line(52:59)
      Line(40:lLine)=BlLine(40:lLine)
      Write(*,Fmt) ' '//Line//' '
      Write(*,*)
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Call qExit('Hello')
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine hello_dirluc
************************************************************************
*                                                                      *
*     Print the program banner, date and time of execution             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*     Modified, Timo Fleig, Dec 2001                                   *
*                           Aug 2004                                   *
*                           Aug 2006                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
#include "clunit.inc"
      Character*8   Fmt
      Character*70  Line,BlLine,StLine
      lpaper = 72
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      Call qEnter('Hello')
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      left=(lPaper-lLine)/2
      Write(Fmt,'(A,I3.3,A)') '(',left,'X,A)'
*----------------------------------------------------------------------*
*     Print the program header                                         *
*----------------------------------------------------------------------*
      Write(LUOUT,*)
      Write(LUOUT,*)
      nLine=35
      Do i=1,nLine
        Line=BlLine
        If ( i.eq.1 .or. i.eq.2 .or. i.eq.nLine )
     &  Line=StLine
        If ( i.eq.5 )
     &  Line='D I R L U C'
        If ( i.eq.6 )
     &  Line='An interface section for LUCIA under DIRAC'
        If ( i.eq.8 )
     &  Line='author:'
        If ( i.eq.9 )
     &  Line='T. Fleig'
        If ( i.eq.10 )
     &  Line='Theoretische Chemie und Computerchemie, '
        If ( i.eq.11)
     &  Line='Heinrich-Heine-Universitaet Duesseldorf '
        If ( i.eq.12)
     &  Line='                            '
        If ( i.eq.13)
     &  Line='   Calling LUCIA version 1999'
        If ( i.eq.14)
     &  Line='     author: J. Olsen, Lund/Aarhus                '
        If ( i.eq.15)
     &  Line='                                                  '
        If ( i.eq.16)
     &  Line='     Traditional sigma vector and density modules.'
        If ( i.eq.17)
     &  Line='                                                  '
        If ( i.eq.18)
     &  Line='Citation:                                         '
        If ( i.eq.19)
     &  Line='  J. Olsen, P. Joergensen, J. Simons,             '
        If ( i.eq.20)
     &  Line='          Chem. Phys. Lett. 169 (1990) 463        '
        If ( i.eq.21)
     &  Line='  T. Fleig, L. Visscher,                          '
        If ( i.eq.22)
     &  Line='          Chem. Phys. 311 (2005) 113              '
        If ( i.eq.25)
     &  Line='Parallelization of LUCIA, Duesseldorf/Odense:     '
        If ( i.eq.26)
     &  Line='S. Knecht                                         '
        If ( i.eq.27)
     &  Line='Theoretische Chemie und Computerchemie,           '
        If ( i.eq.28)
     &  Line='Heinrich-Heine-Universitaet Duesseldorf           '
        If ( i.eq.29)
     &  Line='                                                  '
        If ( i.eq.30)
     &  Line='Citation:                                         '
        If ( i.eq.31)
     &  Line='  S. Knecht, H. J. Aa. Jensen and T. Fleig,       '
        If ( i.eq.32)
     &  Line='          J. Chem. Phys., 128 (2008) 014108       '
        If ( i.eq.33)
     &  Line='                                                  '
        If ( i.eq.34 .or. i.eq.35 )
     &  Line = StLine
        write(LUOUT,*) Line
      End Do
      Write(LUOUT,*)
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Call qExit('Hello')
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine hello_moluc
************************************************************************
*                                                                      *
*     Print the program banner, date and time of execution             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Character*8   Fmt
      Character*120  Line,BlLine,StLine
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      Call qEnter('Hello')
      lPaper=132
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      left=(lPaper-lLine)/2
      Write(Fmt,'(A,I3.3,A)') '(',left,'X,A)'
*----------------------------------------------------------------------*
*     Print the program header                                         *
*----------------------------------------------------------------------*
      nLine=18
      Do i=1,nLine
        Line=BlLine
        If ( i.eq.1 .or. i.eq.nLine )
     &  Line=StLine
c       If ( i.eq.3 )
c    &  Line=_MOLCAS_VERSION_
        If ( i.eq.5 )
     &  Line='M O L U C'
        If ( i.eq.6 )
     &  Line='An input conversion program for LUCIA under MOLCAS'
        If ( i.eq.8 )
     &  Line='author:'
        If ( i.eq.9 )
     &  Line='T. Fleig'
        If ( i.eq.10 )
     &  Line='Dept. of Theoretical Chemistry, '//
     &       'Chemical Centre, Lund (Sweden)'
        If ( i.eq.12 )
     &  Line='(C) copyright, all rights reserved.'
        If ( i.eq.13 )
     &  Line='Permission is hereby granted to use but not to '
        If ( i.eq.14 )
     &  Line='reproduce or distribute any part of this program'
c       If ( i.eq.16 )
c    &  Line='module created at : '//
c    &       _BUILD_DATE_
c       Write(*,Fmt) '*'//Line//'*'
      End Do
*----------------------------------------------------------------------*
*     Print the date and time of execution                             *
*----------------------------------------------------------------------*
      Write(*,*)
      Line=BlLine
      Line(1:31)='Job run on:                 at '
      Line(13:22)=Line(41:50)
      Line(23:27)=Line(60:64)
      Line(32:39)=Line(52:59)
      Line(40:lLine)=BlLine(40:lLine)
      Write(*,Fmt) ' '//Line//' '
      Write(*,*)
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Call qExit('Hello')
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      Function I2EAD(IORB,JORB,KORB,LORB)
*
* Find adress of integral in LUCIA order
*
      IMPLICIT REAL*8           (A-H,O-Z)
*
#include "mxpdim.inc"
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1,
C    &              KSBEVC,KSBEVL,KSBIDT,KSBCNF,KH0,KH0SCR

#include "glbbas.inc"
*
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*
      I2EAD = I2EADS(IORB,JORB,KORB,LORB,WORK(KPINT2))
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION I2EADS(IORB,JORB,KORB,LORB,IJKLOF)
*
* Obtain address of integral (IORB JORB ! KORB LORB) in MOLCAS order
* IORB JORB KORB LORB corresponds to SYMMETRY ordered indeces !!
* Integrals assumed in core
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*
      Dimension IJKLOF(NsmOB,NsmOb,NsmOB)
      Logical iSymj,kSyml,ISYMK,JSYML,ijSymkl,IKSYMJL
      Logical ijklPerm
*.
      NTEST = 000
*
      IABS = IORB
      ISM = ISMFTO(IREOST(IORB))
      IOFF = IBSO(ISM)
*
      JABS = JORB
      JSM = ISMFTO(IREOST(JORB))
      JOFF = IBSO(JSM)
*
      KABS = KORB
      KSM = ISMFTO(IREOST(KORB))
      KOFF = IBSO(KSM)
*
      LABS = LORB
      LSM = ISMFTO(IREOST(LORB))
      LOFF = IBSO(LSM)
*
      If( Ntest.ge. 100) THEN
        write(6,*) ' GMIJKL at your service '
        WRITE(6,*) ' IORB IABS ISM IOFF ',IORB,IABS,ISM,IOFF
        WRITE(6,*) ' JORB JABS JSM JOFF ',JORB,JABS,JSM,JOFF
        WRITE(6,*) ' KORB KABS KSM KOFF ',KORB,KABS,KSM,KOFF
        WRITE(6,*) ' LORB LABS LSM LOFF ',LORB,LABS,LSM,LOFF
      END IF
*
      If ( jSm.gt.iSm .or. ( iSm.eq.jSm .and. JABS.gt.IABS)) Then
        iSym=jSm
        jSym=iSm
        I = JABS - JOFF + 1
        J = IABS - IOFF + 1
      Else
        iSym=iSm
        jSym=jSm
        I = IABS - IOFF + 1
        J = JABS - JOFF + 1
      End If
      ijBlk=jSym+iSym*(iSym-1)/2
      If ( lSm.gt.kSm  .or. ( kSm.eq.lSm .and. LABS.gt.KABS)) Then
        kSym=lSm
        lSym=kSm
        K = LABS -LOFF + 1
        L = KABS - KOFF + 1
      Else
        kSym=kSm
        lSym=lSm
        K = KABS - KOFF + 1
        L = LABS -LOFF + 1
      End If
      klBlk=lSym+kSym*(kSym-1)/2
*
      ijklPerm=.false.
      If ( klBlk.gt.ijBlk ) Then
        iTemp=iSym
        iSym=kSym
        kSym=iTemp
        iTemp=jSym
        jSym=lSym
        lSym=iTemp
        iTemp=ijBlk
        ijBlk=klBlk
        klBlk=iTemp
        ijklPerm=.true.
*
        iTemp = i
        i = k
        k = itemp
        iTemp = j
        j = l
        l = iTemp
      End If
      If(Ntest .ge. 100 ) then
        write(6,*) ' i j k l ',i,j,k,l
        write(6,*) ' Isym,Jsym,Ksym,Lsym',Isym,Jsym,Ksym,Lsym
      End if
*
*  Define offset for given symmetry block
      IBLoff = IJKLof(Isym,Jsym,Ksym)
      If(ntest .ge. 100 )
     &WRITE(6,*) ' IBLoff Isym Jsym Ksym ', IBLoff,ISym,Jsym,Ksym
      iSymj=iSym.eq.jSym
      kSyml=kSym.eq.lSym
      iSymk=iSym.eq.kSym
      jSyml=jSym.eq.lSym
      ikSymjl=iSymk.and.jSyml
      ijSymkl=iSymj.and.kSyml
*
      itOrb=NTOOBS(iSym)
      jtOrb=NTOOBS(jSym)
      ktOrb=NTOOBS(kSym)
      ltOrb=NTOOBS(lSym)
C?    print *,' itOrb,jtOrb,ktOrb,ltOrb',itOrb,jtOrb,ktOrb,ltOrb
      If ( iSymj ) Then
        ijPairs=itOrb*(itOrb+1)/2
        ij=j+i*(i-1)/2
      Else
        ijPairs=itOrb*jtOrb
        ij=j + (i-1)*jtOrb
      End if
*
      IF(KSYML ) THEN
        klPairs=ktOrb*(ktOrb+1)/2
        kl=l+k*(k-1)/2
      ELSE
        klPairs=ktOrb*ltOrb
        kl=l+(k-1)*ltOrb
      End If
C?    print *,' ijPairs,klPairs',ijPairs,klPairs
*
      If ( ikSymjl ) Then
        If ( ij.gt.kl ) Then
          klOff=kl+(kl-1)*(kl-2)/2-1
          ijkl=ij+(kl-1)*ijPairs-klOff
        Else
          ijOff=ij+(ij-1)*(ij-2)/2-1
          ijkl=kl+(ij-1)*klPairs-ijOff
        End If
      Else
        ijkl=ij+(kl-1)*ijPairs
      End If
      If( ntest .ge. 100 )
     & write(6,*) ' ijkl ', ijkl
*
      I2EADS = iblOff-1+ijkl
      If( ntest .ge. 100 ) then
	write(6,*) 'i j k l ', i,j,k,l
	write(6,*) ' ibloff ijkl ',ibloff,ijkl
        write(6,*) ' I2EADS  = ', I2EADS
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      subroutine info_dirac(CARD,WAFFCT,NACRAS,IDBGRP,IRETISH,IPRT)
*
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*                                                                      *
*     Transfer spinor and space information from DIRAC                 *
*     environment to LUCITA common blocks. Decode cards.               *
*                                                                      *
*                      Timo Fleig, Jan. 2002                           *
*                                                                      *
************************************************************************
*
* MPI adaption by Stefan Knecht, Feb. 2006
* only master read in MRCONEE and distributes the information 
* to all nodes
*
      use interface_to_mpi
      Implicit Real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinf_lucita.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "parluci.h"
*
      character*6 WAFFCT
      character*8 FNIPH
      character*14 DUMCHA(MXPIRR)
      character*72 CARD,ITEM
      logical BREIT
      dimension ITEM(MXITEM),IDUMAR(MXPORB),DUMAR(MXPORB),NSPINR(2),
     &          NSTR(2),NCORE(2),NDELE(2),IDUM2(MXPIRR)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRT)
*
* Get info from DIRAC file MRCONEE
C
      IF (MYPROC.EQ.MASTER) THEN
C

      open(LUOINT,FILE='MRCONEE',FORM='UNFORMATTED')
      read(LUOINT) MORB,BREIT,ECORE_DIR,NFSYM,NZ
      read(LUOINT) IDUMMY,(DUMCHA(ID),ID=1,IDUMMY),
     &                    (IDUM2(ID),ID=1,IDUMMY),
     &                    (NSTR(I),I=1,NFSYM),
     &                    (IDUM2(ID),ID=1,NFSYM),
     &                    (NCORE(I),I=1,NFSYM),
     &                    (IDUM2(ID),ID=1,NFSYM),
     &                    (NDELE(I),I=1,NFSYM)
      read(LUOINT)
      read(LUOINT)
      read(LUOINT) (IDUMAR(ID),IDUMAR(ID),DUMAR(ID),ID=1,MORB,1),
     &             (IDUMAR(I),I=1,MORB,1),(IDUM2(ID),ID=1,NFSYM),NBREPD
      close(LUOINT)
*
      if (NTEST.ge.1)
     &  write(6,'(A,1F18.10)') '    DIRAC core energy = ',ECORE_DIR
      if (NTEST.ge.3) then
        write(6,*)
        write(6,*) 'Testing read MRCONEE'
        write(6,*) 'MORB = ',MORB
        write(6,*) 'BREIT = ',BREIT
        write(6,*) 'Number of Dirac boson symmetries: ',NBREPD
        write(6,*) 'Boson symmetries of spinors:'
        do IB=1,MORB,1
          write(6,'(A11,I5,A3,I12)') 'Spinor no. ',IB,' : ',IDUMAR(IB)
        end do
      end if
      END IF
C    /\ read in master
C
C     Distribute the information read from MRCONEE
C
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
        CALL interface_MPI_BCAST(MORB,1,MASTER,global_communicator)
        CALL interface_mpi_bcast_l0(BREIT,1,MASTER,
     &                                  global_communicator)
        CALL interface_MPI_BCAST(ECORE_DIR,1,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NFSYM,1,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NZ,1,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUMMY,1,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUM2,IDUMMY,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NSTR,NFSYM,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUM2,NFSYM,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NCORE,NFSYM,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUM2,NFSYM,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NDELE,NFSYM,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUMAR,2*MORB,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(DUMAR,MORB,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUMAR,MORB,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(IDUM2,NFSYM,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NBREPD,1,MASTER,
     &                 global_communicator)
      END IF
#endif


*
*...  space definitions
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*     Resolve number of frozen orbitals
      if (CARD(1:1).ne.'N') then
        call decode_line(CARD,72,NITEM,ITEM,MXITEM)
        do ISYM=1,NBREPD,1
          call char_to_integer_moluc(ITEM(ISYM),NFRO_INP(ISYM))
        end do
      else
        do ISYM=1,NBREPD,1
          NFRO_INP(ISYM) = 0
        end do
      end if
*
      NSYM = NBREPD
*
      call transf_dirac(IDBGRP,NSYM,NISH,NASH,NOCC,NORB,NEXT,
     &                  NFSYM,NZ,MORB,IDUMAR,NSTR,NCORE,NDELE,
     &                  NACTELD,WAFFCT,IRETISH,IPRT)
*
      NACRAS = NACTELD
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE INFSIR(IWRK,WRK,LWRK)
C
C     Written by Henrik Koch 06-09-89
*     Last modification : Jeppe Olsen, 2 sept 98.
*     NTEST flags added/modified: Timo Fleig, 16 dec 99
C
C     Purpose : Reads in information to interface to R. Harrisons
C               Full CI program and to the subsequent polarization
C               propagator calculation.
C     Argument list :
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "irat.inc"
      PARAMETER (LUINP = 17)
      PARAMETER ( ITAP30 = 16, LUPRI = 6, LUCMDS=17)
      PARAMETER ( LURSP   = 18 , LUONE = 19, LUW6 = 6)
      PARAMETER ( LBUF = 600 )
      DIMENSION BUF(LBUF), IBUF(LBUF)
      DIMENSION TITLE(24),NDEG(8),NBAS(8),NOCC(8),NEJBAS(8)
      DIMENSION IWRK(*),WRK(LWRK)
      LOGICAL LPPOP
      LOGICAL FNDLAB
      COMMON /CIPOL / NBAST, NNBAST, MORB(8), NMORBT, NORB(8),NNORB(8),
     *                NSYMHF, LBINTM, LPPOP(3,3), LSYMOP(3,3),
     *                NCMOT, NNORBX, NLAMDA(8), LUEGVC,NORBT,ISAT(128)
      common /nuclear/ enenuc
C
C
C     Read information on file AONEINT from HERMIT.
C
      ITAP34 = 66
      OPEN (ITAP34,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
      REWIND ITAP34
      READ (ITAP34) TITLE,NST,(NDEG(I),I=1,NST),ENUC
      CLOSE(ITAP34,STATUS='KEEP')
      CORE = ENUC
      enenuc = enuc
C
      NSYMHF = 0
      MXCOEF = 0
      NBFAO  = 0
      DO 100 I = 1,NST
         IF (NDEG(I) .GT. 0) THEN
            NBFAO  = NBFAO  + NDEG(I)
            NSYMHF = NSYMHF + 1
            MXCOEF = MXCOEF + NDEG(I)*NDEG(I)
         END IF
  100 CONTINUE
      WRITE(6,'(//A,2(/12A6)/)')
     *   ' Molecule title from basis set input :',(TITLE(I),I=1,24)
C?    WRITE(6,*) 'NBFAO  : ',NBFAO
C?    WRITE(6,*) 'NSYMHF : ',NSYMHF
C
C     Read information on file SIRIFC written from SIRIUS.
C
      OPEN(ITAP30,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND ITAP30
      CALL MOLLAB('TRCCINT ',ITAP30,6)
      READ (ITAP30) NSYMHF,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYMHF),
     *              (NLAMDA(I),I=1,NSYMHF),(NORB(I),I=1,NSYMHF),
     *              POTNUC,EMCSCF
      CALL ICOPVE(NLAMDA,NEJBAS,8)
C
      KEIGVL = 1
      KEIGSY = (KEIGVL + NORBT)*IRAT
      KEIGVC = (KEIGSY + NORBT - 2)/IRAT + 2
C
      READ (ITAP30) (WRK(KEIGVL+I-1),I=1,NORBT),
     *              (IWRK(KEIGSY+I-1),I=1,NORBT)
      READ (ITAP30) (WRK(KEIGVC+I-1),I=1,NCMOT)
      CLOSE(ITAP30,STATUS='KEEP')
C
      LUEGVC = 60
      OPEN (LUEGVC,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *      FILE='MOEIGVC')
      WRITE(LUEGVC) (WRK(KEIGVC+I-1),I=1,NCMOT)
      CLOSE (LUEGVC,STATUS='KEEP')
C
C     Check information from AONEINT and SIRGEOM.
C
      I_DO_CHECK = 0
      IF(I_DO_CHECK.EQ.0) THEN
        WRITE(6,*) ' Warning : No checks of consistency '
        WRITE(6,*) ' Between AOONEINT and SIRIFC '
      ELSE
      IF ((MXCOEF .NE. NCMOT) .OR. (NBFAO  .NE. NORBT) .OR.
     *    (NBAST  .NE. NBFAO) .OR. (ENUC   .NE. POTNUC)) THEN
         WRITE(LUPRI,*) 'Inconsistency error between AONEINT and'
         WRITE(LUPRI,*) 'SIRGEOM'
         WRITE(LUPRI,*) 'MXCOEF AND NCMOT',MXCOEF,NCMOT
         WRITE(LUPRI,*) 'NBFAO  AND NORBT',NBFAO,NORBT
         WRITE(LUPRI,*) 'NSYMHF AND NSYM',NSYMHF,NSYM
         WRITE(LUPRI,*) 'ENUC   AND POTNUC',ENUC,POTNUC
         Call Abend2( 'INCONSISTENCY ERROR IN LOAD' )
      ELSE
         WRITE(LUPRI,*) 'Input from AONEINT and SIRGEOM was found'
         WRITE(LUPRI,*) 'to be ok, and we thus proceed.'
      ENDIF
      END IF
C
      NTEST = 0
      IF (NTEST .GT. 10) THEN
         DO 120 I = 1,NORBT
            WRITE(LUPRI,'(/A)') 'Orbital number, symmetry and energy'
            WRITE(LUPRI,'(A/)') '-----------------------------------'
            WRITE(LUPRI,'(I3,5X,I1,5X,F16.6)')
     *           I,IWRK(KEIGSY+I-1),WRK(KEIGVL+I-1)
  120    CONTINUE
      ENDIF
C
      IF (NTEST .GT. 2) THEN
         WRITE(LUPRI,*) 'Nuclear repulsion energy : ',POTNUC
         WRITE(LUPRI,*) 'Total SCF energy         : ',EMCSCF
      ENDIF
C
      KONEEL = KEIGVC + NBAST*NBAST
      KMOONE = KONEEL + NBAST*NBAST
      KSCR1  = KMOONE + NBAST*NBAST
      KEND   = KSCR2  + 3*NBAST*NBAST
*
CTF
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) 'WARNING !!!'
      write(6,*) ' KSCR2 is undefined in infsir !! '
      write(6,*) 'KSCR2 = ',KSCR2
CTF
*
      IF ( KEND .GT. LWRK ) THEN
         Call Abend2( 'Insufficient spaces in INFSIR' )
      ENDIF
C
C     ********************************************************
C     * Read one-electron integrals and transform to MO-basis*
C     ********************************************************
C
*. Read
      OPEN (ITAP34,STATUS='OLD',FORM='UNFORMATTED',FILE='AOONEINT')
      REWIND ITAP34
      CALL MOLLAB('ONEHAMIL',ITAP34,6)
 2100 READ (ITAP34) (BUF(I),I=1,LBUF),(IBUF(I),I=1,LBUF),LENGTH
      DO 2200 I = 1,LENGTH
         WRK(KONEEL - 1 + IBUF(I)) = BUF(I)
 2200 CONTINUE
      IF (LENGTH .GE. 0) GO TO 2100
      CLOSE(ITAP34,STATUS='KEEP')
*. Transform
C     TRAN_SYM_BLOC_MAT3(AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      CALL TRAN_SYM_BLOC_MAT3(WRK(KONEEL),WRK(KEIGVC),NST,
     &     NDEG,NEJBAS,WRK(KMOONE),WRK(KSCR1),1)
C
C     WRITE 1-E MOINTS TO DISK TEMPORARILY
C
      NCOEF_MO_MO = 0
      DO ISYM = 1, NSYMHF
        NCOEF_MO_MO = NCOEF_MO_MO + NORB(ISYM)**2
      END DO
*
      OPEN (LUONE,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *      FILE='MOONEINT')
      WRITE(LUONE) NCOEF_MO_MO,(WRK(KMOONE+I-1),I=1,NCOEF_MO_MO)
      CLOSE (LUONE,STATUS='KEEP')
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE INTIM(IPRNT)
*
* Interface to external integrals
*
* If NOINT .ne. 0, only pointers are constructed
* Jeppe Olsen, Winter of 1991
*
* Version : Fall 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "lucinp.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "parluci.h"
*./CINTFO/
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*./CECORE/
      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX

*. : Pointers for symmetry blocks of integrals
*
      CALL INTPNT(WORK(KPINT1),WORK(KLSM1),
     &            WORK(KPINT2),WORK(KLSM2))
*
*. Pointer for orbital indeces for symmetry blocked matrices
C ORBINH1(IORBINH1,NTOOBS,NTOOB,NSMOB)
      CALL ORBINH1(WORK(KINH1),NTOOBS,NTOOB,NSMOB)
*
      IF(INTIMP.EQ.1.AND.NOINT.EQ.0) THEN
*
*  ==============
* . Molcas format
*  ==============
*
        WRITE(6,*) ' Integrals imported from MOLCAS files'
*.Initialize buffers, open
        CALL MKLUNDIO
*. Load one-electron integrals
        CALL GETH0(WORK(KINT1))
*. And two-electron integrals if desired
        IF(INCORE.EQ.1) THEN
          CALL INTIMM(WORK(KINT2),NSMOB)
        END IF
*
* Formatted input of symmetry non-vanishing integrals
*
      ELSE IF (INTIMP.EQ.3.AND.NOINT.EQ.0 ) THEN
        WRITE(6,*) ' Integrals imported formatted (E22.15) '
*.1 : One-electron integrals
        REWIND LUTINT
        READ(LUTINT,'(E22.15)') (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
        IF(INCORE.EQ.1.OR.EXTSPC.EQ.0) THEN
          READ(LUTINT,'(E22.15)') (WORK(KINT2-1+INT2),INT2=1,NINT2)
        END IF
*.3 : Core energy
        WRITE(6,*) ' Core energy will be read '
        READ(LUTINT,'(E22.15)') ECORE
*
*. Proceed in Normal MOLCAS way so
        INTIMP = 1
*
* ===============
* . SIRIUS format
* ===============
*
      ELSE IF (INTIMP .EQ. 5.AND.NOINT.EQ.0  ) THEN
        WRITE(6,*) ' Integrals imported from DALTON/SIRIUS files'
*. Load one-electron integrals
        CALL GETH0S(WORK(KINT1),NOCOB)
*. And two-electron integrals if desired
        IF(INCORE.EQ.1) THEN
	  LSCR = 100000
          IDUM = 2803
          CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIRFAC')
	  CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'KLSCR ')
	  CALL READMO(WORK(KLSCR),LSCR,WORK(KINT2))
          CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'SIRFAC')
        END IF
*
* ===============
* . DIRAC format
* ===============
*
      else if (INTIMP.eq.6.and.NOINT.eq.0) then
        if (myproc.eq.master) then
          write(6,*) 'Integrals imported from DIRAC files'
        else
          write(12,*) 'Integrals imported from DIRAC files'
        end if
        call rdint_dirac(KINT1,KINT2,DOUGRP,IPRNT)
      END IF
*
C     write(6,*) ' Memory check   INTIM 2 '
*. integrals (ij!kk) -(ik!kj)
C       ZIJPNT(IJPNT,IJSM,IGEJ,NOBSM,IOBSM,LOBSM,ISTOB,NTOTOB)
*... Here I am
C     CALL ZIJPNT(WORK(KPNIJ),1,1,NSMOB,IBSO,NTOOBS,IREOST,NTOOB)
C?    write(6,*) ' Home from ZIJPNT '
C    IJKK(XIJKK,IGEJ,IJPNT,NIJ,IJSM,NOBSM,LOBSM,IOBSM,NTOTOB)
C     IF(NOINT.EQ.0)
C    &CALL IJKK(WORK(KIJKK),1,WORK(KPNIJ),1,NSMOB,NTOOBS,
C    &          IBSO,IREOST,NTOOB)
C     WRITE(6,*) ' Memcheck after IJKK '
*. Change one-electron integrals to inactive fock matrix
      IF(NOINT.EQ.0) THEN
C?      WRITE(6,*) ' INTIM : IUSE_PH', IUSE_PH
        CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
        IF(IUSE_PH.EQ.1) THEN
           CALL FI(WORK(KINT1),ECORE_HEX,1)
        ELSE
           ECORE_HEX = 0.0D0
        END IF
      END IF
      ECORE_ORIG = ECORE
      ECORE      = ECORE + ECORE_HEX
      if (myproc.eq.master) then
        WRITE(6,*) 'Updated core energy ',ECORE
      end if
*
C?    WRITE(6,*) ' IDMPIN ', IDMPIN
      IF (IDMPIN.EQ.1 ) THEN
        WRITE(6,*)
     &   ' Integrals written formatted (E22.15) on unit 90'
        LU90 = 90
        REWIND LU90
*.1 : One-electron integrals
        WRITE(LU90,'(E22.15)')
     &   (WORK(KINT1O-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
        WRITE(LU90,'(E22.15)')
     &   (WORK(KINT2-1+INT2),INT2=1,NINT2)
*.3. Core energy
        WRITE(LU90,'(E22.15)')ECORE_ORIG
*.4  Rewind to empty buffer
        REWIND LU90
*.   Symmetry info etc two LU91
        LU91 = 91
        CALL DUMP_1EL_INFO(LU91)
      END IF
C?    WRITE(6,*) ' INTIM : First 10 integrals in WORK(KINT2) '
C?    CALL WRTMAT(WORK(KINT2),1,10,1,10)

C!    Call Abend2( ' Jeppe forced my to stop in INTIM ' )
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE INTIMM(XINT,MAXSYM)
*
* Import all two electron integrals from MOTRA 2e-file
*
*
      IMPLICIT REAL*8(A-H,O-Z)

#include "SysDef.inc"

*.ORBINP
#include "mxpdim.inc"
#include "orbinp.inc"
*.LUNDIO
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Common / LundIO / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*.Output
      DIMENSION XINT(*)
*. For testing
      Ntest = 0
      Ioff = 1
      Do 101 Ism = 1, Maxsym
        Do 102 Jsm = 1,Ism
          Do 103 Ksm = 1, ISm
            If(Ism .eq. Ksm ) Then
             LsmMX = Jsm
            Else
             LsmMX = Ksm
            End if
            DO 104 Lsm = 1, LsmMX
            If ( ieor(iSm-1,jSm-1).ne.ieor(kSm-1,lSm-1) )  goto 104
*
              IJbl = Ism*(Ism-1)/2 + Jsm
              KLbl = Ksm*(Ksm-1)/2 + Lsm
              IJKLbl = IJbl*(IJbl-1)/2 + KLbl
*
              nIorb = Ntoobs(Ism)
              nJorb = Ntoobs(Jsm)
              nKorb = Ntoobs(Ksm)
              nLorb = Ntoobs(Lsm)
*
              If(Ism.Eq.Jsm) Then
                nIJ = NIorb*(NIorb+1)/2
              Else
                nIJ = NIORB*NJORB
              End If
              If(Ksm.Eq.Lsm) Then
               nKL = nKorb*(nKorb+1)/2
              Else
                nKL = nKorb*nLorb
              End If
*
              If(Ism .Eq. Ksm .And. Jsm .Eq. Lsm ) Then
                nIJKL= nIJ*(nIJ+1)/2
              Else
                nIJKL= nIJ*nKL
              End if
*
              If ( Ntest .Ne. 0 ) then
                Write(6,*) ' Ism Jsm Ksm Lsm ', Ism,Jsm,Ksm,Lsm
                Write(6,*) ' Ioff, nIJKL ', Ioff,nIJKL
              End if
              Ibatch = nbatch(IJKLbl)
              iDisk=iDAdr(iBatch)
              If(Ntest.Ne.0) Write(6,*)
     &        ' IJKLbl iBatch iDisk ', IJKLbl,iBatch,iDisk
*. Loop over records
              nRec = nIJKL/lTr2Rec
              If(Nrec*lTr2Rec.Lt. nIJKL) nRec = nRec + 1
              Ioffo = Ioff
              Do 50 IRec = 1, Nrec
                If ( IRec .Ne. Nrec ) Then
                   Nintrc = lTr2Rec
                Else
                   Nintrc = nIJKL -(Nrec-1)*Ltr2Rec
                End if
                Ioff = Ioff + Nintrc
   50        Continue
*
             If( Ntest .Ne. 0 ) then
               Write(6,*) ' Integral block '
               If(Ism.Eq.Ksm .And. Jsm. Eq. Lsm ) Then
                 Call Prsym(Xint(Ioffo),nIJ)
               Else
                 Call Wrtmat(Xint(Ioffo),nIJ,nKL,nIJ,nKL)
               End if
             End if
*
  104       Continue
  103     Continue
  102   Continue
  101 Continue
*
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Create usual input file for LUCITA by writing keywords and
*  values to unit ..
*
      subroutine luciwrt(IPLOCAL,NSHPGS,MXPNGAS,MXPIRR,MXPICI,
     &                   NGAS,NSYM,NCISPC,NSEQ,TYPELUC,MAXIT,
     &                   MXSEQ,IRESDIM,LCSBLK,NACTEL,MSCOMB,SCOMB,
     &                   IRREP,NROOT,MULTIP,MS2,NELGS,MACHINE,
     &                   TITLE,MXTIT,NTIT,CALSIZ,IVFLEV,
     &                   IUSE_PH,IDEMOL,IRSTRT,IDBGRP,
     &                   IPRSTR,IPRCIX,IPRORB,IPRDIA,
     &                   IPRXT,IPRRSP,IPRDEN,IPROCC,IPRNCIV,
     &                   energy_convergence)
*
      implicit real*8 (A-H,O-Z)
*
#include "units.inc"
#include "infpar.h"
#include "parluci.h"
*
      logical MSCOMB
      character*2 TYPELUC(MXSEQ)
      character*3 CALSIZ
      character*6 MACHINE
      character*72 TITLE(MXTIT),HESTR
      CHARACTER LUCIAIN_MPIN*7, NLUCIAIN_MPI*11
      INTEGER LUCIAIN_MPIL
      dimension NSHPGS(MXPNGAS,MXPIRR),NELGS(MXPNGAS,2,MXPICI)
      real(8), intent(in) :: energy_convergence
*
*-----------------------------------------------------------------*
*  Print input for LUCITA input file           
*-----------------------------------------------------------------*
*
      if (IPLOCAL.ge.5) then
         write(6,*) '==========================================='
         write(6,*) '  Input written to file LUCIAIN     '
         write(6,*) '==========================================='
*
         write(*,'(A13)') '* &LUCIA &END'
*
         write(*,'(A5)') 'TITLE'
         NTITMN = NTIT
         do I=1,NTITMN,1
            write(*,'(1X,A70)') TITLE(I)
         end do
*
CTF      write(*,'(A7)') 'MACHINE'
CTF      write(*,'(1X,A6)') MACHINE
*
         if (CALSIZ.eq.'NOR') then
           write(*,'(A6)') 'MEGACI'
         else
           write(*,'(A6)') 'TERACI'
         end if
*
         write(*,'(A6)') 'DIRAC '
*
         write(*,'(1X,I1)') IDBGRP
*
         write(*,'(A3)') 'D2H'
*
         write(*,'(A6)') 'NIRREP'
         write(*,'(1X,I3)') NSYM
*
         write(*,'(A6)') 'IREFSM'
         write(*,'(1X,I3)') IRREP
*
         write(*,'(A5)') 'GASSH'
         write(*,'(1X,I3)') NGAS
*
         HESTR = ' '
         write (HESTR,'(A,I1,A)') '(X,',NSYM-1,'(I3,A1),I3)'
         if (NSYM.eq.1) write (HESTR,'(A)') '(1X,I3)'
         do IGAS = 1,NGAS,1
            write(*,HESTR)
     &      (NSHPGS(IGAS,ISYM),',',ISYM=1,NSYM-1),NSHPGS(IGAS,NSYM)
         end do
*
         write(*,'(A6)') 'NACTEL'
         write(*,'(1X,I3)') NACTEL
*
         write(*,'(A6)') 'GASSPC'
         write(*,'(1X,I3)') NCISPC
         do ISPC=1,NCISPC,1
            write(*,'(1X,2I3)')
     &           (NELGS(IGAS,1,ISPC),NELGS(IGAS,2,ISPC),IGAS=1,NGAS)
         end do
*
         write(*,'(A6)') 'SEQUEN'
         do ICISPC=1,NCISPC,1
            write(*,'(1X,I3)') NSEQ
            do ISEQ=1,NSEQ,1
               if (TYPELUC(ISEQ).eq.'CI') then
                  write(*,'(1X,A3,I3)') 'CI,',MAXIT
               else if (TYPELUC(ISEQ).eq.'PT') then
                  write(*,'(1X,A6,I3)') 'PERTU,',MAXIT
               else if (TYPELUC(ISEQ).eq.'VF') then
                  write(*,'(1X,A8,I3)') 'VECFREE,',IVFLEV
               end if
            end do
         end do
*
         write(*,'(A6)') 'RESTRT'
         write(*,'(1X,I3)') IRSTRT
*
         write(*,'(A5)') 'MULTS'
         write(*,'(1X,I3)') MULTIP
*
         write(*,'(A3)') 'MS2'
         write(*,'(1X,I3)') MS2
*
         if (MSCOMB) then
            write(*,'(A6)') 'MSCOMB'
            write(*,'(1X,F3.1)') SCOMB
         end if
*
         write(*,'(A5)') 'ROOTS'
         write(*,'(1X,I3)') NROOT
*
         write(*,'(A5)') 'DENSI'
         write(*,'(1X,I3)') IDEMOL
*
         write(*,'(A5)') 'NOCSF'
*
         write(*,'(A6)') 'IPRSTR'
         write(*,'(1X,I3)') IPRSTR
         write(*,'(A6)') 'IPRCIX'
         write(*,'(1X,I3)') IPRCIX
         write(*,'(A6)') 'IPRORB'
         write(*,'(1X,I3)') IPRORB
         write(*,'(A6)') 'IPRDIA'
         write(*,'(1X,I3)') IPRDIA
         write(*,'(A5)') 'IPRXT'
         write(*,'(1X,I3)') IPRXT
         write(*,'(A6)') 'IPRRSP'
         write(*,'(1X,I3)') IPRRSP
         write(*,'(A6)') 'IPRDEN'
         write(*,'(1X,I3)') IPRDEN
         write(*,'(A6)') 'IPROCC'
         write(*,'(1X,I3)') IPROCC
         write(*,'(A6)') 'PRNCIV'
         write(*,'(1X,I3)') IPRNCIV
*
         write(*,'(A6)') 'RESDIM'
         write(*,'(1X,I4)') IRESDIM
*
         write(*,'(A6)') 'LCSBLK'
         write(*,'(I20)') LCSBLK

         write(*,'(A6)') 'CICONV'
         write(*,'(1X,E22.15)') energy_convergence
        
*
         if (IUSE_PH.eq.1) then
            write(*,'(A6)') 'USE_PH'
            write(*,'(A6)') 'ADVICE'
         end if
*
         if (CALSIZ.ne.'NOR') then
            write(*,'(A6)') 'CLSSEL'
         end if
*
         write(*,'(A12)') 'END OF INPUT'
*
      end if
*
*--------------------------------------------------
* Write to file LUCIAIN
*--------------------------------------------------
      LUCIAIN_MPIN = "LUCIAIN"
#if defined (VAR_MPI)
      IF (MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NLUCIAIN_MPI,'(A7,A1,I1)') LUCIAIN_MPIN,'.',MYPROC
         LUCIAIN_MPIL=9
      ELSEIF (MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NLUCIAIN_MPI,'(A7,A1,I2)') LUCIAIN_MPIN,'.',MYPROC
         LUCIAIN_MPIL=10
      ELSEIF (MYPROC .LT. 1000) THEN  ! MPI ID has three digits
         WRITE (NLUCIAIN_MPI,'(A7,A1,I3)') LUCIAIN_MPIN,'.',MYPROC
         LUCIAIN_MPIL=11
      ELSE
         CALL QUIT("NMPROC.GT.1000! EXTEND LUCITA_NODE MODULE")
      ENDIF
#else
      NLUCIAIN_MPI=LUCIAIN_MPIN
      LUCIAIN_MPIL=7
#endif

      open(unit=LUINP,file=NLUCIAIN_MPI(1:LUCIAIN_MPIL),
     &     status='UNKNOWN',form='FORMATTED')
         Rewind(LUINP)
         write(LUINP,'(A13)') '* &LUCIA &END'
*
         write(LUINP,'(A5)') 'TITLE'
         NTITMN = NTIT
         do I=1,NTITMN,1
            write(LUINP,'(1X,A70)') TITLE(I)
         end do
*
CTF      write(LUINP,'(A7)') 'MACHINE'
CTF      write(LUINP,'(1X,A6)') MACHINE
*
         if (CALSIZ.eq.'NOR') then
           write(LUINP,'(A6)') 'MEGACI'
         else
           write(LUINP,'(A6)') 'TERACI'
         end if
*
         write(LUINP,'(A6)') 'DIRAC '
*
         write(LUINP,'(1X,I1)') IDBGRP
*
         write(LUINP,'(A3)') 'D2H'
*
         write(LUINP,'(A6)') 'NIRREP'
         write(LUINP,'(1X,I3)') NSYM
*
         write(LUINP,'(A6)') 'IREFSM'
         write(LUINP,'(1X,I3)') IRREP
*
         write(LUINP,'(A5)') 'GASSH'
         write(LUINP,'(1X,I3)') NGAS
*
         HESTR = ' '
         write (HESTR,'(A,I1,A)') '(1X,',NSYM-1,'(I3,A1),I3)'
         if (NSYM.eq.1) write (HESTR,'(A)') '(1X,I3)'
         do IGAS = 1,NGAS,1
            write(LUINP,HESTR)
     &      (NSHPGS(IGAS,ISYM),',',ISYM=1,NSYM-1),NSHPGS(IGAS,NSYM)
         end do
*
         write(LUINP,'(A6)') 'NACTEL'
         write(LUINP,'(1X,I3)') NACTEL
*
         write(LUINP,'(A6)') 'GASSPC'
         write(LUINP,'(1X,I3)') NCISPC
         do ISPC=1,NCISPC,1
            write(LUINP,'(1X,2I3)')
     &           (NELGS(IGAS,1,ISPC),NELGS(IGAS,2,ISPC),IGAS=1,NGAS)
         end do
*
         write(LUINP,'(A6)') 'SEQUEN'
         do ICISPC=1,NCISPC,1
            write(LUINP,'(1X,I3)') NSEQ
            do ISEQ=1,NSEQ,1
               if (TYPELUC(ISEQ).eq.'CI') then
                  write(LUINP,'(1X,A3,I3)') 'CI,',MAXIT
               else if (TYPELUC(ISEQ).eq.'PT') then
                  write(LUINP,'(1X,A6,I3)') 'PERTU,',MAXIT
               else if (TYPELUC(ISEQ).eq.'VF') then
                  write(LUINP,'(1X,A8,I3)') 'VECFREE,',IVFLEV
               end if
            end do
         end do
*
         if (IRSTRT.eq.1) then
           write(LUINP,'(A6)') 'RESTRT'
         end if
*
         write(LUINP,'(A5)') 'MULTS'
         write(LUINP,'(1X,I3)') MULTIP
*
         write(LUINP,'(A3)') 'MS2'
         write(LUINP,'(1X,I3)') MS2
*
         if (MSCOMB) then
            write(LUINP,'(A6)') 'MSCOMB'
            write(LUINP,'(1X,F3.1)') SCOMB
         end if
*
         write(LUINP,'(A5)') 'ROOTS'
         write(LUINP,'(1X,I3)') NROOT
*
         write(LUINP,'(A5)') 'DENSI'
         write(LUINP,'(1X,I3)') IDEMOL
*
         write(LUINP,'(A5)') 'NOCSF'
*
         write(LUINP,'(A6)') 'IPRSTR'
         write(LUINP,'(1X,I3)') IPRSTR
         write(LUINP,'(A6)') 'IPRCIX'
         write(LUINP,'(1X,I3)') IPRCIX
         write(LUINP,'(A6)') 'IPRORB'
         write(LUINP,'(1X,I3)') IPRORB
         write(LUINP,'(A6)') 'IPRDIA'
         write(LUINP,'(1X,I3)') IPRDIA
         write(LUINP,'(A5)') 'IPRXT'
         write(LUINP,'(1X,I3)') IPRXT
         write(LUINP,'(A6)') 'IPRRSP'
         write(LUINP,'(1X,I3)') IPRRSP
         write(LUINP,'(A6)') 'IPRDEN'
         write(LUINP,'(1X,I3)') IPRDEN
         write(LUINP,'(A6)') 'IPROCC'
         write(LUINP,'(1X,I3)') IPROCC
         write(LUINP,'(A6)') 'PRNCIV'
         write(LUINP,'(1X,I3)') IPRNCIV
*
         write(LUINP,'(A6)') 'RESDIM'
         write(LUINP,'(1X,I4)') IRESDIM
*
         write(LUINP,'(A6)') 'LCSBLK'
         write(LUINP,'(I20)') LCSBLK

         write(LUINP,'(A6)') 'CICONV'
         write(LUINP,'(1X,E22.15)') energy_convergence
*
         if (IUSE_PH.eq.1) then
            write(LUINP,'(A6)') 'USE_PH'
            write(LUINP,'(A6)') 'ADVICE'
         end if
*
         if (CALSIZ.ne.'NOR') then
            write(LUINP,'(A6)') 'CLSSEL'
         end if
*
         write(LUINP,'(A12)') 'END OF INPUT'
*
      close (unit=LUINP)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      Subroutine MkLundIO
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Initialize the Common /LundIO/                                   *
*                                                                      *
*     Calling parameters: none                                         *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Parameter ( mxBatch = 106  )
      Parameter ( mxSyBlk = 666  )
      Parameter ( lBlk    = 9600 )
      Parameter ( LuTwo   = 13   )
      Common / LundIO / LuTr2,lTr2Rec,iDAdr(mxBatch),nBatch(mxSyBlk)
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     First set the unit number, record length and open file           *
*     per symmetry element                                             *
*----------------------------------------------------------------------*
      LuTr2=LuTwo
      lTr2Rec=lBlk
*----------------------------------------------------------------------*
*     Load the table of disk adresses                                  *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Generate the symmetry block to batch number translation table    *
*----------------------------------------------------------------------*
      iBatch=0
      Do iSym=1,8
        Do jSym=1,iSym
          Do kSym=1,iSym
            mxlSym=kSym
            If ( kSym.eq.iSym ) mxlSym=jSym
            Do lSym=1,mxlSym
              If ( ieor(iSym-1,jSym-1).eq.ieor(kSym-1,lSym-1) ) Then
                ijPair=jSym+iSym*(iSym-1)/2
                klPair=lSym+kSym*(kSym-1)/2
                iSyBlk=klPair+ijPair*(ijPair-1)/2
                iBatch=iBatch+1
                nBatch(iSyBlk)=iBatch
              End If
            End Do
          End Do
        End Do
      End Do
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
      Return
      End
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Conversion routine from 'dummy' MOLCAS input for LUCITA to
*  true LUCITA input. This routine creates the required
*  default settings and processes the input already given by
*  the MOLCAS dummy.
*************************************************************
*
*   Modified for DIRAC environment
*      Timo Fleig, Jan. 2002
*
*************************************************************
*
      subroutine mol2luc(NTIT,NGAS,NCISPC,IPRT,CARD,CARD2,CARD3,
     &                   energy_convergence)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "luci.inc"
#include "orbinf_lucita.inc"
#include "machine.inc"
#include "irat.inc"
*
      logical MSCOMB
*
      character*2 TYPELUC(MXSEQ)
      character*72 CARD(MXPNGAS),CARD2(MXPNGAS,MXPICI),
     &             ITEM(MXITEM),ITEMX,CARD3(5)
*
      dimension NSHPGS(MXPNGAS,MXPIRR),NELGS(MXPNGAS,2,MXPICI)
      dimension IHELP(4),NINOBS(MXPIRR)
      real(8), intent(in) :: energy_convergence
*
*------------------------------------------------------------------*
*  Create default settings for LUCITA
*------------------------------------------------------------------*
*
*.................................................................
*  0) General settings
*.................................................................
*
*     if (ARCH_TYPE.eq.'IRIX64') MACHINE = '64    '
      call zirat
*
*------------------------------------------------------------------*
*  Import required information from MOLCAS or DIRAC
*------------------------------------------------------------------*
*
      call getmolinf(WAFFCT,IPRT,NACRAS,CARD3,ITEM,MXITEM,
     &               IRETISH,IDBGRP)
*
*.................................................................
*  1) Type of wave function dependencies
*.................................................................
*      none (see getmolinf above)
*
*.................................................................
*  2) Inactive space as GAS 1 (if stated in input)
*.................................................................
*
*
*  INAC
      if (INACOB.eq.1) then
        call decode_line(CARD3(2),72,NITEM,ITEM,MXITEM)
        do ISYM=1,NSYM,1
          call char_to_integer_moluc(ITEM(ISYM),NINOBS(ISYM))
        end do
        NSUM = 0
        do ISYM=1,NSYM,1
          NSUM = NSUM + NINOBS(ISYM)
        end do
      end if
*
*.................................................................
*  3) TYPE dependencies (keywords depending on type of CI chosen)
*.................................................................
*
*         ................
      if (TYPE.eq.'FCI   ') then
*         ................
*
*  NGSSH, NGAS
        if (INACOB.eq.0) then
          NGAS = 1
          do ISYM=1,NSYM,1
            NSHPGS(1,ISYM) = NORB(ISYM)
          end do
        else if (INACOB.eq.1) then
          NGAS = 2
          do ISYM=1,NSYM,1
            NSHPGS(1,ISYM) = NINOBS(ISYM)
            NSHPGS(2,ISYM) = NORB(ISYM) - NINOBS(ISYM)
          end do
        end if
*
*  NACTEL
        NTOCCH = 0
        NTOCCR = 0
        if (NACTEL.eq.-1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NTOCCH = NTOCCH + 2*NOCC(ISYM)
            end do
            NACTEL = NTOCCH
            write(6,'(A42,I3)')
     &      ' Number of active electrons ........... ',NACTEL
            write(6,*)
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NTOCCR = NTOCCR + 2*NISH(ISYM)
            end do
            NOCCIM = NTOCCR + NACRAS
            NACTEL = NOCCIM
            write(6,'(A42,I3)')
     &      ' Number of active electrons ........... ',NACTEL
            write(6,*)
          end if
        end if
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
          write(6,*) 'Number of active electrons exceeds '
          write(6,*) 'the orbital space.'
          write(6,*) "Consider Pauli's famous principle and restart!"
          Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        NCISPC = 1
        ICISPC = 1
        if (INACOB.eq.0) then
          NELGS(1,1,ICISPC) = NACTEL
          NELGS(1,2,ICISPC) = NACTEL
        else if (INACOB.eq.1) then
          NELGS(1,1,ICISPC) = 2*NSUM
          NELGS(1,2,ICISPC) = 2*NSUM
          NELGS(2,1,ICISPC) = NACTEL
          NELGS(2,2,ICISPC) = NACTEL
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Number of active electrons: ',NACTEL
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Number of active electrons: ',NACTEL
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        if (CALSIZ.eq.'NOR'.or.CALSIZ.eq.'LAR') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
        else if (CALSIZ.eq.'HUG') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
CC        TYPELUC(NSEQ) = 'VF'
CC        IVFLEV = 1               ! could be made more flexible
        end if
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  RESDIM
        if (CALSIZ.eq.'NOR') then
          IRESDIM = 1000
        else if (CALSIZ.eq.'LAR') then
          IRESDIM = 10000
        else if (CALSIZ.eq.'HUG') then
          IRESDIM = 100000
        end if
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
*              ..............
      else if (TYPE.eq.'SDCI  ') then
*              ..............
*
*  NGSSH, NGAS
        if (INACOB.eq.0) then
          if (WAFFCT.eq.'HF_SCF') then
            NGAS = 2
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NOCC(ISYM)
              NSHPGS(2,ISYM) = NEXT(ISYM)
              if (NSHPGS(2,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),2,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=3,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do
              end if
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              NGAS = 2
              do ISYM=1,NSYM,1
                NSHPGS(1,ISYM) = NASH(ISYM)
                NSHPGS(2,ISYM) = NEXT(ISYM)
                if (NSHPGS(2,ISYM).gt.MXTSOB) then
                  call gassplit(NSHPGS(1,ISYM),2,MXPNGAS,
     &                          MXTSOB,NGAS,IPRT)
                else
                  do IGAS=3,NGAS,1
                    NSHPGS(IGAS,ISYM) = 0
                  end do
                end if
              end do
            else if (IRETISH.eq.0) then
              NGAS = 3
              do ISYM=1,NSYM,1
                NSHPGS(1,ISYM) = NISH(ISYM)
                NSHPGS(2,ISYM) = NASH(ISYM)
                NSHPGS(3,ISYM) = NEXT(ISYM)
                if (NSHPGS(3,ISYM).gt.MXTSOB) then
                  call gassplit(NSHPGS(1,ISYM),3,MXPNGAS,
     &                          MXTSOB,NGAS,IPRT)
                else
                  do IGAS=4,NGAS,1
                    NSHPGS(IGAS,ISYM) = 0
                  end do
                end if
              end do
            else
              write(6,*) 'IRETISH has unallowed value.'
              write(6,*) 'IRETISH = ',IRETISH
              Call Abend2('Quitting.')
            end if
          end if
*  Inactive orbitals specified:
        else if (INACOB.eq.1) then
          if (WAFFCT.eq.'HF_SCF') then
            NGAS = 3
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NINOBS(ISYM)
              NSHPGS(2,ISYM) = NOCC(ISYM) - NINOBS(ISYM)
              if (NSHPGS(2,ISYM).lt.0) then
                write(6,*) 'Too many inactive orbitals specified'
                write(6,*) 'for symmetry ',ISYM
                write(6,*) 'Reconsider input.'
                Call Abend2( 'Quitting.' )
              end if
              NSHPGS(3,ISYM) = NEXT(ISYM)
              if (NSHPGS(3,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),3,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=4,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do
              end if
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              write(6,*) 'You have specified inactive orbitals.'
              write(6,*) 'But there are no doubly occupied'
              write(6,*) 'orbitals left!      '
              write(6,*) 'Reconsider your input or files.'
              Call Abend2('Quitting.')
            end if
            NGAS = 4
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NINOBS(ISYM)
              NSHPGS(2,ISYM) = NISH(ISYM) - NINOBS(ISYM)
              if (NSHPGS(2,ISYM).lt.0) then
                write(6,*) 'Too many inactive orbitals specified'
                write(6,*) 'for symmetry ',ISYM
                write(6,*) 'Reconsider input.'
                Call Abend2( 'Quitting.' )
              end if
              NSHPGS(3,ISYM) = NASH(ISYM)
              NSHPGS(4,ISYM) = NEXT(ISYM)
              if (NSHPGS(4,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),4,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=5,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do
              end if
            end do
          end if
        end if
*
*  NACTEL
        NTOCCH = 0
        NTOCCR = 0
        if (NACTEL.eq.-1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NTOCCH = NTOCCH + 2*NOCC(ISYM)
            end do
            NACTEL = NTOCCH
            write(6,'(A42,I3)')
     &      ' Number of active electrons ........... ',NACTEL
            write(6,*)
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NTOCCR = NTOCCR + 2*NISH(ISYM)
            end do
            NOCCIM = NTOCCR + NACRAS
            NACTEL = NOCCIM
            write(6,'(A42,I3)')
     &      ' Number of active electrons ........... ',NACTEL
            write(6,*)
          end if
        end if
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
          write(6,*) 'Number of active electrons exceeds '
          write(6,*) 'the orbital space.'
          write(6,*) "Consider Pauli's famous principle and restart!"
          Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        NCISPC = 1                  ! might well be made more flexible
        NOCCEL = 0
        if (INACOB.eq.0) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NOCC(ISYM)
            end do
            NIMN = max(0,NOCCEL-2)
            NIMX = NOCCEL
            NEMN = NOCCEL
            NEMX = NOCCEL
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NIMN
              NELGS(1,2,ICISPC) = NIMX
              do IGS=2,NGAS-1,1
                NELGS(IGS,1,ICISPC) = NEMN - 2
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              NOCCEL = NACRAS
            else
              do ISYM=1,NSYM,1
                NOCCEL = NOCCEL + 2*NISH(ISYM)
              end do
            end if
*
            NIMN = max(0,NOCCEL-2)
            NIMX = NOCCEL
*
            NAMN = max(NACTEL-2,0)
            NAMX = NACTEL
*
            NEMN = NACTEL
            NEMX = NACTEL
*
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NIMN
              NELGS(1,2,ICISPC) = NIMX
              if (IRETISH.eq.0) then
                NELGS(2,1,ICISPC) = NAMN
                NELGS(2,2,ICISPC) = NAMX
                do IGS=3,NGAS-1,1
                  NELGS(IGS,1,ICISPC) = NEMN - 2
                  NELGS(IGS,2,ICISPC) = NEMX
                end do
              end if
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMN
            end do
          end if
*  Inactive orbitals chosen:
        else if (INACOB.eq.1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NOCC(ISYM)
            end do
            NFMN = 2*NSUM
            NFMX = 2*NSUM
            NIMN = max(NFMN,NOCCEL-2)
            NIMX = NOCCEL
            NEMN = NOCCEL
            NEMX = NOCCEL
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NFMN
              NELGS(1,2,ICISPC) = NFMX
              NELGS(2,1,ICISPC) = NIMN
              NELGS(2,2,ICISPC) = NIMX
              do IGS=3,NGAS-1,1
                NELGS(IGS,1,ICISPC) = NEMN - 2
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NISH(ISYM)
            end do
*
            NFMN = 2*NSUM
            NFMX = 2*NSUM
*
            NIMN = max(NFMN,NOCCEL-2)
            NIMX = NOCCEL
*
            ININ1 = max(NACRAS-2,0)
            NAMN = ININ1 + NIMN
            ININ2 = NACRAS + 2
            NASUM = 0
            do ISYM=1,NSYM,1
              NASUM = NASUM + NASH(ISYM)
            end do
            ININ3 = min(ININ2,NASUM)
            NAMX = NOCCEL + ININ3
*
            NEMN = NACTEL
            NEMX = NACTEL
*
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NFMN
              NELGS(1,2,ICISPC) = NFMX
              NELGS(2,1,ICISPC) = NIMN
              NELGS(2,2,ICISPC) = NIMX
              NELGS(3,1,ICISPC) = NAMN
              NELGS(3,2,ICISPC) = NAMX
              do IGS=4,NGAS-1,1
                NELGS(IGS,1,ICISPC) = NEMN - 2
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          end if
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Number of active electrons: ',NACTEL
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Number of active electrons: ',NACTEL
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        if (CALSIZ.eq.'NOR'.or.CALSIZ.eq.'LAR') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
        else if (CALSIZ.eq.'HUG') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
CC        TYPELUC(NSEQ) = 'VF'
CC        IVFLEV = 1               ! could be made more flexible
        end if
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  RESDIM
        if (CALSIZ.eq.'NOR') then
          IRESDIM = 1000
        else if (CALSIZ.eq.'LAR') then
          IRESDIM = 10000
        else if (CALSIZ.eq.'HUG') then
          IRESDIM = 100000
        end if
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
*              ................
      else if (TYPE.eq.'SDTQ  ') then
*              ................
*
*  NGSSH, NGAS
        if (INACOB.eq.0) then
          if (WAFFCT.eq.'HF_SCF') then
            NGAS = 2
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NOCC(ISYM)
              NSHPGS(2,ISYM) = NEXT(ISYM)
              if (NSHPGS(2,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),2,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=3,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do 
              end if
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              NGAS = 2
              do ISYM=1,NSYM,1
                NSHPGS(1,ISYM) = NASH(ISYM)
                NSHPGS(2,ISYM) = NEXT(ISYM)
                if (NSHPGS(2,ISYM).gt.MXTSOB) then
                    call gassplit(NSHPGS(1,ISYM),2,MXPNGAS,
     &                            MXTSOB,NGAS,IPRT)
                else
                  do IGAS=3,NGAS,1
                    NSHPGS(IGAS,ISYM) = 0
                  end do
                end if
              end do
            else if (IRETISH.eq.0) then
              NGAS = 3
              do ISYM=1,NSYM,1
                NSHPGS(1,ISYM) = NISH(ISYM)
                NSHPGS(2,ISYM) = NASH(ISYM)
                NSHPGS(3,ISYM) = NEXT(ISYM)
                if (NSHPGS(3,ISYM).gt.MXTSOB) then
                  call gassplit(NSHPGS(1,ISYM),3,MXPNGAS,
     &                          MXTSOB,NGAS,IPRT)
                else
                  do IGAS=4,NGAS,1
                    NSHPGS(IGAS,ISYM) = 0
                  end do
                end if
              end do
            else
              write(6,*) 'IRETISH has unallowed value.'
              write(6,*) 'IRETISH = ',IRETISH
              Call Abend2('Quitting.')
            end if
          end if
*  Inactive orbitals specified:
        else if (INACOB.eq.1) then
          if (WAFFCT.eq.'HF_SCF') then
            NGAS = 3
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NINOBS(ISYM)
              NSHPGS(2,ISYM) = NOCC(ISYM) - NINOBS(ISYM)
              if (NSHPGS(2,ISYM).lt.0) then
                write(6,*) 'Too many inactive orbitals specified'
                write(6,*) 'for symmetry ',ISYM
                write(6,*) 'Reconsider input.'
                Call Abend2( 'Quitting.' )
              end if
              NSHPGS(3,ISYM) = NEXT(ISYM)
              if (NSHPGS(3,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),3,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=4,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do
              end if
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              write(6,*) 'You have specified inactive orbitals.'
              write(6,*) 'But there are no doubly occupied'
              write(6,*) 'orbitals left!      '
              write(6,*) 'Reconsider your input or files.'
              Call Abend2('Quitting.')
            end if
            NGAS = 4
            do ISYM=1,NSYM,1
              NSHPGS(1,ISYM) = NINOBS(ISYM)
              NSHPGS(2,ISYM) = NISH(ISYM) - NINOBS(ISYM)
              if (NSHPGS(2,ISYM).lt.0) then
                write(6,*) 'Too many inactive orbitals specified'
                write(6,*) 'for symmetry ',ISYM
                write(6,*) 'Reconsider input.'
                Call Abend2( 'Quitting.' )
              end if
              NSHPGS(3,ISYM) = NASH(ISYM)
              NSHPGS(4,ISYM) = NEXT(ISYM)
              if (NSHPGS(4,ISYM).gt.MXTSOB) then
                call gassplit(NSHPGS(1,ISYM),4,MXPNGAS,
     &                        MXTSOB,NGAS,IPRT)
              else
                do IGAS=5,NGAS,1
                  NSHPGS(IGAS,ISYM) = 0
                end do
              end if
            end do
          end if
        end if
*
*  NACTEL
        NTOCCH = 0
        NTOCCR = 0
        if (NACTEL.eq.-1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NTOCCH = NTOCCH + 2*NOCC(ISYM)
            end do
            NACTEL = NTOCCH
            write(6,'(A42,I3)')
     &      ' Number of active electrons ........... ',NACTEL
            write(6,*)
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NTOCCR = NTOCCR + 2*NISH(ISYM)
            end do
            NOCCIM = NTOCCR + NACRAS
            NACTEL = NOCCIM
            write(6,'(A42,I3)')
     &      ' Number of active electrons ........... ',NACTEL
            write(6,*)
          end if
        end if
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
          write(6,*) 'Number of active electrons exceeds '
          write(6,*) 'the orbital space.'
          write(6,*) "Consider Pauli's famous principle and restart!"
          Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        NCISPC = 1               ! might well be made more flexible
        NOCCEL = 0
        if (INACOB.eq.0) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NOCC(ISYM)
            end do
            NIMN = max(0,NOCCEL-4)
            NIMX = NOCCEL
            NEMN = NOCCEL
            NEMX = NOCCEL
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NIMN
              NELGS(1,2,ICISPC) = NIMX
              do IGS=2,NGAS-1,1
                NELGS(IGS,1,ICISPC) = max(NEMN-4,0)
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          else if (WAFFCT.eq.'RASSCF') then
            if (IRETISH.eq.-1) then
              NOCCEL = NACRAS
            else
              do ISYM=1,NSYM,1
                NOCCEL = NOCCEL + 2*NISH(ISYM)
              end do
            end if
*
            NIMN = max(0,NOCCEL-4)
            NIMX = NOCCEL
*
            NAMN = max(NACTEL-4,0)
            NAMX = NACTEL
*
            NEMN = NACTEL
            NEMX = NACTEL
*
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NIMN
              NELGS(1,2,ICISPC) = NIMX
              if (IRETISH.eq.0) then
                NELGS(2,1,ICISPC) = NAMN
                NELGS(2,2,ICISPC) = NAMX
                do IGS=3,NGAS-1,1
                  NELGS(IGS,1,ICISPC) = max(max(NEMN-4,0),NAMN+1)
                  NELGS(IGS,2,ICISPC) = NEMX
                end do
              end if
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          end if
        else if (INACOB.eq.1) then
          if (WAFFCT.eq.'HF_SCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NOCC(ISYM)
            end do
            NFMN = 2*NSUM
            NFMX = 2*NSUM
            NIMN = max(NFMN,NOCCEL-4)
            NIMX = NOCCEL
            NEMN = NOCCEL
            NEMX = NOCCEL
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NFMN
              NELGS(1,2,ICISPC) = NFMX
              NELGS(2,1,ICISPC) = NIMN
              NELGS(2,2,ICISPC) = NIMX
              do IGS=3,NGAS-1,1
                NELGS(IGS,1,ICISPC) = max(NEMN-4,0)
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          else if (WAFFCT.eq.'RASSCF') then
            do ISYM=1,NSYM,1
              NOCCEL = NOCCEL + 2*NISH(ISYM)
            end do
*
            NFMN = 2*NSUM
            NFMX = 2*NSUM
*
            NIMN = max(NFMN,NOCCEL-4)
            NIMX = NOCCEL
*
            ININ1 = max(NACRAS-4,0)
            NAMN = ININ1 + NIMN
            ININ2 = NACRAS + 4
            NASUM = 0
            do ISYM=1,NSYM,1
              NASUM = NASUM + NASH(ISYM)
            end do
            ININ3 = min(ININ2,NASUM)
            NAMX = NOCCEL + ININ3
*
            NEMN = NACTEL
            NEMX = NACTEL
*
            do ICISPC=1,NCISPC,1
              NELGS(1,1,ICISPC) = NFMN
              NELGS(1,2,ICISPC) = NFMX
              NELGS(2,1,ICISPC) = NIMN
              NELGS(2,2,ICISPC) = NIMX
              NELGS(3,1,ICISPC) = NAMN
              NELGS(3,2,ICISPC) = NAMX
              do IGS=4,NGAS-1,1
                NELGS(IGS,1,ICISPC) = max(max(NEMN-4,0),NAMN+1)
                NELGS(IGS,2,ICISPC) = NEMX
              end do
              NELGS(NGAS,1,ICISPC) = NEMN
              NELGS(NGAS,2,ICISPC) = NEMX
            end do
          end if
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Number of active electrons: ',NACTEL
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Number of active electrons: ',NACTEL
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        if (CALSIZ.eq.'NOR'.or.CALSIZ.eq.'LAR') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
        else if (CALSIZ.eq.'HUG') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
CC        TYPELUC(NSEQ) = 'VF'
CC        IVFLEV = 1               ! could be made more flexible
        end if
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  RESDIM
        if (CALSIZ.eq.'NOR') then
          IRESDIM = 1000
        else if (CALSIZ.eq.'LAR') then
          IRESDIM = 10000
        else if (CALSIZ.eq.'HUG') then
          IRESDIM = 100000
        end if
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
*              ...............
      else if (TYPE.eq.'GASCI ') then
*              ...............
*
*  NGSSH, NGAS
*   given in input for moluc; resolve
        IGSFILL = 0
        do IGAS=1,NGAS,1
          call decode_line(CARD(IGAS),72,NITEM,ITEM,MXITEM)
          ITEMX = ITEM(1)
          if (ITEMX(1:4).eq.'NONE') then
            do ISYM=1,NSYM,1
              NSHPGS(IGAS,ISYM) = 0
            end do
          else if (ITEMX(1:3).eq.'ALL'.or.
     &             ITEMX(1:4).eq.'REST') then
            if (IGSFILL.ne.0) then
              write(6,*) 'Several shell spaces defined by'
              write(6,*) 'ALL or REST.'
              write(6,*) 'This is at least confusing.'
              Call Abend2( 'Quitting.' )
            end if
            IGSFILL = IGAS
          else
            if (NITEM.ne.NSYM) then
              write(6,*) 'Erroneous input to GASS.'
              write(6,'(A72)') CARD(IGAS)
              write(6,*) 'Specify either:  NONE'
              write(6,*) '                  ALL'
              write(6,*) '                 REST'
              write(6,*) 'or NSYM integers.'
              write(6,*)
              write(6,*) 'NSYM has the value ',NSYM
              Call Abend2( 'Quitting.' )
            end if
            do ISYM=1,NSYM,1
              call char_to_integer_moluc(ITEM(ISYM),NSHPGS(IGAS,ISYM))
            end do
          end if
          if (IGAS.eq.NGAS.and.ITEMX.eq.'REST') then
            do ISYM=1,NSYM,1
              NSHELLS = 0
              do IGS=1,NGAS-1,1
                NSHELLS = NSHELLS + NSHPGS(IGS,ISYM)
              end do
              NSHPGS(NGAS,ISYM) = NORB(ISYM) - NSHELLS
            end do
          end if
        end do
*
*  NACTEL
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
          write(6,*) 'Number of active electrons exceeds '
          write(6,*) 'the orbital space.'
          write(6,*) "Consider Pauli's famous principle and restart!"
          Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        do ICISPC=1,NCISPC,1
          do IGAS=1,NGAS,1
            read(CARD2(IGAS,ICISPC),*) (NELGS(IGAS,I,ICISPC),I=1,2,1)
          end do
        end do
        if (NACTEL.ne.NELGS(NGAS,2,NCISPC)) then
          write(6,*) 'I am sure that fascinating discussions about'
          write(6,*) 'the energy of such a wave function exist,'
          write(6,*) 'but I am just a dumb program, so I will stop.'
          write(6,*)
          write(6,*) 'Number of active electrons does not match'
          write(6,*) 'total number of electrons in GAS spaces.'
          Call Abend2( 'Quitting.' )
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        if (CALSIZ.eq.'NOR'.or.CALSIZ.eq.'LAR') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
        else if (CALSIZ.eq.'HUG') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
CC        TYPELUC(NSEQ) = 'VF'
CC        IVFLEV = 1               ! could be made more flexible
        end if
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  RESDIM
        if (CALSIZ.eq.'NOR') then
          IRESDIM = 1000
        else if (CALSIZ.eq.'LAR') then
          IRESDIM = 10000
        else if (CALSIZ.eq.'HUG') then
          IRESDIM = 100000
        end if
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
*              ................
      else if (TYPE.eq.'RASCI ') then
*              ................
*
*  NGSSH, NGAS
*   given in input for moluc; resolve
        do ISYM=1,NSYM,1
          NSHPGS(1,ISYM) = NINOBS(ISYM)
        end do
*
        call decode_line(CARD3(3),72,NITEM,ITEM,MXITEM)
        do ISYM=1,NSYM,1
          call char_to_integer_moluc(ITEM(ISYM),NSHPGS(2,ISYM))
        end do
*
        call decode_line(CARD3(4),72,NITEM,ITEM,MXITEM)
        do ISYM=1,NSYM,1
          call char_to_integer_moluc(ITEM(ISYM),NSHPGS(3,ISYM))
        end do
*
        call upcase(CARD3(5))
        if (CARD3(5).eq.'REST') then
          NGAS = 4
          IORBCT = 0
          do ISYM=1,NSYM,1
            NSHPGS(4,ISYM) = NORB(ISYM) - NINOBS(ISYM) -
     &                       NSHPGS(2,ISYM) - NSHPGS(3,ISYM)
            if (NSHPGS(4,ISYM).gt.MXTSOB) then
              call gassplit(NSHPGS(1,ISYM),4,MXPNGAS,
     &                      MXTSOB,NGAS,IPRT)
            end if
          end do
        else
          NGAS = 5
          call decode_line(CARD3(5),72,NITEM,ITEM,MXITEM)
          do ISYM=1,NSYM,1
            call char_to_integer_moluc(ITEM(ISYM),NSHPGS(4,ISYM))
            NSHPGS(5,ISYM) = NORB(ISYM) - NINOBS(ISYM) -
     &                       NSHPGS(2,ISYM) - NSHPGS(3,ISYM)
     &                       - NSHPGS(4,ISYM)
            if (NSHPGS(5,ISYM).gt.MXTSOB) then
              call gassplit(NSHPGS(1,ISYM),5,MXPNGAS,
     &                      MXTSOB,NGAS,IPRT)
            end if
          end do
        end if
*
*  NACTEL
        NTOOB = 0
        do ISYM=1,NSYM,1
          NTOOB = NTOOB + NORB(ISYM)
        end do
        MXTNEL = 2*NTOOB
        if (NACTEL.gt.MXTNEL) then
          write(6,*) 'Number of active electrons exceeds '
          write(6,*) 'the orbital space.'
          write(6,*) "Consider Pauli's famous principle and restart!"
          Call Abend2( 'Quitting.' )
        end if
*
*  GASSPC, NCISPC
        ICISPC = 1
        NCISPC = 1
        MXELR1 = 0
        MXELR2 = 0
        do ISYM=1,NSYM,1
          MXELR1 = MXELR1 + 2*NSHPGS(2,ISYM)
          MXELR2 = MXELR2 + 2*NSHPGS(3,ISYM)
        end do
*
        NFMN = 2*NSUM
        NFMX = 2*NSUM
        NIMX = MXELR1 + NFMX
        NIMN = NIMX - MXHOL1
        NEMN = NACTEL + NFMX
        NEMX = NACTEL + NFMX
        NAMX = NEMX
        NAMN = NAMX - MXELR3
*
        do ICISPC=1,NCISPC,1
          NELGS(1,1,ICISPC) = NFMN
          NELGS(1,2,ICISPC) = NFMX
          NELGS(2,1,ICISPC) = NIMN
          NELGS(2,2,ICISPC) = NIMX
          NELGS(3,1,ICISPC) = NAMN
          NELGS(3,2,ICISPC) = NAMX
          NELGS(4,1,ICISPC) = NEMN
          NELGS(4,2,ICISPC) = NEMX
          if (NGAS.eq.5) then
            do IGAS = 5,NGAS,1
              NELGS(IGAS,1,ICISPC) = NACTEL + NFMX
              NELGS(IGAS,2,ICISPC) = NACTEL + NFMX
            end do
          end if
        end do
        NACTEL = NACTEL+NFMX
        if (NACTEL.ne.NELGS(NGAS,2,NCISPC)) then
          write(6,*) 'I am sure that fascinating discussions about'
          write(6,*) 'the energy of such a wave function exist,'
          write(6,*) 'but I am just a dumb program, so I will stop.'
          write(6,*)
          write(6,*) 'Number of active electrons does not match'
          write(6,*) 'total number of electrons in GAS spaces.'
          Call Abend2( 'Quitting.' )
        end if
*
*  Omit first GAS if no INACTIVE orbitals given
        INUMBER = 0
        do ISYM=1,NSYM,1
          INUMBER = max(INUMBER,NSHPGS(1,ISYM))
        end do
        if (INUMBER.eq.0) then
          NGAS = NGAS - 1
          do IGAS=1,NGAS,1
            do ISYM=1,NSYM,1
              NSHPGS(IGAS,ISYM) = NSHPGS(IGAS+1,ISYM)
              do ICISPC=1,NCISPC,1
                NELGS(IGAS,1,ICISPC) = NELGS(IGAS+1,1,ICISPC)
                NELGS(IGAS,2,ICISPC) = NELGS(IGAS+1,2,ICISPC)
              end do
            end do
          end do
        end if
*
*  MULTS (just a check for consistence)
        call evenodd(IEONAC,NACTEL)
        call evenodd(IEOMUL,MULTIP)
        if (IEONAC.eq.2.and.IEOMUL.eq.2) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
           write(6,*) 'Illegal spin multiplicity given.'
           write(6,*) 'Read a book about fermions.'
           Call Abend2( 'quitting' )
        else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
           if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
           if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
              write(6,*) 'Illegal spin multiplicity given.'
              write(6,*) 'Compare with number of active'
              write(6,*) 'electrons.'
              Call Abend2( 'quitting' )
           end if
        end if
*
*  MS2
        MS2 = MULTIP - 1           ! MS(MAX) = S set by default
*
*  MSCOMB      Spin combinations (sign for)
        MSCOMB = .false.
        SCOMB = 0.
        if (MULTIP.ge.2.and.MS2.eq.0) then
          SCOMB = 1.0
          MSCOMB = .true.
        end if
*
*  SEQUEN
*  We will support CI, PERTUrbation theory, and VECFREE calculations
        NSEQ = 1
        TYPELUC(NSEQ) = 'CI'
        if (CALSIZ.eq.'NOR'.or.CALSIZ.eq.'LAR') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
        else if (CALSIZ.eq.'HUG') then
          TYPELUC(NSEQ) = 'CI'
          MAXIT = IMAXCIITER
CC        TYPELUC(NSEQ) = 'VF'
CC        IVFLEV = 1               ! could be made more flexible
        end if
*
*  DENSIty matrices
*   IDEMOL will be written directly to output
*
*  RESDIM
        if (CALSIZ.eq.'NOR') then
          IRESDIM = 1000
        else if (CALSIZ.eq.'LAR') then
          IRESDIM = 10000
        else if (CALSIZ.eq.'HUG') then
          IRESDIM = 100000
        end if
*
*  LCSBLK
        LCSBLK = 100000
*
*  USE_PH
        IUSE_PH = 1
*
*
      else
        write(6,*) 'Unknown TYPE of CI statement.'
        write(6,*) 'You have chosen  ',TYPE
        write(6,*) 'I will refrain from doing this!'
        Call Abend2( 'Quitting.' )
      end if
*
      
!     disable particle-hole symmetry for the time being - needs further attention, something goes wrong for the As atom
!     in the 4S ground state
!     with the following GAS setup (13 determinants, easy to debug --> compare with DIRRCI) - stefan Nov 2012
!     .NACTEL                                                                                             
!      15                                                                                                 
!     .GASSHE                                                                                             
!      2                                                                                                  
!      2, 0, 0, 1, 0, 1, 1, 0                                                                              
!      1, 1, 1, 0, 1, 0, 0, 0                                                                              
!     .GASSPC                                                                                             
!      1                                                                                                  
!       9 10                                                                                             
!      15 15 
!     .SYMMETRY
!      8
!     .MULTIP
!      4
!
!     FIXME: enable p-h symmetry when bug has been found

      IUSE_PH = 0

*.................................................................
*  4) Number of roots
*.................................................................
*      on common block; passed to luciwrt directly
*
*.................................................................
*  5) State symmetry
*.................................................................
*      on common block; passed to luciwrt directly
*
*.................................................................
*  6) Print command
*.................................................................
*      resolve global print flags to LUCIA system
*
*  Possible print flags
*   IPRSTR (done)
*   IPRCIX (done)
*   IPRORB (done)
*   IPRDIA (done, but a little incomplete)
*   IPRXT  (does not show up anywhere)
*   IPRRSP (for response jobs; omitted)
*   IPRDEN (done)
*   IPROCC (done)
*   IPRNCIV (done)
*
*
      IPRDEN = 0
      if (PRILUC.eq.'NON') then
         IPRSTR = 0
         IPRCIX = 0
         IPRORB = 0
         IPRDIA = 0
         IPRXT  = 0
         IPRRSP = 0
         if (IDEMOL.ge.1) IPRDEN = 1
         IPROCC = 0
         IPRNCIV = 0
      else if (PRILUC.eq.'LOW') then
         IPRSTR = 0
         IPRCIX = 0
         IPRORB = 0
         IPRDIA = 0
         IPRXT  = 0
         IPRRSP = 0
         if (IDEMOL.ge.1) IPRDEN = 1
         IPROCC = 0
         IPRNCIV = 1
      else if (PRILUC.eq.'MED') then
         IPRSTR = 10
         IPRCIX = 3
         IPRORB = 1
         IPRDIA = 2
         IPRXT  = 0
         IPRRSP = 0
         if (IDEMOL.ge.1) IPRDEN = 5
         IPROCC = 0
         IPRNCIV = 1
      else if (PRILUC.eq.'HIG') then
         IPRSTR = 20
         IPRCIX = 5
         IPRORB = 2
         IPRDIA = 10
         IPRXT  = 0
         IPRRSP = 0
         if (IDEMOL.ge.1) IPRDEN = 5
         IPROCC = 0
         IPRNCIV = 2
      else if (PRILUC.eq.'VHI') then
         IPRSTR = 500
         IPRCIX = 100
         IPRORB = 10
         IPRDIA = 200
         IPRXT  = 0
         IPRRSP = 0
         if (IDEMOL.ge.1) IPRDEN = 500
         IPROCC = 00
         IPRNCIV = 5
      end if
*
      if (IDEMOL.ge.1.and.IPRDEN.lt.1) IPRDEN = 1
*
*--------------------------------------------------------------*
*  Create input file for LUCIA
*--------------------------------------------------------------*
      call luciwrt(IPRT,NSHPGS,MXPNGAS,MXPIRR,MXPICI,
     &             NGAS,NSYM,NCISPC,NSEQ,TYPELUC,MAXIT,
     &             MXSEQ,IRESDIM,LCSBLK,NACTEL,MSCOMB,SCOMB,
     &             IRREP,NROOT,MULTIP,MS2,NELGS,MACHINE,
     &             TITLE,MXTIT,NTIT,CALSIZ,IVFLEV,
     &             IUSE_PH,IDEMOL,IRSTRT,IDBGRP,
     &             IPRSTR,IPRCIX,IPRORB,IPRDIA,IPRXT,IPRRSP,
     &             IPRDEN,IPROCC,IPRNCIV,energy_convergence)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Determine max. number of 2-el. integrals of the Dirac-Coulomb
* operator depending on symmetry information (no index permutation symmetry!)
*
      subroutine nsymint_dc(NCT,NTOOBS,NIRREP,NORB,IPRINT)
*
*-----------------------------------------
*  Timo Fleig, May 2001
*-----------------------------------------
*
      implicit real*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "multd2h.inc"
*
      dimension NTOOBS(NIRREP),ISYMOB(MXPORB)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      ICT = 0
      do ISM=1,NIRREP,1
        do IOB=ICT+1,ICT+NTOOBS(ISM),1
          ISYMOB(IOB) = ISM
        end do
        ICT = ICT + NTOOBS(ISM)
      end do
      if (NTEST.ge.3) then
        write(6,*) 'ISYMOB array:'
        call iwrtma(ISYMOB,NORB,1,NORB,1)
      end if
*
      NCT = 0
      do JORB=1,NORB,1
        do IORB=JORB,NORB,1
          ISYM = ISYMOB(IORB)
          JSYM = ISYMOB(JORB)
          KLSYM = MULTD2H(ISYM,JSYM)
          do KSYM=1,NIRREP,1
            LSYM = MULTD2H(KLSYM,KSYM)
            if (NTEST.ge.3) then
              write(6,'(A,4I3)')
     &              'Index symmetries : ',ISYM,JSYM,KSYM,LSYM
            end if
            do KORB=1,NTOOBS(KSYM),1
              do LORB=1,NTOOBS(LSYM),1
                NCT = NCT + 1
              end do
            end do
          end do
        end do
      end do
*
      if (NTEST.ge.2) then
        write(6,'(/i14,a)') NCT,
     &             ' 2-el. integrals need to be allocated for DIRAC!'
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Calculate total number of integrals on record
* of a symmetry reduced list (r/i, particle symmetry)
*
      subroutine numints_red(IGLSUM,NSMBLK,SYDI,NIRR_PN,IPRINT)
*
*------------------------------------
*  Timo Fleig, October 2000
*
*   modified for scalar relativistic runs with LUCITA
*   Nov. 2000
*------------------------------------
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "multd2h.inc"
#include "files.inc"
#include "clunit.inc"
#include "parluci.h"
*
      integer SYDI
      dimension SYDI(NIRR_PN)
      dimension IDEKS(MXPORB),IQISYM(MXPIRR)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.1) then
        write(LUWRT,*)
        write(LUWRT,*) '=========================='
        write(LUWRT,*) ' numints_red speaking:    '
        write(LUWRT,*) '=========================='
      end if

!     initialize IGLSUM
      IGLSUM = 0
!
*
* Set up triangular array
      II = 0
      do I=1,MXPORB,1
        IDEKS(I) = II
        II = II + I
      end do
*
* Initialize arrays
      do IS=1,NIRR_PN,1
        IQISYM(IS) = 0
      end do
*
* Symmetry packing labels
      NOCC = 0
      do IS=1,NIRR_PN,1
        if (SYDI(IS).eq.0) goto 10
        NOCC = NOCC + 1
        IQISYM(NOCC) = IS
   10   continue
      end do
      if (NTEST.ge.5) then
        write(LUWRT,*) 'Symmetry packing labels:'
        call iwrtma(IQISYM,NOCC,1,NOCC,1)
      end if
*
* Number of symmetry blocks and total number of integrals
      INIT = 1
      do ISP=1,NIRR_PN,1
        if (NTEST.ge.10) write(LUWRT,*) 'ISP =',ISP
        NSP=IQISYM(ISP)
        if (NSP.eq.0)then
          write(LUWRT,*) 'Number of symm for P is  ',NSP
          write(LUWRT,*) 'Skipping two-el. int. read.'
          goto 333
        end if
        if (NTEST.ge.10) write(LUWRT,*) 'NSP =',NSP
        NOP=SYDI(NSP)
        do ISQ=1,ISP,1
          if (NTEST.ge.10) write(LUWRT,*) 'ISQ =',ISQ
          NSQ = IQISYM(ISQ)
          NOQ = SYDI(NSQ)
          ISPQ = MULTD2H(NSP,NSQ)
          if (NTEST.ge.10) write(LUWRT,*) 'NSQ =',NSQ
          ntest = 00
          do ISR=1,ISP,1
            if (NTEST.ge.10) write(LUWRT,*) 'ISR =',ISR
            NSR = IQISYM(ISR)
            NOR = SYDI(NSR)
            ISPQR = MULTD2H(ISPQ,NSR)
            ISSM=ISR
            if (ISR.EQ.ISP) ISSM=ISQ
            if (NTEST.ge.10) write(LUWRT,*) 'NSR =',NSR
            do 40 ISS=1,ISSM,1
              if (NTEST.ge.10) write(LUWRT,*) 'ISS =',ISS
              NSS=IQISYM(ISS)
              NOS=SYDI(NSS)
              if (NTEST.ge.10) then
                write(LUWRT,'(A,4I6)') 
     &          'NOP,NOQ,NOR,NOS',NOP,NOQ,NOR,NOS
                write(LUWRT,'(A,2I3)') 'ISPQ,ISPQR',ISPQ,ISPQR
              end if
              if (NTEST.ge.10) write(LUWRT,*) 'NSS =',NSS
              if(NSS.NE.ISPQR) goto 40
              if(NOS*NOR*NOQ*NOP.eq.0) goto 40
*
*  NUMBER OF INTEGRALS IN THIS SYMMETRYBLOCK
*
              if (NSP.EQ.NSQ) then
*
*  TYPES <PP|PP> AND <PP|QQ>
*
                NPQ = IDEKS(NOP+1)
                NRS = IDEKS(NOR+1)
              else
*
*  TYPES <PQ|PQ> AND <PQ|RS>
*
                NPQ = NOP*NOQ
                NRS = NOR*NOS
              end if
              if (NSP.EQ.NSR) then
*
*  TYPES <PP|PP> AND <PQ|PQ>
*
                NIGL = NPQ*(NPQ+1)/2
              else
*
*  TYPES <PP|QQ> AND <PQ|RS>
*
                NIGL = NPQ*NRS
              end if
              INIT = INIT + 1
              IGLSUM = IGLSUM + NIGL
              if (NTEST.ge.1) then
                 write(LUWRT,601) NSP,NSQ,NSR,NSS,NIGL
              end if
40          continue
          end do
        end do
* Number of non-vanishing symmetry blocks
        NSMBLK = INIT - 1
333     continue
      end do
601   FORMAT(' NUMBER OF INTEGRALS IN SYMMETRYBLOCK '
     &             ,4I2,' IS : ',I10)
*
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PUTINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM)
*
* Put integrals in permanent integral list
*
* Jeppe Olsen, Jan. 1999
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*. Specific input
      DIMENSION XINT(*)
*
      CALL QENTER('PUTIN')
*. Offset and number of integrals
*
      IF(ITP.EQ.0) THEN
        NI = NTOOBS(ISM)
      ELSE
        NI = NOBPTS(ITP,ISM)
      END IF
*
      IOFF = IBSO(ISM)
      DO IITP = 1, ITP -1
        IOFF = IOFF + NOBPTS(IITP,ISM)
      END DO
*
      IF(JTP.EQ.0) THEN
        NJ = NTOOBS(JSM)
      ELSE
        NJ = NOBPTS(JTP,JSM)
      END IF
*
      JOFF = IBSO(JSM)
      DO JJTP = 1, JTP -1
        JOFF = JOFF + NOBPTS(JJTP,JSM)
      END DO
*
      IF(KTP.EQ.0) THEN
        NK = NTOOBS(KSM)
      ELSE
        NK = NOBPTS(KTP,KSM)
      END IF
*
      KOFF = IBSO(KSM)
      DO KKTP = 1, KTP -1
        KOFF = KOFF + NOBPTS(KKTP,KSM)
      END DO
*
      IF(LTP.EQ.0) THEN
        NL = NTOOBS(LSM)
      ELSE
        NL = NOBPTS(LTP,LSM)
      END IF
*
      LOFF = IBSO(LSM)
      DO LLTP = 1, LTP -1
        LOFF = LOFF + NOBPTS(LLTP,LSM)
      END DO
*
      INT_IN = 0
      DO LOB = LOFF,LOFF+NL-1
       DO KOB = KOFF,KOFF+NK-1
        DO JOB = JOFF,JOFF+NJ-1
         DO IOB = IOFF,IOFF+NI-1
C?         WRITE(6,*) ' IOB, JOB, KOB, LOB', IOB,JOB,KOB,LOB
           INT_OUT = I2EAD(IOB,JOB,KOB,LOB)
           INT_IN = INT_IN + 1
C?         WRITE(6,*) ' INT_OUT, INT_IN ', INT_OUT,INT_IN
C?         WRITE(6,*) ' KINT2-1+INT_OUT = ',KINT2-1+INT_OUT
           WORK(KINT2-1+INT_OUT) = XINT(INT_IN)
         END DO
        END DO
       END DO
      END DO
*
      CALL QEXIT('PUTIN')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PUTMOAO(CMOAO)
*
* SAVE   MOAO matrix CMOAO on LUMOUT
*
* A sunny day in April 96
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "crun.inc"
*
      IF(INTIMP.EQ.1) THEN
*. MOLCAS environment
        WRITE(6,*) ' PUTMOAO : MOLCAS environment'
        CALL PUTMOAO_MOLCAS(CMOAO,LUMOUT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PUTMOAO_MOLCAS(CMOAO,LU)
*
* WRITE MOAO matrix CMOAO on file LU in MOLCAS LUMORB format

*
* GETOBS assumed called to define /MOLOBS/
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 TITLEMO
      COMMON/MOLOBS/
     & IOList(64),iToc(64),nBas(8),nOrb(8),nFro(8),nDel(8),
     & Nsym
*
      LOCC = 0
*. Full NBAS X NBAS assumed
* Routine with the same name exists in DIRAC.
C     CALL WRVEC('CIAORB',LU,NSYM,NBAS,NBAS,CMOAO,OCC,LOCC,
C    &           ' MO orbitals obtained from LUCIA ')
      WRITE(6,*) ' Mo coefficients written to ', LU
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine rdcom
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Read integrals in DIRAC format
* routines from LUCIAREL interface
*
      subroutine rdint_dirac(KINT1,KINT2,DOUGRP,IPRINT)
*
      implicit real*8(A-H,O-Z)
*
      character*3 DOUGRP
      character*8 DIRACTWO
      character*7 DIRACONE
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "clunit.inc"
#include "parluci.h"
*
      dimension NBASD(MXPIRR),NSTR(2)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      if (myproc.eq.master) then
        write(6,*)
        write(6,*) 'Using spinfree Dirac Hamiltonian!'
      else
        write(12,*)
        write(12,*) 'Using spinfree Dirac Hamiltonian!'
      end if
*     
      ISPINFREE = 1
      ISYDISUM = 0
      DIRACONE = 'MRCONEE'
      DIRACTWO = 'MDCINT  '
      do ISYM = 1,NIRREP,1
        ISYDISUM = ISYDISUM + NTOOBS(ISYM)
      end do
*
* Interface to fully relativistic environment necessitates
* this classification:
*     Point group if of quaternionic type:
      if (DOUGRP.eq.'C1 '.or.DOUGRP.eq.'Ci ') then
        IPGTYP = 1
        if (DOUGRP.eq.'Ci ') then
          INVERSM = 1
        else
          INVERSM = 0
        end if
      else if (DOUGRP.eq.'C2 '.or.DOUGRP.eq.'Cs '.or.
     &         DOUGRP.eq.'C2h') then
        IPGTYP = 2
        if (DOUGRP.eq.'C2h') then
          INVERSM = 1
        else
          INVERSM = 0
        end if
      else if (DOUGRP.eq.'D2 '.or.DOUGRP.eq.'C2v'.or.
     &         DOUGRP.eq.'D2h') then
        IPGTYP = 3
        if (DOUGRP.eq.'D2h') then
          INVERSM = 1
        else
          INVERSM = 0
        end if
      end if
*
* Static memory
      call memman(KBOSYM,2*ISYDISUM,'ADDS  ',2,'IBOSUM')
*
* Set mark
      call memman(IDUM,IDUM,'MARK  ',IDUM,'RDIRAC')
*
* Local memory
      NUM_ONEI = (2*ISYDISUM)**2
      call memman(KDUMAR,ISYDISUM*2,'ADDL  ',2,'DUMARR')
      call memman(KIDUMAR,ISYDISUM*2,'ADDL  ',2,'IDUMAR')
      call memman(KDIONEE,NUM_ONEI,'ADDL  ',2,'IDIONE')
*
* This is a VERY preliminary solution for running bigger
* calculations. Partial symmetry reduction of allocated space.
!     NTEST_save = NTEST
!     NTEST  = 2000

      call nsymint_dc(NUM_TWEI,NTOOBS,NIRREP,ISYDISUM,IPRINT)
      if (NTEST.ge.1)
     & write(6,*) 'Allocating space for ',NUM_TWEI,' 2-el. ints'
      call memman(KDITWOE,NUM_TWEI,'ADDL  ',2,'IDITWO')
      call dzero(WORK(KDITWOE),NUM_TWEI)
      call memman(KINVREO,ISYDISUM,'ADDL  ',2,'INVREO')
      call memman(KKR,(2*ISYDISUM+1),'ADDL  ',2,'IKR   ')
      call memman(KINDK,num_twei,'ADDL  ',1,'INDK  ')
      call memman(KINDL,num_twei,'ADDL  ',1,'INDL  ')
      call memman(KIOFF12,(2*ISYDISUM)**2,'ADDL  ',2,'IOFF12')
      call memman(KNUMINT,(2*ISYDISUM)**2,'ADDL  ',2,'NUMINT')
      call findhigh(NTOOBS(1),NIRREP,IMXDIM)
      call memman(KIREO,IMXDIM*NIRREP,'ADDL  ',2,'IIREO')
      IMXDDIM = IMXDIM*2
      call memman(KPICHE2,(2*IMXDDIM)**2,'ADDL  ',2,'IPIHE2')
*
      call rdone_dirac(WORK(KDIONEE),LUOINT,DIRACONE,MORB,
     &                 ISYDISUM,NBASD,ISPINFREE,WORK(KBOSYM),
     &                 WORK(KIDUMAR),WORK(KDUMAR),NIRREP,
     &                 NGSSH,NGAS,NSTR,
     &                 IPGTYP,INVERSM,INVERSM+1,NTEST)
*
      call rdtwo_dirac(WORK(KDITWOE),
     &                 LUTINT,DIRACTWO,
     &                 ISYDISUM,ISPINFREE,MORB,
     &                 WORK(KKR),WORK(KIOFF12),WORK(KNUMINT),
     &                 WORK(KINDK),WORK(KINDL),NTEST)
*
      IHAM12 = 2
      ISYDISUM2 = ISYDISUM*2
      call dist_ints_dirac2(WORK(KPICHE2),WORK(KDIONEE),
     &                     WORK(KINT1),WORK(KDITWOE),
     &                     WORK(KKR),WORK(KIOFF12),WORK(KNUMINT),
     &                     WORK(KINDK),WORK(KINDL),
     &                     WORK(KINT2),
     &                     WORK(KBOSYM),WORK(KIREO),
     &                     WORK(KINVREO),NTOOBS,NSTR,
     &                     ISYDISUM,ISYDISUM2,
     &                     INVERSM,
     &                     IHAM12,ISPINFREE,IMXDIM,NIRREP,NTEST)
*
!     NTEST = NTEST_save

      call memman(IDUM,IDUM,'FLUSM ',IDUM,'RDIRAC')
*
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Input routine for LUCIA under MOLCAS 4
*  A simplified MOLCAS-type input is read and converted to
*  the usual LUCIA requirements.
*
*  Note concerning EXPERT mode:
*   If you want to use it, the first 2 lines of your input for
*   MOLUC-LUCIA (LUCITA) should be
*       &LUCITA &END
*      EXPERT
*      * &LUCIA &END
*       ...    (first line of regular input for LUCIA)
*       ...    (regular input for LUCIA)
*       ...
*       END OF INPUT
*   And ONE blank line after 'END OF INPUT'. This is essential,
*   otherwise the program will crash.
*
*  Feb-Mar 1999,   Timo Fleig
*  Feb     2000,   Timo Fleig
*
      subroutine rdmoluc
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "SysDef.inc"
#include "luci.inc"
#include "units.inc"
*
      parameter (NCMD = 21)
      dimension IMOKW(NCMD)
      character*3 SCR3
      character*4 COMMAND,CMD(NCMD)
      character*6 SCR6,MOLUCENV
      character*72 LINE,CAREXP(1000),
     &             CARD(MXPNGAS),CARD2(MXPNGAS,MXPICI),CARD3(5)
*
      common/MLCENV/MOLUCENV
*
*  CMD(1) = 'TITL', CMD(2) = 'CITP' ... etc.
      data CMD /'TITL','INWF','CITP','NROO','SSYM',
     &          'NACT','MULT','PRIN','PRLO','SIZE',
     &          'INAC','GASS','GSSP','FRMO','RAS1',
     &          'RAS2','RAS3','DENS','ENVI','EXPE',
     &          'END '/
*
* Set logical units
      call setunits
*
      do ICMD = 1,NCMD,1
         IMOKW(ICMD) = 0
      end do
*
      NTIT = 0
*
*
      call qenter('rdmoluc')
*
*---  Read input from standard input ----------------------------------*
*  Read Namelist
      call rdnlst(5,'LUCITA')
10    read(*,'(A72)',End=991) LINE
*  Put input line to very left
      COMMAND = LINE(1:4)
*  Convert possible lower case characters to upper case
      call upcase(COMMAND)
*  I have read a comment line
      if ( COMMAND(1:1).eq.'*' ) goto 10
      JCMD=0
*  Compare read commands with hardwired character strings
*  If found, store in JCMD
20    continue
      do ICMD=1,NCMD
         If ( COMMAND.eq.CMD(ICMD) ) then
           JCMD = ICMD
           IMOKW(JCMD) = 1
         end if
      End Do
*  If not found, call error message
      if ( JCMD.eq.0 ) Call ErrMsg('rdmoluc',001)
*  Go to predefined number of this keyword (if found)
      Goto ( 100,  200,  300,  400,  500,
     &       600,  700,  800,  900,  1000,
     &       1100, 1200, 1300, 1400, 1500,
     &       1600, 1700, 1800, 1900, 2000,
     &       3000 ) JCMD
*
*---  process TITL command --------------------------------------------*
 100  continue
      read(*,'(A72)',End=991) LINE
      COMMAND = LINE(1:4)
      call UpCase(COMMAND)
      if ( COMMAND(1:1).eq.'*' ) goto 100
*  Check whether LINE is keyword:
      JCMD=0
      do ICMD=1,NCMD
         if ( COMMAND.eq.CMD(ICMD) ) JCMD = ICMD
      End Do
*  I have found the next keyword. Process it:
      If ( JCMD.ne.0 ) Goto 20
*  Still title line. Store it:
      NTIT = NTIT + 1
      if ( NTIT.le.MXTIT ) TITLE(NTIT) = LINE
      goto 100
*
*---  process INWF command --------------------------------------------*
*     Type of initial orbitals to be used
 200  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 200
      SCR6 = LINE(1:6)
      read(SCR6,'(A6)',Err=992) WAFFCT
      call upcase(WAFFCT)
      goto 10
*
*---  process CITP command --------------------------------------------*
*     Type of CI calculation to be performed
 300  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 300
      SCR6 = LINE(1:6)
      read(SCR6,'(A6)',Err=992) TYPE
      call upcase(TYPE)
      goto 10
*
*---  process NROO command --------------------------------------------*
 400  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 400
      read(LINE,*,Err=992) NROOT
      goto 10
*
*---  process SSYM command --------------------------------------------*
 500  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 500
      read(LINE,*,Err=992) IRREP
      goto 10
*
*---  process NACT command --------------------------------------------*
*     Number of active electrons
 600  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 600
      read(LINE,*,Err=992) NACTEL
      goto 10
*
*---  process MULT command --------------------------------------------*
*     Spin multiplicity
 700  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 700
      read(LINE,*,Err=992) MULTIP
      goto 10
*
*---  process PRIN command --------------------------------------------*
*     Global print flag for subsequent LUCIA run
 800  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 800
      SCR3 = LINE(1:3)
      read(SCR3,'(A3)',Err=992) PRILUC
      call upcase(PRILUC)
      goto 10
*
*---  process PRLO command --------------------------------------------*
*     Local print flag for moluc
 900  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 900
      read(LINE,*,Err=992) IPLOCAL
      goto 10
*
*---  process SIZE command --------------------------------------------*
*     Set approximate size of CI calculation
1000  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 800
      SCR3 = LINE(1:3)
      read(SCR3,'(A3)',Err=992) CALSIZ
      call upcase(CALSIZ)
      goto 10
*
*---  process INAC command --------------------------------------------*
*     Number of orbitals in the inactive space
*     (can be used with FCI, SDCI, and SDTQ;
*      not with GASCI;
*      must with RASCI)
1100  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1100
      read(LINE,'(A72)',Err=992) CARD3(2)
      goto 10
*
*---  process GASS command --------------------------------------------*
*     Number of orbitals per GAS per symmetry
1200  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1200
      read(LINE,*,Err=992) NGAS
      do IGAS=1,NGAS,1
        read(*,'(A72)',End=991) LINE
        call upcase(LINE)
        read(LINE,'(A72)',Err=992) CARD(IGAS)
      end do
      goto 10
*
*---  process GSSP command --------------------------------------------*
*     Cumulative min. and max. number of electrons per GAS
1300  continue
      if (IMOKW(12).eq.0) then
        write(6,*) 'GASS has to be specified before GSSP.'
        write(6,*) 'I do not know what to do.'
        Call Abend2('Quitting.')
      end if
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1300
      read(LINE,*,Err=992) NCISPC
      do ISPC=1,NCISPC,1
        do IGAS=1,NGAS,1
          read(*,'(A72)',End=991) LINE
          read(LINE,'(A72)',Err=992) CARD2(IGAS,ISPC)
        end do
      end do
      goto 10
*
*---  process FRMO command --------------------------------------------*
*     Number of frozen orbitals in MOTRA transformation
1400  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1400
      read(LINE,'(A72)',Err=992) CARD3(1)
      goto 10
*
*---  process RAS1 command --------------------------------------------*
*     Number of orbitals in RAS1 per symmetry
*     Maximum number of holes in RAS1
1500  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1500
      read(LINE,'(A72)',Err=992) CARD3(3)
      read(*,'(A72)',End=991) LINE
      read(LINE,*,Err=992) MXHOL1
      goto 10
*
*---  process RAS2 command --------------------------------------------*
*     Number of orbitals in RAS2 per symmetry
1600  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1600
      read(LINE,'(A72)',Err=992) CARD3(4)
      goto 10
*
*---  process RAS3 command --------------------------------------------*
*     Number of orbitals in RAS3 per symmetry
*     Maximum number of electrons in RAS3
1700  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1700
      read(LINE,'(A72)',Err=992) CARD3(5)
      read(*,'(A72)',End=991) LINE
      read(LINE,*,Err=992) MXELR3
      goto 10
*
*---  process DENS command --------------------------------------------*
*     Calculation of one- and/or two-particle density matrix
1800  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1800
      read(LINE,*,Err=992) IDEMOL
      goto 10
*
*---  process ENVI command --------------------------------------------*
*     Interface to DIRAC spinfree formalism (scalar relativistic)
1900  continue
      read(*,'(A72)',End=991) LINE
      if ( LINE(1:1).eq.'*' ) goto 1900
      read(LINE,'(A6)',Err=992) MOLUCENV
      goto 10
*
*---  process EXPE command --------------------------------------------*
*     Expert mode for LUCIA; skip all further input; no output
2000  continue
      goto 6666
*
3000  Continue
*----------------------------------------------------------------------*
*  Insert defaults and print (error) messages where appropriate
*----------------------------------------------------------------------*
*
*  Code for keyword status vector:
*
*     IMOKW(N) = 0           keyword not provided in input,
*                            default or error
*     IMOKW(N) = 1           keyword given in input
*     IMOKW(N) = 2           keyword has been set to default value
*     IMOKW(N) = 3           keyword ignored, default is used
*
*  1: Default Title
*
      if (IMOKW(1).eq.0) then
         TITLE(1) = ' Running LUCIA under MOLCAS (DIRAC interface)'
         TITLE(2) = ' Version of Nov / 2000 '
         TITLE(3) = ' (You should AT LEAST supply a title !) '
         IMOKW(1) = 2
      end if
*
*  2: Type of initial wave function (no default)
*
      if (IMOKW(2).eq.0) then
         write(6,*)
     & ' Keyword for type of initial wave function missing. '
         write(6,*) ' This keyword is mandatory. '
         Call Abend2('Input error')
      else
         if (WAFFCT.ne.'HF_SCF'.and.WAFFCT.ne.'RASSCF') then
            write(6,*)
     &  ' Type of initial wave function not implemented. '
            write(6,*) ' You have chosen: ',WAFFCT
            write(6,*) ' Allowed types: '
            write(6,*) ' HF_SCF, RASSCF '
            Call Abend2('Input error')
         end if
      end if
*
*  3: CI type for LUCIA (no default)
*
      if (IMOKW(3).eq.0) then
         write(6,*) ' Keyword for type of CI calculation missing. '
         write(6,*) ' This keyword is mandatory. '
         Call Abend2('Input error')
      else
         if (TYPE.ne.'FCI   '.and.TYPE.ne.'SDCI  '.and.
     &       TYPE.ne.'GASCI '.and.TYPE.ne.'RASCI '.and.
     &       TYPE.ne.'SDTQ  ') then
            write(6,*) ' Type of CI calculation not specified.'
            write(6,*) ' You have chosen: ',TYPE
            write(6,*) ' Allowed types: '
            write(6,*) ' FCI, SDCI, GASCI, RASCI, SDTQ '
            Call Abend2('Input error')
         end if
      end if
*
*  4: Default for number of roots to be obtained
*
      if (IMOKW(4).eq.0) then
         NROOT = 1
         IMOKW(4) = 2
      end if
*
*  5: Default for state symmetry in point group
*
      if (IMOKW(5).eq.0) then
         IRREP = 1
         IMOKW(5) = 2
      end if
*
*  6: Default number of active electrons
*
      if (IMOKW(6).eq.0) then
         if (TYPE.eq.'GASCI '.or.TYPE.eq.'RASCI ') then
            write(6,*) 'Number of active electrons NACTEL'
            write(6,*) 'has to be specified in this type of'
            write(6,*) 'calculation. Quitting.'
            Call Abend()
         else
            NACTEL = -1
            IMOKW(6) = 2
         end if
      else if (IMOKW(6).eq.1) then
         if (TYPE.eq.'FCI   '.or.TYPE.eq.'SDCI  '.or.
     &       TYPE.eq.'SDTQ  ') then
            write(6,*) 'In CI types FCI, SDCI, and SDTQ'
            write(6,*) 'I am using default numbers of'
            write(6,*) 'active electrons. Your input will'
            write(6,*) 'be ignored. Program will continue.'
            write(6,*)
            write(6,*) 'Use GASCI or RASCI for ionic compounds'
            write(6,*) 'or provide the orbitals accordingly.'
            NACTEL = -1
            IMOKW(6) = 2
         end if
      end if
*
*  7: Spin multiplicity (no default)
*
      if (IMOKW(7).eq.1) then
         if (IMOKW(6).eq.1) then
            call evenodd(IEONAC,NACTEL)
            call evenodd(IEOMUL,MULTIP)
            if (IEONAC.eq.2.and.IEOMUL.eq.2) then
               write(6,*) 'Illegal spin multiplicity given.'
               write(6,*) 'Read a book about fermions.'
               Call Abend2('quitting')
            else if (IEONAC.eq.1.and.IEOMUL.eq.1) then
               write(6,*) 'Illegal spin multiplicity given.'
               write(6,*) 'Read a book about fermions.'
               Call Abend2('quitting')
            else if (IEONAC.eq.1.and.IEOMUL.eq.2) then
               if (MULTIP.lt.2.or.MULTIP.gt.(NACTEL+1)) then
                  write(6,*) 'Illegal spin multiplicity given.'
                  write(6,*) 'Compare with number of active'
                  write(6,*) 'electrons.'
                  Call Abend2('quitting')
               end if
            else if (IEONAC.eq.2.and.IEOMUL.eq.1) then
               if (MULTIP.lt.1.or.MULTIP.gt.(NACTEL+1)) then
                  write(6,*) 'Illegal spin multiplicity given.'
                  write(6,*) 'Compare with number of active'
                  write(6,*) 'electrons.'
                  Call Abend2('quitting')
               end if
            end if
         end if
      else if (IMOKW(7).eq.0) then
         write(6,*) 'Spin multiplicity is a MANDATORY keyword.'
         write(6,*) 'Specify and restart.'
         Call Abend2('Quitting.')
      end if
*
*  8: Default for global print parameter
*
      if (IMOKW(8).eq.0) then
         PRILUC = 'LOW'
         IMOKW(8) = 2
      end if
*
*  9: Default for local print parameter
*
      if (IMOKW(9).eq.0) then
         IPLOCAL = 0
         IMOKW(9) = 2
      end if
*
* 10: Default for approximate size of calculation
*
      if (IMOKW(10).eq.0) then
         CALSIZ = 'NOR'
         IMAXCIITER = 100
         IMOKW(10) = 2
      else
         if (CALSIZ.ne.'NOR'.and.CALSIZ.ne.'LAR'.and.
     &       CALSIZ.ne.'HUG') then
           write(6,'(A,A3)') 'SIZE specified to ',CALSIZ
           write(6,*) 'This is an unknown type.'
           Call Abend2('Quitting.')
         end if
         if (CALSIZ.eq.'NOR'.or.CALSIZ.eq.'LAR') IMAXCIITER = 100
         if (CALSIZ.eq.'HUG') IMAXCIITER = 50
      end if
*
* 11: Orbital distribution in inactive space (no defaults if RASCI)
*
      if (IMOKW(11).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'Number of inactive orbitals per sym'
            write(6,*) 'has to be specified. This is mandatory'
            write(6,*) 'in a RASCI calculation.'
            Call Abend2('Quitting.')
         else
            INACOB = 0
            IMOKW(11) = 2
         end if
      else if (IMOKW(11).eq.1) then
         if (TYPE.eq.'GASCI ') then
            write(6,*) 'INACtive orbitals not allowed in GASCI.'
            write(6,*) 'Use GAS 1 for these orbitals and restart.'
            Call Abend2('Quitting.')
         end if
         INACOB = 1
      end if
*
* 12: Orbital distribution in GAS spaces (no defaults if GASCI)
*
      if (IMOKW(12).eq.0) then
         if (TYPE.eq.'GASCI ') then
            write(6,*) 'GASCI type requires GASS to be specified.'
            write(6,*) 'Else, I do not know what to do.'
            Call Abend2('Quitting.')
         end if
         IMOKW(12) = 2
      end if
      if (IMOKW(12).eq.1.and.TYPE.ne.'GASCI ') then
         write(6,*) 'GASS specified and not doing GASCI.'
         write(6,*) 'This input will be ignored.'
         write(6,*) 'Program will continue with ',TYPE,'.'
         IMOKW(12) = 3
      end if
*
* 13: Cumulative min. and max. numbers of electrons in GAS spaces
*
      if (IMOKW(13).eq.0) then
         if (TYPE.eq.'GASCI ') then
            write(6,*) 'GSSP has to be specified in GASCI calculation.'
            write(6,*) 'No defaults in this type of CI.'
            Call Abend2('Quitting.')
         end if
         IMOKW(13) = 2
      end if
      if (IMOKW(13).eq.1.and.TYPE.ne.'GASCI ') then
         write(6,*) 'GSSP specified and not doing GASCI.'
         write(6,*) 'This input will be ignored.'
         write(6,*) 'Program will continue with ',TYPE,'.'
         write(6,*) 'Defaults will be set.'
         IMOKW(13) = 3
      end if
*
* 14: FRMO (frozen orbitals from MOTRA) specification; no default
*     can be used with any CI or WF type
*
      if (IMOKW(14).eq.0) then
         IMOKW(14) = 2
      end if
*
* 15: RAS1 specification; no default
*
      if (IMOKW(15).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'RAS1 keyword missing. Add and restart.'
            Call Abend2('Quitting.')
         end if
         IMOKW(15) = 2
      else if (IMOKW(15).eq.1) then
         if (TYPE.ne.'RASCI ') then
            write(6,*) 'Keyword RAS1 not compatible with TYPE ',TYPE
            write(6,*) 'I will ignore your this input.'
            write(6,*) 'Program will continue.'
            IMOKW(15) = 3
         end if
      end if
*
* 16: RAS2 specification; no default
*
      if (IMOKW(16).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'RAS2 keyword missing. Add and restart.'
            Call Abend2('Quitting.')
         end if
         IMOKW(16) = 2
      else if (IMOKW(16).eq.1) then
         if (TYPE.ne.'RASCI ') then
            write(6,*) 'Keyword RAS2 not compatible with TYPE ',TYPE
            write(6,*) 'I will ignore your this input.'
            write(6,*) 'Program will continue.'
            IMOKW(16) = 3
         end if
      end if
*
* 17: RAS3 specification; no default
*
      if (IMOKW(17).eq.0) then
         if (TYPE.eq.'RASCI ') then
            write(6,*) 'RAS3 keyword missing. Add and restart.'
            Call Abend2('Quitting.')
         end if
         IMOKW(17) = 2
      else if (IMOKW(17).eq.1) then
         if (TYPE.ne.'RASCI ') then
            write(6,*) 'Keyword RAS3 not compatible with TYPE ',TYPE
            write(6,*) 'I will ignore your this input.'
            write(6,*) 'Program will continue.'
            IMOKW(17) = 3
         end if
      end if
*
* 18: Default setting for DENSity matrices
*
      if (IMOKW(18).eq.0) then
         IDEMOL = 0
         IMOKW(18) = 2
      end if
*
* 19: Default setting for MOLUCENVironment variable
*
      if (IMOKW(19).eq.0) then
         MOLUCENV = 'MOLCAS'
         IMOKW(19) = 2
      end if
*
*-----------------------------------------------------------------*
*  Print section
*-----------------------------------------------------------------*
*
*  1: Title
*
      write(*,*)
c     call XFLUSH(6)
      write(*,'(6X,120A1)') ('*',I=1,120)
c     call XFLUSH(6)
      write(*,'(6X,120A1)') '*',(' ',I=1,118),'*'
c     call XFLUSH(6)
      write(*,'(6X,57A1,A6,57A1)')
     &        '*',(' ',I=1,56),'Title:',(' ',I=1,56),'*'
c     call XFLUSH(6)
      do I=1,NTIT
         write(*,'(6X,24A1,A72,24A1)')
     &        '*',(' ',J=1,23),Title(I),(' ',J=1,23),'*'
c     call XFLUSH(6)
      end do
      write(*,'(6X,120A1)') '*',(' ',I=1,118),'*'
c     call XFLUSH(6)
      write(*,'(6X,120A1)') ('*',I=1,120)
c     call XFLUSH(6)
      write(*,*)
*
*  2: Initial wave function
*
      write(6,'(A42,A6)')
     &  ' Orbitals as initial wave function .... ',WAFFCT
      write(6,*)
*
*  3: Type of CI calculation
*
      write(6,'(A42,A6)')
     &  ' Type of calculation .................. ',TYPE
      write(6,*)
*
*  4: Number of roots to be treated
*
      write(6,'(A42,I3)')
     &  ' Number of roots to be obtained ....... ',NROOT
      write(6,*)
*
*  5: State symmetry
*
      write(6,'(A42,I3)')
     &  ' Calculation carried out in irrep ..... ',IRREP
      write(6,*)
*
*  6: Number of active electrons
*
      if (IMOKW(6).eq.1) then
        write(6,'(A42,I3)')
     &  ' Number of active electrons ........... ',NACTEL
        write(6,*)
      end if
*
*  7: Spin multiplicity
*
      write(6,'(A42,I3)')
     &  ' Spin multiplicity .................... ',MULTIP
      write(6,*)
*
*  8: LUCIA global print parameter
*
      write(6,'(A42,A3)')
     &  ' Global print level is ................ ',PRILUC
      write(6,*)
*
*  9: MOLUC local print parameter
*
      write(6,'(A42,1I3)')
     &  ' Local print level is ................. ',IPLOCAL
      write(6,*)
*
* 10: SIZE of CI calculation
*
      write(6,'(A42,A3)')
     &  ' Approximate size of CI calculation ... ',CALSIZ
      write(6,*)
*
* 11: INAC specification of inactive orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 12: GASS specification of orbitals per GAS per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 13: GSSP specification of electrons in GAS Spaces
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 14: FRMO specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 15: RAS1 specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 16: RAS2 specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 17: RAS3 specification of orbitals per symmetry
*
*     No printing. Will be directed to luciwrt after resolution
*
*
* 18: DENSity matrix level
*
      if (IDEMOL.ge.1) then
        write(6,'(A42,1I3)')
     &    ' Calculating density matrices at level. ',IDEMOL
        write(6,*)
      end if
*
* 19: DIRAc environment setting
*
      if (MOLUCENV(1:5).eq.'DIRAC') then
        write(6,'(A42,A6)')
     &  ' Moluc environment variable set to .... ',MOLUCENV
        write(6,*)
        write(6,'(A50)')
     &  ' This interfacing might be done in the future. '
        write(6,'(A50)')
     &  ' For the moment, run in EXPERT mode under DIRAC. '
        call Abend2('Quitting.')
        write(6,*)
      end if
*
*-----------------------------------------------------------------*
*  Call conversion routine for creating LUCIA input
*-----------------------------------------------------------------*
      call mol2luc(NTIT,NGAS,NCISPC,IPLOCAL,
     &             CARD,CARD2,CARD3,1.0d-8)
*
6666  continue
      if (IMOKW(20).eq.1) then
*-----------------------------------------------------------------*
*  Expert mode
*-----------------------------------------------------------------*
         NLIN = 1
5555     continue
         read(*,'(A72)',End=991) CAREXP(NLIN)
         NLIN = NLIN + 1
         if (CAREXP(NLIN-1).ne.
     &'                                              ') goto 5555
*modifications ???* sk 
         open(unit=LUINP,file='LUCIAIN',
     &        status='NEW',form='FORMATTED')
           do I=1,NLIN-2,1
             write(LUINP,'(A72)') CAREXP(I)
           end do
         close (unit=LUINP)
      else
        do ICMD=1,NCMD-2,1
           if (IMOKW(ICMD).eq.0) then
              write(6,*) 'Problem with keyword input'
              write(6,*) 'IMOKW(',ICMD,') = ',IMOKW(ICMD)
              Call Abend2('stopping')
           end if
        end do
      end if
*
      Call qExit('rdmoluc')
      RETURN
991   Call StkA(0,COMMAND)
      Call ErrMsg('rdmoluc',002)
      Call Quit(20)
992   Call StkA(0,COMMAND)
      Call ErrMsg('rdmoluc',003)
      Call Quit(20)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
*  Import one-electron integrals from DIRAC environment
*
*  MPI adaption by Stefan Knecht, Feb. 2006
*
*  Master reads all the information in and braodcasts 
*  everything to all nodes in global_communicator
*
      subroutine rdone_dirac(ONEINT,LUONE,DIRACONE,MORB,NDIM,
     &                       NBAS,ISPINFREE,IBOSYM,IDUMAR,DUMAR,
     &                       NIRREP,NGSSH,NGAS,NSTR,
     &                       IPGTYP,INVERSM,NFSYM,IPRINT)
*
      use interface_to_mpi
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "parluci.h"
*
      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX
*
      character*7 DIRACONE
      character*14 DUMCHA(MXPIRR)
      logical BREIT
*
      dimension    NBAS(NIRREP),ICT(MXPIRR),
     &             NGSSH(MXPIRR,MXPNGAS),IORBF(2),NSTR(NFSYM)
      dimension    ONEINT((2*NDIM)**2) 
      real(8)   :: DUMAR(*)
      integer   :: IBOSYM(MORB), IDUMAR(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
*
      IZERO = 0
      call isetvc(ICT,IZERO,NIRREP)
*
*
      IF (MYPROC.EQ.MASTER) THEN
C
      open(LUONE,FILE=DIRACONE,FORM='UNFORMATTED')
      read(LUONE) MORB,BREIT,ECORE_DIR
      read(LUONE) IDUMMY,(DUMCHA(ID),ID=1,IDUMMY),
     &                   (IDUMAR(ID),ID=1,IDUMMY),
     &                   (NSTR(I),I=1,NFSYM)
      read(LUONE)
      read(LUONE)
      read(LUONE) (IDUMAR(ID),IDUMAR(ID),DUMAR(ID),ID=1,MORB,1),
     &            (IBOSYM(I),I=1,MORB,1)
*
      if (NTEST.ge.3) then
        write(6,*)
        write(6,*) 'Testing read MRCONEE'
        write(6,*) 'MORB = ',MORB
        write(6,*) 'BREIT = ',BREIT
        write(6,*) 'Boson symmetries of spinors:'
        do IB=1,MORB,1
          write(6,'(A11,I5,A3,I12)') 'Spinor no. ',IB,' : ',IBOSYM(IB)
        end do
        do IFS=1,NFSYM,1
          write(6,'(A,I3,A,I5)')
     &          'Dim. of ferm. sym. ',IFS,' : ',NSTR(IFS)
        end do
      end if
*
      END IF
C    /\ read in master
C
C     Distribute the information read from MRCONEE
C
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
         CALL interface_MPI_BCAST(MORB,1,MASTER,global_communicator)
         CALL interface_mpi_bcast_l0(BREIT,1,MASTER,
     &                                   global_communicator)
         CALL interface_MPI_BCAST(ECORE_DIR,1,MASTER,
     &                            global_communicator)
         CALL interface_MPI_BCAST(NSTR,NFSYM,MASTER,
     &                            global_communicator)

         CALL interface_mpi_bcast_i1_work_f77(IBOSYM,MORB,MASTER,
     &                                        global_communicator)
      END IF
#endif

      IERROR = 0
      NORBSM = 0
      do ISPINOR=1,MORB,1
        ISYM = IBOSYM(ISPINOR) + 1
        ICT(ISYM) = ICT(ISYM) + 1
      end do
      do IRREP=1,NIRREP,1
        do IGAS=1,NGAS,1
          NORBSM = NORBSM + NGSSH(IRREP,IGAS)
        end do
        ICT(IRREP) = ICT(IRREP)/2
        if (ICT(IRREP).ne.NORBSM) IERROR = IERROR + 1
        NORBSM = 0
      end do
*
      if (2*NDIM.ne.MORB) then
        write(6,*) 'Orbital dimensions of '
        write(6,*) 'DIRAC and input do not match.'
        write(6,*) 'DIRAC : ',MORB,' spinors'
        write(6,*) 'Input : ',NDIM,' orbitals'
        call abend2('Quitting.')
      end if
      if (IERROR.ge.1) then
        write(6,*) IERROR,' orbital dimensions not matching.'
        write(6,*) 'Symmetry dimensions should be:'
        call iwrtma(ICT,1,NIRREP,1,MXPIRR)
        call abend2('Quitting.')
      end if
*
      IF (MYPROC.EQ.MASTER) THEN
* Read the integrals and close file
      read(LUONE) (ONEINT(I),TDUM,I=1,MORB*MORB)
      close(LUONE)
*

#if defined LUCI_DEBUG
      if (NTEST.ge.10) then
        do IND=1,MORB,1
          INOFF = (IND-1)*MORB
          write(6,*)
     &          IND,(ONEINT(INOFF+JND),JND=1,MORB,1)
        end do
      end if
#endif
      END IF
C     /\ read in master
C
C     Distribute the information read from MRCONEE
C
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
         CALL interface_MPI_BCAST(ONEINT,MORB*MORB,MASTER,
     &                            global_communicator)
         CALL interface_MPI_BCAST(TDUM,1,MASTER,
     &                            global_communicator)
         CALL interface_MPI_BCAST(ECORE_DIR,1,MASTER,
     &                            global_communicator)
      END IF
#endif
*
C     ...on here!
C
      if (myproc.eq.master) then
        write(6,'(1X,A,1F18.10)') 'DIRAC core energy = ',ECORE_DIR
      end if
C
      ECORE = ECORE_DIR
*
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine rdtwo_dirac(RKLR,
     &                       LUTWO,DIRACTWO,
     &                       NDIM,ISPINFREE,MORB,
     &                       KR,IOFF12,NUMINT,INDK,INDL,IPRINT)
C**********************************************************************
*     Read DIRAC-COULOMB two-electron integrals
*     (Later also BREIT integrals)
*     (GAUNT integrals can be included in COULOMB type ints.,
*      because they have the same symmetry.)
*
*     parallel calculation:
*     - broadcasting of general data from MASTER to nodes
*     - integrals are either read in by every node or distributed via
*       broadcast from MASTER to nodes
*     Stefan Knecht, Dec 08
*
C**********************************************************************
      use interface_to_mpi
      implicit real*8 (A-H,O-Z)
*
#include "parluci.h"
#if defined (VAR_MPI)
#include "infpar.h"
#endif
#include "files.inc"
*
      character*8 DIRACTWO
      character*10 DATEX,TIMEX*8
C     logical BREIT
      dimension RKLR(*)
      dimension KR(-NDIM:NDIM),INDK(*),INDL(*),
     &          IOFF12(MORB,MORB),NUMINT(MORB,MORB)
      INTEGER IDIST_IJKL, IDIST_IJKL_A(NMPROC)
      LOGICAL EX, IJKL_AV
C
      NTESTL = 00
      NTEST = max(NTESTL,IPRINT)
C
C
      NINT       = 0
      N2EFIL     = 0
      ISYMOFF    = 0
      IDIST_IJKL = 0
      IJKL_AV    = .TRUE.
      CALL IZERO(IDIST_IJKL_A,NMPROC)

!     initialize offset and integral counter arrays
      CALL IZERO(NUMINT,    MORB**2)
      CALL IZERO(IOFF12,    MORB**2)
C
C     check for integral file on local disks
      INQUIRE(FILE='MDCINT',EXIST=EX)
      IF( EX ) IDIST_IJKL = -1
      IDIST_IJKL_A(MYPROC+1) = IDIST_IJKL
#if defined (VAR_MPI)
      call interface_mpi_gather(idist_ijkl,1,idist_ijkl_a,1,master,
     &                          global_communicator)
      IF(MYPROC .eq. MASTER) THEN
        IDIST_IJKL_TEMP = IDIST_IJKL_A(MYPROC+1)
        DO I = 1, NMPROC
          IF( IDIST_IJKL_TEMP .ne. IDIST_IJKL_A(I)) IJKL_AV = .FALSE.
        END DO
      END IF
      CALL interface_mpi_bcast_l0(IJKL_AV,1,MASTER,
     &                                global_communicator)
#endif
C
C     check for availability of integral file
      IF( MYPROC .eq. MASTER .and. ( IDIST_IJKL .eq. 0 ) )THEN
        WRITE(LUWRT,*) ' *** ERROR in RDTWO_DIRAC: no integral file'//
     &                 ' present! ***'
        CALL QUIT( ' *** ERROR in RDTWO_DIRAC: no integral file 
     &  present! ***' )
      END IF
C
C     determine name of file (default is 'MDCINT  ')
   30 if (N2EFIL.ge.1.and.N2EFIL.le.9) then
        write(DIRACTWO,'(A7,I1)') 'MDCINT0',N2EFIL
      else if (N2EFIL.ge.10.and.N2EFIL.le.99) then
        write(DIRACTWO,'(A6,I2)') 'MDCINT',N2EFIL
      end if
C
      IF( MYPROC .eq. MASTER .or. IJKL_AV )THEN
C
C       Open file (using same logical unit for all files)
        open(LUTWO,FILE=DIRACTWO,FORM='UNFORMATTED')
        rewind (LUTWO)
      END IF
      if (N2EFIL.eq.0) then
        if( myproc .eq. master .or. IJKL_AV ) then
          read (LUTWO) datex,timex,nkr,
     &          (kr(i),kr(-i),i=1,nkr)
        end if
#if defined (VAR_MPI)
        if( nmproc .gt. 1 .and. (.not. IJKL_AV)) then
           call interface_mpi_bcast(nkr,1,master,
     &                              global_communicator)
           do j = 1, nkr
             call interface_mpi_bcast(kr(-j),1,master,
     &                                global_communicator)
             call interface_mpi_bcast(kr(j),1,master,
     &                                global_communicator)
           end do
        end if
#endif
      end if
*
      if(2*NKR.ne.MORB) then
        WRITE(LUWRT,*) 'MRCONEE and MDCINT are incompatible.'
        CALL QUIT( '*** ERROR in RDTWO_DIRAC: incompatible intgral 
     & files. ***')
      end if
      if (myproc.eq.master) then
        write(LUWRT,*) 'Integrals have been calculated on the'
        write(LUWRT,'(A12,1X,A)') '            ',DATEX,' at ',TIMEX
      end if
#if defined LUCI_DEBUG
      if (NTEST.ge.5) then
        write(LUWRT,*)
        write(LUWRT,*) 'Info from MDCINT file'
        write(LUWRT,*) 'NKR is (myproc)',NKR, myproc
        write(LUWRT,*) 'KR array:(myproc)', myproc
        do I=1,NKR,1
          write(LUWRT,'(2I8)') KR(-I),KR(I)
        end do
      end if
#endif
C
C     -------------------------
C     READ IN COULOMB INTEGRALS
C     -------------------------
      IKR = 0
      JKR = 0
      IF( MYPROC .eq. MASTER .or. IJKL_AV )THEN
        REWIND(LUTWO)
        IF(N2EFIL.eq.0) READ(LUTWO)
      END IF
C
  10  IKR_SAVE = IKR
      JKR_SAVE = JKR
      IF( (MYPROC .eq. MASTER) .or. IJKL_AV )THEN
        read (LUTWO) IKR,JKR,NZ,
     &               (indk(inz+ISYMOFF),indl(inz+ISYMOFF),inz=1,NZ,1),
     &               (rklr(inz+ISYMOFF),inz=1,NZ,1)
      END IF
#if defined (VAR_MPI)
C
C     broadcast if file MDCINTxxx is not available on all local disks
      IF( NMPROC .gt. 1 .and. (.not. IJKL_AV ))THEN
        CALL interface_MPI_BCAST(IKR,1,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(JKR,1,MASTER,
     &                 global_communicator)
        CALL interface_MPI_BCAST(NZ,1,MASTER,
     &                 global_communicator)
        DO J = 1, NZ, 1
           CALL interface_MPI_BCAST(INDK(J+ISYMOFF),1,MASTER,
     &                    global_communicator)
           CALL interface_MPI_BCAST(INDL(J+ISYMOFF),1,MASTER,
     &                    global_communicator)
           CALL interface_MPI_BCAST(RKLR(J+ISYMOFF),1,MASTER,
     &                    global_communicator)
        END DO
      END IF
#endif
C
      IF( IKR .eq. 0 )THEN
        IF( (MYPROC .eq. MASTER) .or. (IJKL_AV) )THEN
          close (LUTWO)
          write(LUWRT,*) 'End of file MDCINT.'
        END IF
        IF( JKR .eq. 0 )THEN
          write(LUWRT,*) 'Only or last 2-el. integral file.'
          goto 20
        ELSE IF(JKR.ne.0) then
          N2EFIL = N2EFIL + 1
          if (N2EFIL.le.9) then
          write(LUWRT,'(A,A7,I1)') 'Reading next file ','MDCINT0',N2EFIL
          else if (N2EFIL.gt.9.and.N2EFIL.le.99) then
          write(LUWRT,'(A,A6,I2)') 'Reading next file ','MDCINT',N2EFIL
          end if
          goto 30
        END IF
      END IF
C
      if (IKR.ne.IKR_SAVE.or.JKR.ne.JKR_SAVE) then
        IOFF12(KR(IKR),KR(JKR)) = ISYMOFF
        NUMINT(KR(IKR),KR(JKR)) = NZ
      else
C       redundant Kramers block has been read (symmetric)
        write(LUWRT,*) 'Redundant integral block detected.'
        write(LUWRT,*) 'IKR,IKR_SAVE,JKR,JKR_SAVE ',
     &                  IKR,IKR_SAVE,JKR,JKR_SAVE
        CALL QUIT( ' *** ERROR in RDTWO_DIRAC: redundant integral 
     & block detected. *** ')
      end if
*
#if defined LUCI_DEBUG
      if (NTEST.ge.10) then
        write(LUWRT,'(A,2I6)') 'IKR,JKR : ',IKR,JKR
        write(LUWRT,'(A,2I4,A,I6)')
     &        'IOFF12(',KR(IKR),KR(JKR),') = ',IOFF12(KR(IKR),KR(JKR))
        write(LUWRT,'(A,2I4,A,I6)')
     &        'NUMINT(',KR(IKR),KR(JKR),') = ',NUMINT(KR(IKR),KR(JKR))
        write(LUWRT,'(A,1X,2I6)') 'Integral block: ',KR(IKR),KR(JKR)
        if (NTEST.ge.100) then
          do INT=1,NZ,1
            write(LUWRT,'(A,2I8,1X,1D15.6)')
     &      'Indices K,L ',KR(INDK(INT+ISYMOFF)),KR(INDL(INT+ISYMOFF)),
     &                     RKLR(INT+ISYMOFF)
          end do
        end if
      end if
#endif
      ISYMOFF = ISYMOFF + NZ
      GOTO 10
C
 20   CONTINUE
#if defined LUCI_DEBUG
      if (NTEST.ge.50) then
        write(LUWRT,'(/a)') '===================='
        write(LUWRT,'(a )') 'rdtwo_dirac speaking'
        write(LUWRT,'(a )') '===================='
        write(LUWRT,'(a/)') 'Full two-electron integral list:'
        do IJ=1,MORB,1
          do II=1,MORB,1
            write(LUWRT,*) 'Block :',II,IJ
            write(LUWRT,*) 'number of elements :',NUMINT(II,IJ)
            do INT=1,NUMINT(II,IJ),1
              write(LUWRT,'(2I4,1D15.6)')
     &              KR(INDK(INT+IOFF12(II,IJ))),
     &              KR(INDL(INT+IOFF12(II,IJ))),
     &              RKLR(INT+IOFF12(II,IJ))
            end do
          end do
        end do
      end if
#endif
*
#if defined (VAR_MPI)
      if( NMPROC .gt. 1) CALL interface_MPI_BARRIER(global_communicator)
#endif
*
      RETURN
 1000 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Coulomb integrals written to MRCTWOE')
 1001 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Unique Breit integrals read from MDBINT (',A10,1X,A8,')'
     &/I8,' Combined integrals written to MRCTWOE')
 1002 FORMAT (/' Real arithmetic will be used in EXPAND and TWOLIN')
 1003 FORMAT (/' Complex arithmetic will be used')
 1010 FORMAT (//' CAUTION : No non-zero integrals found on ',A//)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE READMO(WRK,LWRK,XIJKL)
C
C     Written by Henrik Koch 27-Mar-1990
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LUINT = 27)
      PARAMETER (LUPRI = 6)
*. Modified for LUCIA : XIJKL added to list
C
C     -------------------------------------------------------
C     IRAT  = (real word length) / (integer word length)
C     IRAT2 = (real word length) / (half-integer word length)
C             if available and used, otherwise IRAT2 = IRAT
C     PARAMETER (IRAT = 2, IRAT2 = 2)
#include "irat.inc"
      DIMENSION WRK(LWRK)
      PARAMETER ( MAXRHF = 30, MAXVIR = 225)
      INTEGER P,Q,R,S,A,B,C,D,E,F,G
      COMMON /CIPOL / NBAST, NNBAST, MORB(8), NMORBT, NORB(8),NNORB(8),
     *                NSYMHF, LBINTM, LPPOP(3,3), LSYMOP(3,3),
     *                NCMOT, NNORBX, NLAMDA(8), LUEGVC,NORBT,ISAT(128)
C
      OPEN (LUINT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     *      FILE='MOTWOINT')
      REWIND LUINT
C
      REWIND LUINT
      READ (LUINT)
      READ (LUINT) LBINTM, JTRLVL
C?    WRITE(6,*) 'LBINTM   :  ',LBINTM,JTRLVL
C
      KONEMO = 1
      KTWOMO = KONEMO + NORBT*NORBT
C     KBUF   = KTWOMO + NORBT*NORBT*NORBT*NORBT
      KBUF   = KTWOMO
      KIBUF  = KBUF   + LBINTM
      KEND   = KIBUF  + LBINTM/IRAT + 1
C
      IF ( KEND .GT. LWRK ) THEN
         Call Abend2( 'Insufficient spaces in READMO' )
      ENDIF
C
C-----------------------------------
C     Read MO integrals into memory.
C-----------------------------------
C
      CALL REDMO1(XIJKL,WRK(KBUF),WRK(KIBUF),LBINTM,NORBT)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE REDMO1(TWOMO,BUF,IBUF,LBUF,NORBT)
C
C     Written by Henrik Koch 27-Mar-1990.
C
*
*. Modified to be SIRIUS-LUCIA interface, february 1993
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LUINT = 27)
      PARAMETER (LUPRI = 6)
C
C     -------------------------------------------------------
C     IRAT  = (real word length) / (integer word length)
C     IRAT2 = (real word length) / (half-integer word length)
C             if available and used, otherwise IRAT2 = IRAT
C     PARAMETER (IRAT = 2, IRAT2 = 2)
#include "irat.inc"
      DIMENSION TWOMO(*),
     *          BUF(LBUF),IBUF(LBUF)
      PARAMETER ( MAXRHF = 30, MAXVIR = 225)
      INTEGER P,Q,R,S,A,B,C,D,E,F,G
C----- bit manipulation definitions
      PARAMETER (IBT02=3, IBT08=255, IBT10=1023, IBT16=65535)
      PARAMETER (MYSHF=16,IBTMY=IBT16,MAXLN=16)
      IBTAND(I,J) = IAND(I,J)
C     IBTOR(I,J)  = IOR(I,J)
C     IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
C     IBTXOR(I,J) = IEOR(I,J)
C-----
C
      REWIND LUINT
      CALL MOLLAB('MOLTWOEL',LUINT,LUPRI)
C
      INDCD = 0
      NINTR = 0
  200 READ (LUINT) BUF,IBUF,LENGTH
      IF (LENGTH .EQ. 0) GOTO 200
      IF (LENGTH .EQ. -1) GOTO 9500
      INDCDN = IBTAND(IBTSHR(IBUF(1),16),IBT16)
      IF ( INDCDN .NE. INDCD ) THEN
         INDCD = INDCDN
         IC    = IBTAND(IBTSHR(INDCD,8),IBT08)
         ID    = IBTAND(       INDCD,   IBT08)
      ENDIF
      DO 280 I = 1,LENGTH
         IA = IBTAND(IBTSHR(IBUF(I),8),IBT08)
         IB = IBTAND(       IBUF(I),   IBT08)
         NINTR = NINTR + 1
	 IABCD = I2EAD(IA,IB,IC,ID)
C         write(6,*) ' IA IB IC ID IABCD ',IA,IB,IC,ID,IABCD
	 TWOMO(IABCD) = BUF(I)
c
  280 CONTINUE
      GOTO 200
 9500 CONTINUE
*
      WRITE(6,*) ' Number of integrals read ', NINTR
      WRITE(6,*) ' Indeces of last integral ', IA,IB,IC,ID
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLH2(XLAMBDA)
*
*. Scale two electron integrals
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      CALL SCALVE(WORK(KINT2),XLAMBDA,NINT2)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
* Compute boson symmetry information from DIRAC
* common blocks.
*
*    Timo Fleig, Jan. 2002
*
      subroutine transf_dirac(IDBGRP,NBSYM,
     &                        NISH_L,NASH_L,NOCC_L,NORB_L,NSSH_L,
     &                        NFSYM,NZ,MORB,IBOSYM,NSTR,NCORE,
     &                        NDELE,NACTELD,WAFFCT,IRETISH,IPRNT)
*
      use interface_to_mpi
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "dcbbas.h"
#include "dcborb.h"
#include "infpar.h"
#include "parluci.h"
*
      character*6 WAFFCT
      dimension IBOSYM(MORB),ICT(MXPIRR),NSTR(NFSYM),NSSH_L(NBSYM),
     &          NORB_L(NBSYM),NISH_L(NBSYM),NASH_L(NBSYM),NOCC_L(NBSYM),
     &          NCORE(NFSYM),NDELE(NFSYM),
     &          IFSOFF(2),NVIR(2)
*
C
C     arrange for the MPI stuff and correct node number
C     to the total number of running invocations.
C
#if defined (VAR_MPI)
      MASTER = MPARID
      MYPROC = MYTID
C     Add the master node, NUMNOD = number of slaves
      NMPROC = NUMNOD + 1
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
C
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
* Number of functions in each boson irrep
      do I=1,MXPIRR,1
        ICT(I) = 0
      end do
      do ISPINOR=1,MORB,1
        ISYM = IBOSYM(ISPINOR) + 1
        ICT(ISYM) = ICT(ISYM) + 1
      end do
*
* Determine number of orbitals per irrep in subspaces
      do ISM = 1,NBSYM,1
        NISH_L(ISM) = 0
        NASH_L(ISM) = 0
        NSSH_L(ISM) = 0
      end do
*
      do IFS = 1,NFSYM,1
        if (NTEST.ge.2) then
          write(6,*) 'Info from DIRAC common block dcborb.h :'
          write(6,*) '(In terms of Kramers pairs)'
          write(6,*) 'Fermion symmetry        : ',IFS
          write(6,*) 'Number of Kramers pairs : ',NSTR(IFS)
          write(6,*) 'Moltra core Kra pai     : ',NCORE(IFS)
          write(6,*) 'Dirac Inactive orbitals : ',NISH(IFS)
          write(6,*) 'Dirac Active   orbitals : ',NASH(IFS)
          write(6,*) 'Dirac External orbitals : ',NSSH(IFS)
          write(6,*) 'Moltra deleted Kra pai  : ',NDELE(IFS)
        end if
*
* Offset for start of fermion irrep on IBOSYM array
        ISUM = 0
        do IIFS = 1,IFS-1,1
          ISUM = ISUM + 2*NSTR(IIFS)
        end do
        IFSOFF(IFS) = ISUM
*
* Virtual set in MOLTRA
        NVIR(IFS) = NSSH(IFS) - NDELE(IFS)
*
*
        ISTAI = IFSOFF(IFS) + 1
        IENDI = ISTAI + NISH(IFS) - NCORE(IFS) - 1
        if (NTEST.ge.10) write(6,*) 'ISTAI,IENDI ',ISTAI,IENDI
        do IISH = ISTAI,IENDI,1
          ISYM = IBOSYM(IISH) + 1
          NISH_L(ISYM) = NISH_L(ISYM) + 1
        end do
*
        ISTAA = IENDI + 1
        IENDA = ISTAA + NASH(IFS) - 1
        if (NTEST.ge.10) write(6,*) 'ISTAA,IENDA ',ISTAA,IENDA
        do IIASH = ISTAA,IENDA,1
          ISYM = IBOSYM(IIASH) + 1
          NASH_L(ISYM) = NASH_L(ISYM) + 1
        end do
*
        ISTAE = IENDA + 1
        IENDE = ISTAE + NVIR(IFS) - 1
        if (NTEST.ge.10) write(6,*) 'ISTAE,IENDE ',ISTAE,IENDE
        do ISSH = ISTAE,IENDE,1
          ISYM = IBOSYM(ISSH) + 1
          NSSH_L(ISYM) = NSSH_L(ISYM) + 1
        end do
      end do
*     ^ End loop over fermion symmetries.
*
      if (NTEST.ge.5) then
        write(6,*) 'transf_dirac information:'
        write(6,'(A,8I3)') 'NISH_L ',(NISH_L(I),I=1,NBSYM,1)
        write(6,'(A,8I3)') 'NASH_L ',(NASH_L(I),I=1,NBSYM,1)
        write(6,'(A,8I3)') 'NSSH_L ',(NSSH_L(I),I=1,NBSYM,1)
      end if
*
C       this is necessary! - sk -
C 
#if defined (VAR_MPI)
      IF (NMPROC.GT.1) THEN
      CALL interface_MPI_BCAST(NBSYM,1,
     &                         MASTER,global_communicator)
      DO I = 1,NBSYM,1 
      CALL interface_MPI_BCAST(NISH_L(I),1,
     &                         MASTER,global_communicator)
      CALL interface_MPI_BCAST(NASH_L(I),1,
     &                         MASTER,global_communicator)
      CALL interface_MPI_BCAST(NSSH_L(I),1,
     &                         MASTER,global_communicator)
      END DO
      END IF
#endif
C
C
      IERROR = 0
      IRETISH = -1
      do ISM = 1,NBSYM,1
        NOCC_L(ISM) = NISH_L(ISM) + NASH_L(ISM)
        if (NOCC_L(ISM).gt.MXTSOB) then
          write(6,*) ' Too many occupied orbitals per symmetry.'
          write(6,*) ' NOCC_L(ISM=',ISM,') = ', NOCC_L(ISM)
          write(6,*) ' Use GASCI setup (user option) '
          write(6,*) '  or increase parameter MXTSOB from: ',MXTSOB
          stop 'Quitting.'
        end if
        NORB_L(ISM) = NOCC_L(ISM) + NSSH_L(ISM)
        if ((ICT(ISM)/2).ne.NORB_L(ISM)) then
          IERROR = 1
        end if
        if (NISH_L(ISM).gt.0) IRETISH = 0
        if (NISH_L(ISM).lt.0) then
          write(6,*) 'Impossible number of inactive shells.'
          write(6,*) 'NISH_L for symm. ',ISM,' is ',NISH_L(ISM)
          stop 'Quitting.'
        end if
      end do
C
      if (IERROR.ne.0) then
        write(6,*) 'Problem in transf_dirac.'
        do ISM = 1,NBSYM,1
          write(6,*) 'ICT(ISM)/2 for ISM =',ISM,' is ',ICT(ISM)/2
          write(6,*) ' NISH_L(ISM) is ',NISH_L(ISM)
          write(6,*) ' NOCC_L(ISM) is ',NOCC_L(ISM)
          write(6,*) ' NORB_L(ISM) is ',NORB_L(ISM)
        end do
        stop 'Quitting.'
      end if
*
* Determine double group label
      if (NFSYM.eq.1) then
        if (NZ.eq.4) IDBGRP = 8
        if (NZ.eq.2) IDBGRP = 5
        if (NZ.eq.1) IDBGRP = 2
      else if (NFSYM.eq.2) then
        if (NZ.eq.4) IDBGRP = 7
        if (NZ.eq.2) IDBGRP = 4
        if (NZ.eq.1) IDBGRP = 1
      end if
*
      if (NTEST.ge.3) then
        write(6,'(A,8I4)') 'Inactive orbs per symm : ',
     &                      (NISH_L(I),I=1,NBSYM,1)
        write(6,'(A,8I4)') 'Active   orbs per symm : ',
     &                      (NASH_L(I),I=1,NBSYM,1)
        write(6,'(A,8I4)') 'Occupied orbs per symm : ',
     &                      (NOCC_L(I),I=1,NBSYM,1)
        write(6,'(A,8I4)') 'External orbs per symm : ',
     &                      (NSSH_L(I),I=1,NBSYM,1)
        write(6,'(A,8I4)') 'Total    orbs per symm : ',
     &                      (NORB_L(I),I=1,NBSYM,1)
        write(6,'(A,1X,I3)') ' Double group ',IDBGRP
      end if
*
      if (WAFFCT.eq.'HF_SCF'.and.NAELEC.gt.0) then
        write(6,*) 'You have specified open-shell electrons in DIRAC.'
        write(6,*) 'Im am not allowed to run WAFFCT = DHFSCF'
        write(6,*) 'in this case. Use OSHSCF instead.'
        stop 'Quitting.'
      end if
      NACTELD = NAELEC
*
      return
      end
