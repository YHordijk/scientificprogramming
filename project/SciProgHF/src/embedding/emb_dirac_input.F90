
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

module fde_input_dirac

   use fde_types
   use fde_cfg
   use fde_mag_cfg
   use fde_data
   use fde_nadd_derv
   use fde_input
   
   private

   public fde_dirac_input
   
   contains

!   ----------------------------------------------------------------------------
      subroutine fde_dirac_input(word,echo_input)
!   ----------------------------------------------------------------------------
         use dirac_cfg
#include "implicit.h"
#include "priunit.h" 
#include "dcbgen.h"
      PARAMETER (NTABLE = 8)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, line*80
 
      logical :: Echo_Input,newdef
      integer :: pl
      integer :: lunit
      type(fde_files)  :: ftmp
      type(fde_export) :: etmp
      type(fde_import) :: itmp
      character(len=80):: nadd_fun
!
      call fde_get_files_info(ftmp)
      call fde_get_export_info(etmp)
      call fde_get_import_info(itmp)
      
      lunit = ftmp%logfile%unit
      
      ICHANG = 0
      IF (WORD(1:4).eq.'*FDE') THEN
        WORD1 = WORD
        DO
          READ (LUCMD, '(A7)') WORD
          CALL UPCASE(WORD)
          SELECT CASE (WORD(1:1))
          CASE('!','#')
            CYCLE
          CASE('.')
            ICHANG = ICHANG + 1
            SELECT CASE(WORD(1:7))
            CASE ('.OPTION')
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',lunit)
              CYCLE

            CASE ('.EMBPOT')
               read(lucmd,'(A60)') ftmp%embpot%name
               itmp%im_vemb        = .true.
               itmp%im_update_vemb = .false.
               
            CASE ('.PRINT')
              read(lucmd,*) pl
              call fde_set_print_level(pl)

            CASE ('.FRDENS')
               read(lucmd,'(A60)') ftmp%frozen%name
               itmp%im_frozen   = .true.

            CASE ('.UPDATE')
               itmp%im_update_vemb = .true.

            CASE ('.NUCREP')
               read(lucmd,*) fde_intersub_nucrep

            CASE ('.ENER_B')
               read(lucmd,*) fde_frozen_energy 

            CASE ('.OLDESP')
               fde_export_old_esp = .true. 

            CASE ('.SKIPK')
               call fde_set_skip_nadd_ke(.true.)
            CASE ('.SKIPX')
               call fde_set_skip_nadd_xc(.true.)

! read kinetic energy functional for non-additive term
            CASE ('.NAXCF')
               read(lucmd,'(A80)') nadd_fun
               call fde_set_nadd_xcf(nadd_fun)

! read kinetic energy functional for non-additive term
            CASE ('.NAKEF')
               read(lucmd,*) nadd_fun
               call fde_set_nadd_kef(nadd_fun)

! flags for response read kinetic energy functional for non-additive term
            CASE ('.RSP')
               itmp%im_frozen = .true.
               dirac_cfg_fde_response = .true.

            CASE ('.RSPLDA')
               itmp%im_frozen = .true.
               call fde_set_nadd_all_alda
               dirac_cfg_fde_response = .true.

            CASE ('.LEVEL')
               read(lucmd,*) etmp%ex_level 

            CASE ('.GRIDOU')
               read(lucmd,*) ftmp%export%name
               dirac_cfg_fde_export = .true.
               etmp%do_grid_out = .true.
! doprp is part of dcbgen.h
               doprp = .true.

            CASE ('.WRFRMT')
               read(lucmd,*) etmp%ex_format

            CASE ('.RDFRMT')
               read(lucmd,*) itmp%im_format

            CASE ('.EXONLY')
               read(lucmd,*) ftmp%export%name
               dirac_cfg_fde_export = .true.
               etmp%do_grid_out = .true.
! doprp is part of dcbgen.h
               skip_all_imports = .true.
               doprp = .true.

            CASE ('.NOSDFT')
!              no spin density contributions
               fde_cfg_no_sdft = .true.

! gosia: keywords for properties in LAO basis sets:
!           FDE-LAO contributions to property gradient:
!           -------------------------------------------
            CASE ('.LDPT')
!              Direct-LAO contribution to property gradient, involving embedding PoTential
               fde_rsp_propgrad_lao = .true.
               fde_rsp_propgrad_lao_direct_embpot = .true.
               fde_exclude_rsp_propgrad_lao_direct_embpot = .false.
            CASE ('.NOLDPT')
!              Exclude Direct-LAO contribution to property gradient, involving embedding PoTential
               fde_rsp_propgrad_lao_direct_embpot = .false.
               fde_exclude_rsp_propgrad_lao_direct_embpot = .true.
            CASE ('.LDKR')
!              Direct-LAO contribution to property gradient, involving embedding KeRnel (nonzero if SDFT)
               fde_rsp_propgrad_lao = .true.
               fde_rsp_propgrad_lao_direct_embker = .true.
               fde_exclude_rsp_propgrad_lao_direct_embker = .false.
            CASE ('.NOLDKR')
!              Exclude Direct-LAO contribution to property gradient, involving embedding KeRnel
               fde_rsp_propgrad_lao_direct_embker = .false.
               fde_exclude_rsp_propgrad_lao_direct_embker = .true.
            CASE ('.LRKR')
!              Reorthonormalization-LAO contribution to property gradient, involving embedding KeRnel
               fde_rsp_propgrad_lao = .true.
               fde_rsp_propgrad_lao_reorth_embker = .true.
               fde_exclude_rsp_propgrad_lao_reorth_embker = .false.
            CASE ('.NOLRKR')
!              Exclude Reorthonormalization-LAO contribution to property gradient, involving embedding KeRnel
               fde_rsp_propgrad_lao_reorth_embker = .false.
               fde_exclude_rsp_propgrad_lao_reorth_embker = .true.
            CASE ('.L11KR')
!              Group all LAO contributions to property gradient involving uncoupled embedding kernel (w11)
               fde_rsp_propgrad_lao = .true.
               fde_rsp_propgrad_lao_w11 = .true.
               fde_rsp_propgrad_lao_direct_embker = .true.
               fde_rsp_propgrad_lao_reorth_embker = .true.
               fde_exclude_rsp_propgrad_lao_w11 = .false.
               fde_exclude_rsp_propgrad_lao_direct_embker = .false.
               fde_exclude_rsp_propgrad_lao_reorth_embker = .false.
            CASE ('.NL11KR')
!              Exclude all LAO contributions to property gradient involving uncoupled embedding kernel (w11)
               fde_rsp_propgrad_lao_w11 = .false.
               fde_rsp_propgrad_lao_direct_embker = .false.
               fde_rsp_propgrad_lao_reorth_embker = .false.
               fde_exclude_rsp_propgrad_lao_w11 = .true.
               fde_exclude_rsp_propgrad_lao_direct_embker = .true.
               fde_exclude_rsp_propgrad_lao_reorth_embker = .true.
            CASE ('.LAO11')
!              Group all LAO contributions to property gradient involving embedding PoTential (v) and uncoupled embedding kernel (w11)
               fde_rsp_propgrad_lao = .true.
               fde_rsp_propgrad_lao_vw11 = .true.
               fde_rsp_propgrad_lao_direct_embpot = .true.
               fde_rsp_propgrad_lao_direct_embker = .true.
               fde_rsp_propgrad_lao_reorth_embker = .true.
               fde_exclude_rsp_propgrad_lao_direct_embpot = .false. 
               fde_exclude_rsp_propgrad_lao_direct_embker = .false.
               fde_exclude_rsp_propgrad_lao_reorth_embker = .false. 
!           keywords related to coupled-kernel contributions (to do):
            CASE ('.LAOFRZ')
!              LAO contribution to property gradient from the FRoZen density (read from file)
!              it depends on embedding kernel and does NOT include the Coulomb term...
               fde_rsp_propgrad_lao = .true.
               fde_lao_frozen_embker_nonadd = .true.
            CASE ('.LFCOUL')
!              ...this is the additional Coulomb term that contributes to the coupling-kernel-like term (.LAOFRZ)
               fde_rsp_propgrad_lao = .true.
               fde_lao_frozen_embker_coulomb = .true.
            CASE ('.LAO12')
!              Group all LAO contributions to property gradient involving coupling kerned terms (w12)
               fde_rsp_propgrad_lao = .true.
               fde_lao_frozen_embker_nonadd = .true.
               fde_lao_frozen_embker_coulomb = .true.
!           import/export of perturbed density:
!           -----------------------------------
            CASE ('.PERTIM')
!              import perturbed density (additional columns on 'GRIDOUT' file)
               read(lucmd,'(A60)') ftmp%frozen_pert_direct%name
               read(lucmd,'(A60)') ftmp%frozen_pert_reorth%name
               itmp%im_frozen_pert_direct = .true.
               itmp%im_frozen_pert_reorth = .true.
               itmp%im_frozen = .true.
               fde_rsp_mag_lao_import = .true.
            CASE ('.PERTEX')
!              export perturbed density (additional columns on 'GRIDOUT' file)
               dirac_cfg_fde_export = .true.
               etmp%do_grid_out = .true.
               fde_rsp_mag_lao_export = .true.
!           embedding contributions to expectation value of magnetizability:
!           ----------------------------------------------------------------
            CASE ('.EMAFDE')
!              Expectation value contribution to MAgnetizability - FDE (all terms, excluding terms by keywords below)
!              it has to be called with PERTIM which imports perturbed density
               fde_magn_expval_lao = .true.
            CASE ('.MNOPOT')
!              Magnetizability NO embedding POTential
               fde_lao_magn_expval_no_embpot = .true.
            CASE ('.MNOUKE')
!              Magnetizability NO Uncoupled KErnel
               fde_lao_magn_expval_no_uncoup_embker = .true.
            CASE ('.MNONKE')
!              Magnetizability NO coupled Nonadditive KErnel
               fde_lao_magn_expval_no_coupl_nonadd_embker = .true.
 

            CASE DEFAULT
              WRITE (lunit,'(/,3A,/)') ' Keyword "',WORD, &
                 '" not recognized in FDEINP.'
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',lunit)
              CALL QUIT('Illegal keyword in FDEINP.')
            END SELECT

          CASE('*')
            EXIT

          CASE DEFAULT
            WRITE (lunit,'(/,3A,/)') ' Prompt "',WORD, &
               '" not recognized in FDEINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',lunit)
            CALL QUIT('Illegal prompt in FDEINP.')

          END SELECT
        END DO
      END IF

! saving information read 
      call fde_set_files_info(ftmp)
      call fde_set_export_info(etmp)
      call fde_set_import_info(itmp)
      
! now these will be validated, and changed if needed
      call fde_input_validate_and_postprocess

!     options are echoed
         if (Echo_Input) then

           write (lunit,*) ' '
           CALL PRSYMB(lunit,'=',75,0)
           write (lunit,*) ' '

           write (lunit,'(1X,A)') &
              ' Frozen Density Embedding (FDE) Calculation setup'  

           call fde_input_show_setup 

           CALL PRSYMB(lunit,'=',75,0)
           write (lunit,*) ' '
         endif 

      END SUBROUTINE

end module fde_input_dirac
