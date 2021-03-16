module fde_mag_cfg


   implicit none

   public fde_rsp_propgrad_lao
   public fde_rsp_propgrad_lao_direct_embpot
   public fde_rsp_propgrad_lao_direct_embker
   public fde_rsp_propgrad_lao_reorth_embker
   public fde_rsp_propgrad_lao_w11
   public fde_rsp_propgrad_lao_vw11

   public fde_exclude_rsp_propgrad_lao_direct_embpot
   public fde_exclude_rsp_propgrad_lao_direct_embker
   public fde_exclude_rsp_propgrad_lao_reorth_embker
   public fde_exclude_rsp_propgrad_lao_w11


   public fde_rsp_mag_lao_import
   public fde_rsp_mag_lao_export
   public fde_lao_frozen_embker_nonadd
   public fde_lao_frozen_embker_coulomb

   public fde_magn_expval_lao
   public fde_lao_magn_expval_no_embpot
   public fde_lao_magn_expval_no_uncoup_embker
   public fde_lao_magn_expval_no_coupl_nonadd_embker



!  import/export of perturbed density:
   logical, save :: fde_rsp_mag_lao_import        = .false.
   logical, save :: fde_rsp_mag_lao_export        = .false.

!  flags related to FDE-LAO contributions to property gradient:
   logical, save :: fde_rsp_propgrad_lao                           = .false.
   logical, save :: fde_rsp_propgrad_lao_direct_embpot             = .false.
   logical, save :: fde_rsp_propgrad_lao_direct_embker             = .false.
   logical, save :: fde_rsp_propgrad_lao_reorth_embker             = .false.
   logical, save :: fde_rsp_propgrad_lao_w11                       = .false.
   logical, save :: fde_rsp_propgrad_lao_vw11                      = .false.
   logical, save :: fde_exclude_rsp_propgrad_lao_direct_embpot     = .false.
   logical, save :: fde_exclude_rsp_propgrad_lao_direct_embker     = .false.
   logical, save :: fde_exclude_rsp_propgrad_lao_reorth_embker     = .false.
   logical, save :: fde_exclude_rsp_propgrad_lao_w11               = .false.

   logical, save :: fde_lao_frozen_embker_nonadd  = .false.
   logical, save :: fde_lao_frozen_embker_coulomb = .false.

!  flags related to FDE-LAO contributions to expectation value part of the magnetizability tensor:
   logical, save :: fde_magn_expval_lao = .false.
   logical, save :: fde_lao_magn_expval_no_embpot = .false.
   logical, save :: fde_lao_magn_expval_no_uncoup_embker = .false.
   logical, save :: fde_lao_magn_expval_no_coupl_nonadd_embker = .false.


end module fde_mag_cfg
