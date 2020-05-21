module dynSubgridAdjustmentsMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Holds the methods for adjusting state variables at each sub-grid level,
  ! based on weight changes and other information stored in the
  ! state_updater structures.
  ! These subroutines were originally contained in the sub-grid level
  ! data type definitions, but moved here to keep all the dynamic subgrid
  ! functionality together.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use decompMod              , only : bounds_type
  use clm_varpar             , only : ndecomp_pools, nlevdecomp
  use clm_varctl             , only : use_crop
  use clm_varcon             , only : dzsoi_decomp
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use dynColumnStateUpdaterMod, only : column_state_updater_type
  use ColumnDataType         , only : column_carbon_state, column_nitrogen_state
  use ColumnDataType         , only : column_phosphorus_state
  use VegetationDataType     , only : vegetation_carbon_state, vegetation_nitrogen_state
  use VegetationDataType     , only : vegetation_phosphorus_state
  use SpeciesMod             , only : CN_SPECIES_N, CN_SPECIES_P
  !---Currently contains module class methods for gpu
  use dynUpdateModAcc
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  real(r8) , parameter :: npool_seed_param     = 0.1_r8
  real(r8) , parameter :: ppool_seed_param     = 0.01_r8
  !$acc declare copyin(npool_seed_param, ppool_seed_param)

  public :: dyn_veg_cs_Adjustments
  public :: dyn_col_cs_Adjustments
  public :: dyn_veg_ns_Adjustments
  public :: dyn_col_ns_Adjustments
  public :: dyn_veg_ps_Adjustments
  public :: dyn_col_ps_Adjustments
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dyn_veg_cs_Adjustments(        &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafc_seed,                      &
       dwt_deadstemc_seed,                  &
       conv_cflux,                          &
       dwt_frootc_to_litter,                &
       dwt_livecrootc_to_litter,            &
       dwt_deadcrootc_to_litter,            &
       prod10_cflux,                        &
       prod100_cflux,                       &
       crop_product_cflux,                  &
       veg_cs                               &
       )
    !
    ! !DESCRIPTION:
    ! Adjust veg-level state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:) ! patch filter that includes inactive points
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafc_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemc_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_cflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootc_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_cflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_cflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_cflux       (bounds%begp:)
    type(vegetation_carbon_state)  , intent(inout) :: veg_cs
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical                     :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafc_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_storage_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafc_xfer_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemc_patch(bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_cflux(bounds%begp:bounds%endp)
    real(r8)                    :: deadstemc_patch_temp(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'CStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp


    old_weight_was_zero = old_weight_was_zeroAcc(patch_state_updater, bounds)
    patch_grew          = patch_grewAcc(patch_state_updater,bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = veg_cs%species                          , &
         leaf_patch                 = veg_cs%leafc           , &
         leaf_storage_patch         = veg_cs%leafc_storage   , &
         leaf_xfer_patch            = veg_cs%leafc_xfer      , &
         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew               , &
         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero      , &
         seed_leaf_patch            = seed_leafc_patch         , &
         seed_leaf_storage_patch    = seed_leafc_storage_patch , &
         seed_leaf_xfer_patch       = seed_leafc_xfer_patch    , &
         seed_deadstem_patch        = seed_deadstemc_patch   )

    ! 1) LEAFC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%leafc              , &
         flux_out_grc_area = conv_cflux                  , &
         seed              = seed_leafc_patch            , &
         seed_addition     = dwt_leafc_seed   )


    ! 2) LEAFC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%leafc_storage      , &
         flux_out_grc_area = conv_cflux                  , &
         seed              = seed_leafc_storage_patch    , &
         seed_addition     = dwt_leafc_seed           )

    ! 3) LEAF_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%leafc_xfer         , &
         flux_out_grc_area = conv_cflux                  , &
         seed              = seed_leafc_xfer_patch       , &
         seed_addition     = dwt_leafc_seed        )

    ! 4) FROOTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%frootc             , &
         flux_out_col_area = dwt_frootc_to_litter)

    ! 5) FROOTC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%frootc_storage     , &
         flux_out_grc_area = conv_cflux)

    ! 6) FROOTC_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%frootc_xfer        , &
         flux_out_grc_area = conv_cflux)

    ! 7) LIVESTEMC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livestemc          , &
         flux_out_grc_area = conv_cflux)

    ! 8) LIVESTEMC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livestemc_storage  , &
         flux_out_grc_area = conv_cflux)

    ! 9) LIVESTEMC_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livestemc_xfer     , &
         flux_out_grc_area = conv_cflux)

    ! 10) PROD10_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = veg_cs%deadstemc(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(&
         patch_state_updater    ,&
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod10                          , &
         var                        = deadstemc_patch_temp     , &
         flux1_out                  = prod10_cflux             , &
         flux2_out                  = wood_product_cflux       , &
         seed                       = seed_deadstemc_patch     )

    ! 11) PROD100_CFLUX
    wood_product_cflux(begp:endp)      = 0._r8
    deadstemc_patch_temp(begp:endp)    = veg_cs%deadstemc(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(&
         patch_state_updater    ,&
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pprod100                         , &
         var                        = deadstemc_patch_temp     , &
         flux1_out                  = prod100_cflux             , &
         flux2_out                  = wood_product_cflux       , &
         seed                       = seed_deadstemc_patch    )

    ! 12) DEADSTEMC_PATCH
    wood_product_cflux(begp:endp)      = 0._r8
    call update_patch_state_partition_flux_by_typeAcc(&
         patch_state_updater    ,&
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         flux1_fraction_by_pft_type = pconv                          , &
         var                        = veg_cs%deadstemc , &
         flux1_out                  = conv_cflux               , &
         flux2_out                  = wood_product_cflux        , &
         seed                       = seed_deadstemc_patch     , &
         seed_addition              = dwt_deadstemc_seed     )

    ! 13) DEADSTEMC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadstemc_storage  , &
         flux_out_grc_area = conv_cflux)

    ! 14) DEADSTEMC_XFER_PATCH
     call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadstemc_xfer     , &
         flux_out_grc_area = conv_cflux)

    ! 15) LIVECROOTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livecrootc         , &
         flux_out_col_area = dwt_livecrootc_to_litter)

    ! 16) LIVECROOTC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livecrootc_storage , &
         flux_out_grc_area = conv_cflux)

    ! 17) LIVECROOTC_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%livecrootc_xfer    , &
         flux_out_grc_area = conv_cflux)

    ! 18) DEADCROOTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadcrootc         , &
         flux_out_col_area = dwt_deadcrootc_to_litter)

    ! 19) DEADCROOTC_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadcrootc_storage , &
         flux_out_grc_area = conv_cflux)

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%deadcrootc_xfer    , &
         flux_out_grc_area = conv_cflux)

    ! 21) GRESP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%gresp_storage      , &
         flux_out_grc_area = conv_cflux)

    ! 22) GRESP_XFER_STORAGE
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%gresp_xfer         , &
         flux_out_grc_area = conv_cflux)

    ! 23) CPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%cpool              , &
         flux_out_grc_area = conv_cflux)

    ! 24) XSMRPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%xsmrpool           , &
         flux_out_grc_area = conv_cflux)

    ! 25) CTRUNC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%ctrunc             , &
         flux_out_grc_area = conv_cflux)

    ! 26) DISPVEGC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%dispvegc)

    ! 27) STORVEGC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%storvegc)

    ! 28) TOTVEGC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%totvegc)

    ! 29) TOTPFTC_PATCH
    call update_patch_stateAcc(patch_state_updater,                 &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_cs%totpftc)

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_stateAcc(patch_state_updater,                 &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = veg_cs%cropseedc_deficit , &
            flux_out_grc_area = conv_cflux)
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp, endp
       dwt_frootc_to_litter(p)     = -1._r8 * dwt_frootc_to_litter(p)
       dwt_livecrootc_to_litter(p) = -1._r8 * dwt_livecrootc_to_litter(p)
       dwt_deadcrootc_to_litter(p) = -1._r8 * dwt_deadcrootc_to_litter(p)
    end do

  end subroutine dyn_veg_cs_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_cs_Adjustments(bounds, clump_index, column_state_updater, col_cs)
    !
    ! !DESCRIPTION:
    ! Adjust column-level carbon state variables and compute associated fluxes
    ! when patch areas change due to dynamic landuse
    !
    ! !USES:
    !
    ! !ARGUMENTS:
      !$acc routine seq
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_carbon_state)       , intent(inout) :: col_cs
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    real(r8)     ::   col_slice(bounds%begc:bounds%endc)
    !character(len=*), parameter :: subname = 'CStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    col_cs%dyn_cbal_adjustments(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
         col_slice = col_cs%decomp_cpools_vr(begc:endc, j, l)
          call update_column_state_no_special_handling_acc( column_state_updater , &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = col_slice,     &
               adjustment  = adjustment_one_level)

          col_cs%dyn_cbal_adjustments(begc:endc) = &
               col_cs%dyn_cbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

          col_cs%decomp_cpools_vr(begc:endc, j, l) = col_slice
       end do
    end do

    do j = 1, nlevdecomp
      col_slice = col_cs%ctrunc_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater , &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,     &
            adjustment  = adjustment_one_level)

       col_cs%dyn_cbal_adjustments(begc:endc) = &
            col_cs%dyn_cbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

        col_cs%ctrunc_vr(begc:endc,j) = col_slice
    end do

  end subroutine dyn_col_cs_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_veg_ns_Adjustments(        &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafn_seed,                      &
       dwt_deadstemn_seed,                  &
       dwt_npool_seed,                      &
       conv_nflux,                          &
       dwt_frootn_to_litter,                &
       dwt_livecrootn_to_litter,            &
       dwt_deadcrootn_to_litter,            &
       prod10_nflux,                        &
       prod100_nflux,                       &
       crop_product_nflux,                  &
       veg_ns                               &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafn_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemn_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_npool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_nflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootn_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_nflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_nflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_nflux       (bounds%begp:)
    type(vegetation_nitrogen_state), intent(inout) :: veg_ns
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafn_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemn_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_npool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_nflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemn_patch_temp     (bounds%begp:bounds%endp)

    !character(len=*), parameter :: subname = 'NStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp


    old_weight_was_zero = old_weight_was_zeroAcc(patch_state_updater,bounds)
    patch_grew          = patch_grewAcc(patch_state_updater,bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_N                        , &
         leaf_patch                 = veg_ns%leafn           , &
         leaf_storage_patch         = veg_ns%leafn_storage   , &
         leaf_xfer_patch            = veg_ns%leafn_xfer      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafn_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafn_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafn_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemn_patch(begp:endp)     , &
         pool_seed_param            = npool_seed_param                    , &
         pool_seed_patch            = seed_npool_patch(begp:endp))

    ! 1) LEAFN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = veg_ns%leafn    , &
         flux_out_grc_area = conv_nflux        , &
         seed              = seed_leafn_patch  , &
         seed_addition     = dwt_leafn_seed   )

    ! 2) LEAFN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = veg_ns%leafn_storage    , &
         flux_out_grc_area = conv_nflux                , &
         seed              = seed_leafn_storage_patch  , &
         seed_addition     = dwt_leafn_seed           )

    ! 3) LEAFN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = veg_ns%leafn_xfer   , &
         flux_out_grc_area = conv_nflux            , &
         seed              = seed_leafn_xfer_patch , &
         seed_addition     = dwt_leafn_seed        )

    ! 4) FROOTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%frootn             , &
         flux_out_col_area = dwt_frootn_to_litter)

    ! 5) FROOTN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%frootn_storage     , &
         flux_out_grc_area = conv_nflux)

    ! 6) FROOTN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%frootn_xfer        , &
         flux_out_grc_area = conv_nflux)

    ! 7) LIVESTEMN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%livestemn          , &
         flux_out_grc_area = conv_nflux)

    ! 8) LIVESTEMN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%livestemn_storage  , &
         flux_out_grc_area = conv_nflux)

    ! 9) LIVESTEMN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%livestemn_xfer     , &
         flux_out_grc_area = conv_nflux)

    ! 10) PROD10_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = veg_ns%deadstemn(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater , &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemn_patch_temp     , &
         flux1_out                  = prod10_nflux             , &
         flux2_out                  = wood_product_nflux       , &
         seed                       = seed_deadstemn_patch     )

    ! 11) PROD100_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = veg_ns%deadstemn(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater , &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemn_patch_temp     , &
         flux1_out                  = prod100_nflux            , &
         flux2_out                  = wood_product_nflux       , &
         seed                       = seed_deadstemn_patch    )

    ! 12) DEADSTEMN_PATCH
    wood_product_nflux(begp:endp)      = 0._r8
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater , &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = veg_ns%deadstemn       , &
         flux1_out                  = conv_nflux               , &
         flux2_out                  = wood_product_nflux       , &
         seed                       = seed_deadstemn_patch     , &
         seed_addition              = dwt_deadstemn_seed   )

    ! 13) DEADSTEMN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%deadstemn_storage  , &
         flux_out_grc_area = conv_nflux)

    ! 14) DEADSTEMN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%deadstemn_xfer     , &
         flux_out_grc_area = conv_nflux)

    ! 15) LIVECROOTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%livecrootn         , &
         flux_out_col_area = dwt_livecrootn_to_litter)

    ! 16) LIVECROOTN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%livecrootn_storage , &
         flux_out_grc_area = conv_nflux)

    ! 17) LIVECROOTN_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%livecrootn_xfer    , &
         flux_out_grc_area = conv_nflux)

    ! 18) DEADCROOTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%deadcrootn         , &
         flux_out_col_area = dwt_deadcrootn_to_litter)

    ! 19) DEADCROOTN_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%deadcrootn_storage , &
         flux_out_grc_area = conv_nflux)

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%deadcrootn_xfer    , &
         flux_out_grc_area = conv_nflux)

    ! 21) RETRANSN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%retransn           , &
         flux_out_grc_area = conv_nflux)

    ! 22) NTRUNC_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%ntrunc             , &
         flux_out_grc_area = conv_nflux)

    ! 23) CPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%npool              , &
         flux_out_grc_area = conv_nflux)

    ! 24) DISPVEGN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%dispvegn)

    ! 25) STORVEGN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%storvegn)

    ! 26) TOTVEGN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%totvegn)

    ! 27) TOTPFTN_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ns%totpftn)

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_stateAcc(patch_state_updater,       &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = veg_ns%cropseedn_deficit , &
            flux_out_grc_area = conv_nflux)
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootn_to_litter(p)     = -1._r8 * dwt_frootn_to_litter(p)
       dwt_livecrootn_to_litter(p) = -1._r8 * dwt_livecrootn_to_litter(p)
       dwt_deadcrootn_to_litter(p) = -1._r8 * dwt_deadcrootn_to_litter(p)
    end do

  end subroutine dyn_veg_ns_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_ns_Adjustments(bounds, clump_index, column_state_updater, col_ns)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
        !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_nitrogen_state)     , intent(inout) :: col_ns
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    real(r8)  :: col_slice(bounds%begc:bounds%endc)
    !character(len=*), parameter :: subname = 'NStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc
    col_ns%dyn_nbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
         col_slice = col_ns%decomp_npools_vr(begc:endc, j, l)
          call update_column_state_no_special_handling_acc(column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = col_slice,     &
               adjustment  = adjustment_one_level)

          col_ns%dyn_nbal_adjustments(begc:endc) = &
               col_ns%dyn_nbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

          col_ns%decomp_npools_vr(begc:endc, j, l) = col_slice
       end do
    end do

    do j = 1, nlevdecomp
      col_slice = col_ns%ntrunc_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc(column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,     &
            adjustment  = adjustment_one_level)

       col_ns%dyn_nbal_adjustments(begc:endc) = &
            col_ns%dyn_nbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       col_ns%ntrunc_vr(begc:endc,j) = col_slice

       col_slice = col_ns%sminn_vr(begc:endc, j)
       call update_column_state_no_special_handling_acc(column_state_updater, &
           bounds      = bounds                          , &
           clump_index = clump_index                     , &
           var         =  col_slice, &
           adjustment  = adjustment_one_level)

       col_ns%dyn_nbal_adjustments(begc:endc) = &
           col_ns%dyn_nbal_adjustments(begc:endc) + &
           adjustment_one_level(begc:endc) * dzsoi_decomp(j)

        col_ns%sminn_vr(begc:endc, j) = col_slice
    end do

  end subroutine dyn_col_ns_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_veg_ps_Adjustments(        &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafp_seed,                      &
       dwt_deadstemp_seed,                  &
       dwt_ppool_seed,                      &
       conv_pflux,                          &
       dwt_frootp_to_litter,                &
       dwt_livecrootp_to_litter,            &
       dwt_deadcrootp_to_litter,            &
       prod10_pflux,                        &
       prod100_pflux,                       &
       crop_product_pflux,                  &
       veg_ps                               &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: num_filterp_with_inactive
    integer                          , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)         , intent(in)    :: prior_weights
    type(patch_state_updater_type)   , intent(in)    :: patch_state_updater
    real(r8)                         , intent(inout) :: dwt_leafp_seed           (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_deadstemp_seed       (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_ppool_seed           (bounds%begp:)
    real(r8)                         , intent(inout) :: conv_pflux               (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_frootp_to_litter     (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_livecrootp_to_litter (bounds%begp:)
    real(r8)                         , intent(inout) :: dwt_deadcrootp_to_litter (bounds%begp:)
    real(r8)                         , intent(inout) :: prod10_pflux             (bounds%begp:)
    real(r8)                         , intent(inout) :: prod100_pflux            (bounds%begp:)
    real(r8)                         , intent(inout) :: crop_product_pflux       (bounds%begp:)
    type(vegetation_phosphorus_state), intent(inout) :: veg_ps
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafp_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemp_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_ppool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_pflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemp_patch_temp     (bounds%begp:bounds%endp)

    !character(len=*), parameter :: subname = 'dyn_veg_ps_Adjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp


    old_weight_was_zero = old_weight_was_zeroAcc(patch_state_updater,bounds)
    patch_grew          = patch_grewAcc(patch_state_updater,bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_P                        , &
         leaf_patch                 = veg_ps%leafp           , &
         leaf_storage_patch         = veg_ps%leafp_storage   , &
         leaf_xfer_patch            = veg_ps%leafp_xfer      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafp_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafp_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafp_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemp_patch(begp:endp)     , &
         pool_seed_param            = ppool_seed_param                    , &
         pool_seed_patch            = seed_ppool_patch(begp:endp))

    ! 1) LEAFP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = veg_ps%leafp    , &
         flux_out_grc_area = conv_pflux        , &
         seed              = seed_leafp_patch  , &
         seed_addition     = dwt_leafp_seed   )

    ! 2) LEAFP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = veg_ps%leafp_storage    , &
         flux_out_grc_area = conv_pflux                , &
         seed              = seed_leafp_storage_patch  , &
         seed_addition     = dwt_leafp_seed           )

    ! 3) LEAFP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = veg_ps%leafp_xfer   , &
         flux_out_grc_area = conv_pflux            , &
         seed              = seed_leafp_xfer_patch , &
         seed_addition     = dwt_leafp_seed        )

    ! 4) FROOTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%frootp             , &
         flux_out_col_area = dwt_frootp_to_litter)

    ! 5) FROOTP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%frootp_storage     , &
         flux_out_grc_area = conv_pflux)

    ! 6) FROOTP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%frootp_xfer        , &
         flux_out_grc_area = conv_pflux)

    ! 7) PROD10_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = veg_ps%deadstemp(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemp_patch_temp     , &
         flux1_out                  = prod10_pflux             , &
         flux2_out                  = wood_product_pflux       , &
         seed                       = seed_deadstemp_patch     )

    ! 8) PROD100_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = veg_ps%deadstemp(begp:endp)
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemp_patch_temp     , &
         flux1_out                  = prod100_pflux            , &
         flux2_out                  = wood_product_pflux       , &
         seed                       = seed_deadstemp_patch    )

    ! 9) DEADSTEMP_PATCH
    wood_product_pflux(begp:endp)      = 0._r8
    call update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = veg_ps%deadstemp       , &
         flux1_out                  = conv_pflux               , &
         flux2_out                  = wood_product_pflux       , &
         seed                       = seed_deadstemp_patch     , &
         seed_addition              = dwt_deadstemp_seed   )

    ! 10) DEADSTEMP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%deadstemp_storage  , &
         flux_out_grc_area = conv_pflux)

    ! 11) DEADSTEMP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%deadstemp_xfer     , &
         flux_out_grc_area = conv_pflux)

    ! 12) LIVESTEMP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%livestemp          , &
         flux_out_grc_area = conv_pflux)

    ! 13) LIVESTEMP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%livestemp_storage  , &
         flux_out_grc_area = conv_pflux)

    ! 14) LIVESTEMP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%livestemp_xfer     , &
         flux_out_grc_area = conv_pflux)

    ! 15) LIVECROOTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%livecrootp         , &
         flux_out_col_area = dwt_livecrootp_to_litter)

    ! 16) LIVECROOTP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%livecrootp_storage , &
         flux_out_grc_area = conv_pflux)

    ! 17) LIVECROOTP_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%livecrootp_xfer    , &
         flux_out_grc_area = conv_pflux)

    ! 18) DEADCROOTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%deadcrootp         , &
         flux_out_col_area = dwt_deadcrootp_to_litter)

    ! 19) DEADCROOTP_STORAGE_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%deadcrootp_storage , &
         flux_out_grc_area = conv_pflux)

    ! 20) DEADCROOT_XFER_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%deadcrootp_xfer    , &
         flux_out_grc_area = conv_pflux)

    ! 21) RETRANSP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%retransp           , &
         flux_out_grc_area = conv_pflux)

    ! 22) PTRUNC_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%ptrunc             , &
         flux_out_grc_area = conv_pflux)

    ! 23) PPOOL_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%ppool              , &
         flux_out_grc_area = conv_pflux)

    ! 24) DISPVEGP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%dispvegp)

    ! 25) STORVEGP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%storvegp)

    ! 26) TOTVEGP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%totvegp)

    ! 27) TOTPFTP_PATCH
    call update_patch_stateAcc(patch_state_updater,       &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = veg_ps%totpftp)

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call update_patch_stateAcc(patch_state_updater,       &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = veg_ps%cropseedp_deficit , &
            flux_out_grc_area = conv_pflux)
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootp_to_litter(p)     = -1._r8 * dwt_frootp_to_litter(p)
       dwt_livecrootp_to_litter(p) = -1._r8 * dwt_livecrootp_to_litter(p)
       dwt_deadcrootp_to_litter(p) = -1._r8 * dwt_deadcrootp_to_litter(p)
    end do

  end subroutine dyn_veg_ps_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_ps_Adjustments(bounds, clump_index, column_state_updater, col_ps)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
      !$acc routine seq
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_phosphorus_state)   , intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    integer           :: l, j
    integer           :: begc, endc
    real(r8)          :: adjustment_one_level(bounds%begc:bounds%endc)
    real(r8)  :: col_slice(bounds%begc:bounds%endc)
    !character(len=*), parameter :: subname = 'NStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    col_ps%dyn_pbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
         col_slice = col_ps%decomp_ppools_vr(begc:endc, j, l)
          call update_column_state_no_special_handling_acc( column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = col_slice,     &
               adjustment  = adjustment_one_level)

          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

          col_ps%decomp_ppools_vr(begc:endc, j, l) = col_slice
       end do
    end do

    do j = 1, nlevdecomp
      col_slice = col_ps%ptrunc_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,                &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       col_ps%ptrunc_vr(begc:endc,j) = col_slice
       !
       col_slice = col_ps%solutionp_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,             &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       col_ps%solutionp_vr(begc:endc,j) = col_slice
       !
       col_slice = col_ps%labilep_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,               &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       col_ps%labilep_vr(begc:endc,j) = col_slice
       !!
       col_slice = col_ps%secondp_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,               &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       col_ps%secondp_vr(begc:endc,j) = col_slice
    end do

  end subroutine dyn_col_ps_Adjustments

  !-----------------------------------------------------------------------
  subroutine dyn_col_ps_Adjustments_acc(bounds, clump_index, column_state_updater, col_ps)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !$acc routine seq
    ! !USES:
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    type(column_phosphorus_state)   , intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)
    real(r8)  :: col_slice(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------
    begc = bounds%begc
    endc = bounds%endc
    !allocate(col_slice(begc:endc))

    col_ps%dyn_pbal_adjustments(begc:endc) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp

         col_slice = col_ps%decomp_ppools_vr(begc:endc,j,l)

          call update_column_state_no_special_handling_acc( column_state_updater, &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = col_slice, &!col_ps%decomp_ppools_vr(begc:endc, j, l),     &
               adjustment  = adjustment_one_level )

          col_ps%dyn_pbal_adjustments(begc:endc) =      &
               col_ps%dyn_pbal_adjustments(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

          col_ps%decomp_ppools_vr(begc:endc,j,l) = col_slice
       end do
    end do

    do j = 1, nlevdecomp
      col_slice = col_ps%ptrunc_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,                &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       col_ps%ptrunc_vr(begc:endc,j) = col_slice
       !!
       col_slice = col_ps%solutionp_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,             &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       col_ps%solutionp_vr(begc:endc,j) = col_slice
        !!
       col_slice = col_ps%labilep_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,               &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
       col_ps%labilep_vr(begc:endc,j) = col_slice
       !!
       col_slice = col_ps%secondp_vr(begc:endc,j)
       call update_column_state_no_special_handling_acc( column_state_updater, &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = col_slice,               &
            adjustment  = adjustment_one_level)

       col_ps%dyn_pbal_adjustments(begc:endc) =      &
            col_ps%dyn_pbal_adjustments(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
        col_ps%secondp_vr(begc:endc,j) = col_slice
    end do
    !deallocate(col_slice)
  end subroutine dyn_col_ps_Adjustments_acc


end module dynSubgridAdjustmentsMod
