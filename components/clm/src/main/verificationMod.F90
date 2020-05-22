module verificationMod 
contains 
subroutine update_vars_UpdateDaylength(gpu)
     use GridcellType, only : grc_pp 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_UpdateDaylength.txt"
     else
          file='cpu_UpdateDaylength.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc grc_pp%dayl, & 
     !$acc grc_pp%prev_dayl )
     end if 
     !! CPU print statements !! 
     write(10,*) 'grc_pp%dayl' 
     write(10,*) grc_pp%dayl
     write(10,*) 'grc_pp%prev_dayl' 
     write(10,*) grc_pp%prev_dayl
     close(10)
end subroutine 
subroutine update_vars_clm_drv_init(gpu)
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : canopystate_vars 
     use clm_instMod, only : photosyns_vars 
     use ColumnDataType, only : col_ef 
     use ColumnDataType, only : col_wf 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_clm_drv_init.txt"
     else
          file='cpu_clm_drv_init.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%h2osno_old, & 
     !$acc col_ws%do_capsnow, & 
     !$acc col_ws%frac_iceold )
     !$acc update self(& 
     !$acc canopystate_vars%frac_veg_nosno_patch )
     !$acc update self(& 
     !$acc photosyns_vars%cisha_z_patch, & 
     !$acc photosyns_vars%cisun_z_patch )
     !$acc update self(& 
     !$acc col_ef%eflx_bot )
     !$acc update self(& 
     !$acc col_wf%qflx_glcice )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%h2osno_old' 
     write(10,*) col_ws%h2osno_old
     write(10,*) 'col_ws%do_capsnow' 
     write(10,*) col_ws%do_capsnow
     write(10,*) 'col_ws%frac_iceold' 
     write(10,*) col_ws%frac_iceold
     write(10,*) 'canopystate_vars%frac_veg_nosno_patch' 
     write(10,*) canopystate_vars%frac_veg_nosno_patch
     write(10,*) 'photosyns_vars%cisha_z_patch' 
     write(10,*) photosyns_vars%cisha_z_patch
     write(10,*) 'photosyns_vars%cisun_z_patch' 
     write(10,*) photosyns_vars%cisun_z_patch
     write(10,*) 'col_ef%eflx_bot' 
     write(10,*) col_ef%eflx_bot
     write(10,*) 'col_wf%qflx_glcice' 
     write(10,*) col_wf%qflx_glcice
     close(10)
end subroutine 
subroutine update_vars_downscale_forcings(gpu)
     use clm_instMod, only : atm2lnd_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_downscale_forcings.txt"
     else
          file='cpu_downscale_forcings.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc atm2lnd_vars%forc_rain_downscaled_col, & 
     !$acc atm2lnd_vars%forc_th_downscaled_col, & 
     !$acc atm2lnd_vars%forc_q_downscaled_col, & 
     !$acc atm2lnd_vars%forc_t_downscaled_col, & 
     !$acc atm2lnd_vars%forc_rho_downscaled_col, & 
     !$acc atm2lnd_vars%forc_pbot_downscaled_col, & 
     !$acc atm2lnd_vars%forc_snow_downscaled_col, & 
     !$acc atm2lnd_vars%forc_lwrad_downscaled_col, & 
     !$acc atm2lnd_vars%forc_lwrad_not_downscaled_grc, & 
     !$acc atm2lnd_vars%forc_lwrad_downscaled_col, & 
     !$acc atm2lnd_vars%forc_lwrad_not_downscaled_grc, & 
     !$acc atm2lnd_vars%forc_rain_downscaled_col, & 
     !$acc atm2lnd_vars%forc_th_downscaled_col, & 
     !$acc atm2lnd_vars%forc_q_downscaled_col, & 
     !$acc atm2lnd_vars%forc_t_downscaled_col, & 
     !$acc atm2lnd_vars%forc_rho_downscaled_col, & 
     !$acc atm2lnd_vars%forc_pbot_downscaled_col, & 
     !$acc atm2lnd_vars%forc_snow_downscaled_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'atm2lnd_vars%forc_rain_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_rain_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_th_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_th_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_q_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_q_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_t_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_t_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_rho_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_rho_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_pbot_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_pbot_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_snow_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_snow_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_lwrad_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_lwrad_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_lwrad_not_downscaled_grc' 
     write(10,*) atm2lnd_vars%forc_lwrad_not_downscaled_grc
     write(10,*) 'atm2lnd_vars%forc_lwrad_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_lwrad_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_lwrad_not_downscaled_grc' 
     write(10,*) atm2lnd_vars%forc_lwrad_not_downscaled_grc
     write(10,*) 'atm2lnd_vars%forc_rain_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_rain_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_th_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_th_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_q_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_q_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_t_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_t_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_rho_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_rho_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_pbot_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_pbot_downscaled_col
     write(10,*) 'atm2lnd_vars%forc_snow_downscaled_col' 
     write(10,*) atm2lnd_vars%forc_snow_downscaled_col
     close(10)
end subroutine 
subroutine update_vars_CanopyHydrology(gpu)
     use ColumnType, only : col_pp 
     use ColumnDataType, only : col_es 
     use VegetationDataType, only : veg_ws 
     use ColumnDataType, only : col_ws 
     use VegetationDataType, only : veg_wf 
     use ColumnDataType, only : col_wf 
     use clm_instMod, only : aerosol_vars 
     use clm_instMod, only : atm2lnd_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_CanopyHydrology.txt"
     else
          file='cpu_CanopyHydrology.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z )
     !$acc update self(& 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc veg_ws%h2ocan, & 
     !$acc veg_ws%fdry, & 
     !$acc veg_ws%fwet )
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%frac_iceold, & 
     !$acc col_ws%snow_depth, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%swe_old, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%int_snow, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%frac_sno, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%h2osfc, & 
     !$acc col_ws%frac_h2osfc, & 
     !$acc col_ws%frac_sno )
     !$acc update self(& 
     !$acc veg_wf%qflx_surf_irrig_patch, & 
     !$acc veg_wf%qflx_grnd_irrig_patch, & 
     !$acc veg_wf%qflx_prec_intr, & 
     !$acc veg_wf%n_irrig_steps_left, & 
     !$acc veg_wf%irrig_rate, & 
     !$acc veg_wf%qflx_rain_grnd, & 
     !$acc veg_wf%qflx_real_irrig_patch, & 
     !$acc veg_wf%qflx_over_supply_patch, & 
     !$acc veg_wf%qflx_irrig_patch, & 
     !$acc veg_wf%qflx_snwcp_ice, & 
     !$acc veg_wf%qflx_supply_patch, & 
     !$acc veg_wf%qflx_snwcp_liq, & 
     !$acc veg_wf%qflx_dirct_rain, & 
     !$acc veg_wf%qflx_leafdrip, & 
     !$acc veg_wf%qflx_snow_grnd, & 
     !$acc veg_wf%qflx_prec_grnd )
     !$acc update self(& 
     !$acc col_wf%qflx_snow_h2osfc, & 
     !$acc col_wf%qflx_floodc, & 
     !$acc col_wf%qflx_over_supply, & 
     !$acc col_wf%qflx_surf_irrig, & 
     !$acc col_wf%qflx_grnd_irrig, & 
     !$acc col_wf%qflx_snow_grnd )
     !$acc update self(& 
     !$acc aerosol_vars%mss_dsttot_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bc_top_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_octot_col, & 
     !$acc aerosol_vars%mss_dst_top_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst_col_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_oc_top_col, & 
     !$acc aerosol_vars%mss_bctot_col, & 
     !$acc aerosol_vars%mss_bc_col_col, & 
     !$acc aerosol_vars%mss_oc_col_col, & 
     !$acc aerosol_vars%mss_dst2_col )
     !$acc update self(& 
     !$acc atm2lnd_vars%supply_grc )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'veg_ws%h2ocan' 
     write(10,*) veg_ws%h2ocan
     write(10,*) 'veg_ws%fdry' 
     write(10,*) veg_ws%fdry
     write(10,*) 'veg_ws%fwet' 
     write(10,*) veg_ws%fwet
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%frac_iceold' 
     write(10,*) col_ws%frac_iceold
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%swe_old' 
     write(10,*) col_ws%swe_old
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%int_snow' 
     write(10,*) col_ws%int_snow
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%frac_sno' 
     write(10,*) col_ws%frac_sno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_ws%frac_h2osfc' 
     write(10,*) col_ws%frac_h2osfc
     write(10,*) 'col_ws%frac_sno' 
     write(10,*) col_ws%frac_sno
     write(10,*) 'veg_wf%qflx_surf_irrig_patch' 
     write(10,*) veg_wf%qflx_surf_irrig_patch
     write(10,*) 'veg_wf%qflx_grnd_irrig_patch' 
     write(10,*) veg_wf%qflx_grnd_irrig_patch
     write(10,*) 'veg_wf%qflx_prec_intr' 
     write(10,*) veg_wf%qflx_prec_intr
     write(10,*) 'veg_wf%n_irrig_steps_left' 
     write(10,*) veg_wf%n_irrig_steps_left
     write(10,*) 'veg_wf%irrig_rate' 
     write(10,*) veg_wf%irrig_rate
     write(10,*) 'veg_wf%qflx_rain_grnd' 
     write(10,*) veg_wf%qflx_rain_grnd
     write(10,*) 'veg_wf%qflx_real_irrig_patch' 
     write(10,*) veg_wf%qflx_real_irrig_patch
     write(10,*) 'veg_wf%qflx_over_supply_patch' 
     write(10,*) veg_wf%qflx_over_supply_patch
     write(10,*) 'veg_wf%qflx_irrig_patch' 
     write(10,*) veg_wf%qflx_irrig_patch
     write(10,*) 'veg_wf%qflx_snwcp_ice' 
     write(10,*) veg_wf%qflx_snwcp_ice
     write(10,*) 'veg_wf%qflx_supply_patch' 
     write(10,*) veg_wf%qflx_supply_patch
     write(10,*) 'veg_wf%qflx_snwcp_liq' 
     write(10,*) veg_wf%qflx_snwcp_liq
     write(10,*) 'veg_wf%qflx_dirct_rain' 
     write(10,*) veg_wf%qflx_dirct_rain
     write(10,*) 'veg_wf%qflx_leafdrip' 
     write(10,*) veg_wf%qflx_leafdrip
     write(10,*) 'veg_wf%qflx_snow_grnd' 
     write(10,*) veg_wf%qflx_snow_grnd
     write(10,*) 'veg_wf%qflx_prec_grnd' 
     write(10,*) veg_wf%qflx_prec_grnd
     write(10,*) 'col_wf%qflx_snow_h2osfc' 
     write(10,*) col_wf%qflx_snow_h2osfc
     write(10,*) 'col_wf%qflx_floodc' 
     write(10,*) col_wf%qflx_floodc
     write(10,*) 'col_wf%qflx_over_supply' 
     write(10,*) col_wf%qflx_over_supply
     write(10,*) 'col_wf%qflx_surf_irrig' 
     write(10,*) col_wf%qflx_surf_irrig
     write(10,*) 'col_wf%qflx_grnd_irrig' 
     write(10,*) col_wf%qflx_grnd_irrig
     write(10,*) 'col_wf%qflx_snow_grnd' 
     write(10,*) col_wf%qflx_snow_grnd
     write(10,*) 'aerosol_vars%mss_dsttot_col' 
     write(10,*) aerosol_vars%mss_dsttot_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bc_top_col' 
     write(10,*) aerosol_vars%mss_bc_top_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_octot_col' 
     write(10,*) aerosol_vars%mss_octot_col
     write(10,*) 'aerosol_vars%mss_dst_top_col' 
     write(10,*) aerosol_vars%mss_dst_top_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst_col_col' 
     write(10,*) aerosol_vars%mss_dst_col_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_oc_top_col' 
     write(10,*) aerosol_vars%mss_oc_top_col
     write(10,*) 'aerosol_vars%mss_bctot_col' 
     write(10,*) aerosol_vars%mss_bctot_col
     write(10,*) 'aerosol_vars%mss_bc_col_col' 
     write(10,*) aerosol_vars%mss_bc_col_col
     write(10,*) 'aerosol_vars%mss_oc_col_col' 
     write(10,*) aerosol_vars%mss_oc_col_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'atm2lnd_vars%supply_grc' 
     write(10,*) atm2lnd_vars%supply_grc
     close(10)
end subroutine 
subroutine update_vars_CanopySunShadeFractions(gpu)
     use clm_instMod, only : solarabs_vars 
     use clm_instMod, only : canopystate_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_CanopySunShadeFractions.txt"
     else
          file='cpu_CanopySunShadeFractions.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc solarabs_vars%parsha_z_patch, & 
     !$acc solarabs_vars%parsun_z_patch )
     !$acc update self(& 
     !$acc canopystate_vars%fsun_patch, & 
     !$acc canopystate_vars%laisun_patch, & 
     !$acc canopystate_vars%laisha_patch, & 
     !$acc canopystate_vars%laisha_z_patch, & 
     !$acc canopystate_vars%laisun_z_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'solarabs_vars%parsha_z_patch' 
     write(10,*) solarabs_vars%parsha_z_patch
     write(10,*) 'solarabs_vars%parsun_z_patch' 
     write(10,*) solarabs_vars%parsun_z_patch
     write(10,*) 'canopystate_vars%fsun_patch' 
     write(10,*) canopystate_vars%fsun_patch
     write(10,*) 'canopystate_vars%laisun_patch' 
     write(10,*) canopystate_vars%laisun_patch
     write(10,*) 'canopystate_vars%laisha_patch' 
     write(10,*) canopystate_vars%laisha_patch
     write(10,*) 'canopystate_vars%laisha_z_patch' 
     write(10,*) canopystate_vars%laisha_z_patch
     write(10,*) 'canopystate_vars%laisun_z_patch' 
     write(10,*) canopystate_vars%laisun_z_patch
     close(10)
end subroutine 
subroutine update_vars_SurfaceRadiation(gpu)
     use ColumnType, only : col_pp 
     use clm_instMod, only : canopystate_vars 
     use clm_instMod, only : solarabs_vars 
     use clm_instMod, only : surfrad_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SurfaceRadiation.txt"
     else
          file='cpu_SurfaceRadiation.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pp%snl )
     !$acc update self(& 
     !$acc canopystate_vars%fsun_patch )
     !$acc update self(& 
     !$acc solarabs_vars%fsr_nir_d_patch, & 
     !$acc solarabs_vars%fsa_r_patch, & 
     !$acc solarabs_vars%fsr_patch, & 
     !$acc solarabs_vars%fsr_nir_i_patch, & 
     !$acc solarabs_vars%fsr_nir_d_ln_patch, & 
     !$acc solarabs_vars%sabg_pen_patch, & 
     !$acc solarabs_vars%fsa_patch, & 
     !$acc solarabs_vars%sabv_patch, & 
     !$acc solarabs_vars%sabg_soil_patch, & 
     !$acc solarabs_vars%sabg_snow_patch, & 
     !$acc solarabs_vars%sub_surf_abs_SW_col, & 
     !$acc solarabs_vars%sabg_patch, & 
     !$acc solarabs_vars%fsds_nir_d_ln_patch, & 
     !$acc solarabs_vars%fsds_nir_d_patch, & 
     !$acc solarabs_vars%fsds_nir_i_patch, & 
     !$acc solarabs_vars%sabg_lyr_patch )
     !$acc update self(& 
     !$acc surfrad_vars%sfc_frc_oc_patch, & 
     !$acc surfrad_vars%fsr_vis_d_ln_patch, & 
     !$acc surfrad_vars%fsr_vis_d_patch, & 
     !$acc surfrad_vars%sfc_frc_dst_patch, & 
     !$acc surfrad_vars%fsds_sno_nd_patch, & 
     !$acc surfrad_vars%sfc_frc_bc_patch, & 
     !$acc surfrad_vars%fsds_vis_d_patch, & 
     !$acc surfrad_vars%fsr_sno_vd_patch, & 
     !$acc surfrad_vars%fsds_vis_d_ln_patch, & 
     !$acc surfrad_vars%fsds_sno_vd_patch, & 
     !$acc surfrad_vars%fsds_sno_vi_patch, & 
     !$acc surfrad_vars%sfc_frc_dst_sno_patch, & 
     !$acc surfrad_vars%fsds_vis_i_ln_patch, & 
     !$acc surfrad_vars%fsr_sno_ni_patch, & 
     !$acc surfrad_vars%fsr_vis_i_patch, & 
     !$acc surfrad_vars%sfc_frc_oc_sno_patch, & 
     !$acc surfrad_vars%fsds_sno_ni_patch, & 
     !$acc surfrad_vars%sfc_frc_aer_sno_patch, & 
     !$acc surfrad_vars%fsr_sno_vi_patch, & 
     !$acc surfrad_vars%fsr_sno_nd_patch, & 
     !$acc surfrad_vars%sfc_frc_aer_patch, & 
     !$acc surfrad_vars%fsds_vis_i_patch, & 
     !$acc surfrad_vars%parveg_ln_patch, & 
     !$acc surfrad_vars%sfc_frc_bc_sno_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'canopystate_vars%fsun_patch' 
     write(10,*) canopystate_vars%fsun_patch
     write(10,*) 'solarabs_vars%fsr_nir_d_patch' 
     write(10,*) solarabs_vars%fsr_nir_d_patch
     write(10,*) 'solarabs_vars%fsa_r_patch' 
     write(10,*) solarabs_vars%fsa_r_patch
     write(10,*) 'solarabs_vars%fsr_patch' 
     write(10,*) solarabs_vars%fsr_patch
     write(10,*) 'solarabs_vars%fsr_nir_i_patch' 
     write(10,*) solarabs_vars%fsr_nir_i_patch
     write(10,*) 'solarabs_vars%fsr_nir_d_ln_patch' 
     write(10,*) solarabs_vars%fsr_nir_d_ln_patch
     write(10,*) 'solarabs_vars%sabg_pen_patch' 
     write(10,*) solarabs_vars%sabg_pen_patch
     write(10,*) 'solarabs_vars%fsa_patch' 
     write(10,*) solarabs_vars%fsa_patch
     write(10,*) 'solarabs_vars%sabv_patch' 
     write(10,*) solarabs_vars%sabv_patch
     write(10,*) 'solarabs_vars%sabg_soil_patch' 
     write(10,*) solarabs_vars%sabg_soil_patch
     write(10,*) 'solarabs_vars%sabg_snow_patch' 
     write(10,*) solarabs_vars%sabg_snow_patch
     write(10,*) 'solarabs_vars%sub_surf_abs_SW_col' 
     write(10,*) solarabs_vars%sub_surf_abs_SW_col
     write(10,*) 'solarabs_vars%sabg_patch' 
     write(10,*) solarabs_vars%sabg_patch
     write(10,*) 'solarabs_vars%fsds_nir_d_ln_patch' 
     write(10,*) solarabs_vars%fsds_nir_d_ln_patch
     write(10,*) 'solarabs_vars%fsds_nir_d_patch' 
     write(10,*) solarabs_vars%fsds_nir_d_patch
     write(10,*) 'solarabs_vars%fsds_nir_i_patch' 
     write(10,*) solarabs_vars%fsds_nir_i_patch
     write(10,*) 'solarabs_vars%sabg_lyr_patch' 
     write(10,*) solarabs_vars%sabg_lyr_patch
     write(10,*) 'surfrad_vars%sfc_frc_oc_patch' 
     write(10,*) surfrad_vars%sfc_frc_oc_patch
     write(10,*) 'surfrad_vars%fsr_vis_d_ln_patch' 
     write(10,*) surfrad_vars%fsr_vis_d_ln_patch
     write(10,*) 'surfrad_vars%fsr_vis_d_patch' 
     write(10,*) surfrad_vars%fsr_vis_d_patch
     write(10,*) 'surfrad_vars%sfc_frc_dst_patch' 
     write(10,*) surfrad_vars%sfc_frc_dst_patch
     write(10,*) 'surfrad_vars%fsds_sno_nd_patch' 
     write(10,*) surfrad_vars%fsds_sno_nd_patch
     write(10,*) 'surfrad_vars%sfc_frc_bc_patch' 
     write(10,*) surfrad_vars%sfc_frc_bc_patch
     write(10,*) 'surfrad_vars%fsds_vis_d_patch' 
     write(10,*) surfrad_vars%fsds_vis_d_patch
     write(10,*) 'surfrad_vars%fsr_sno_vd_patch' 
     write(10,*) surfrad_vars%fsr_sno_vd_patch
     write(10,*) 'surfrad_vars%fsds_vis_d_ln_patch' 
     write(10,*) surfrad_vars%fsds_vis_d_ln_patch
     write(10,*) 'surfrad_vars%fsds_sno_vd_patch' 
     write(10,*) surfrad_vars%fsds_sno_vd_patch
     write(10,*) 'surfrad_vars%fsds_sno_vi_patch' 
     write(10,*) surfrad_vars%fsds_sno_vi_patch
     write(10,*) 'surfrad_vars%sfc_frc_dst_sno_patch' 
     write(10,*) surfrad_vars%sfc_frc_dst_sno_patch
     write(10,*) 'surfrad_vars%fsds_vis_i_ln_patch' 
     write(10,*) surfrad_vars%fsds_vis_i_ln_patch
     write(10,*) 'surfrad_vars%fsr_sno_ni_patch' 
     write(10,*) surfrad_vars%fsr_sno_ni_patch
     write(10,*) 'surfrad_vars%fsr_vis_i_patch' 
     write(10,*) surfrad_vars%fsr_vis_i_patch
     write(10,*) 'surfrad_vars%sfc_frc_oc_sno_patch' 
     write(10,*) surfrad_vars%sfc_frc_oc_sno_patch
     write(10,*) 'surfrad_vars%fsds_sno_ni_patch' 
     write(10,*) surfrad_vars%fsds_sno_ni_patch
     write(10,*) 'surfrad_vars%sfc_frc_aer_sno_patch' 
     write(10,*) surfrad_vars%sfc_frc_aer_sno_patch
     write(10,*) 'surfrad_vars%fsr_sno_vi_patch' 
     write(10,*) surfrad_vars%fsr_sno_vi_patch
     write(10,*) 'surfrad_vars%fsr_sno_nd_patch' 
     write(10,*) surfrad_vars%fsr_sno_nd_patch
     write(10,*) 'surfrad_vars%sfc_frc_aer_patch' 
     write(10,*) surfrad_vars%sfc_frc_aer_patch
     write(10,*) 'surfrad_vars%fsds_vis_i_patch' 
     write(10,*) surfrad_vars%fsds_vis_i_patch
     write(10,*) 'surfrad_vars%parveg_ln_patch' 
     write(10,*) surfrad_vars%parveg_ln_patch
     write(10,*) 'surfrad_vars%sfc_frc_bc_sno_patch' 
     write(10,*) surfrad_vars%sfc_frc_bc_sno_patch
     close(10)
end subroutine 
subroutine update_vars_UrbanRadiation(gpu)
     use LandunitType, only : lun_pp 
     use UrbanParamsType, only : urbanparams_vars 
     use clm_instMod, only : solarabs_vars 
     use VegetationDataType, only : veg_ef 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_UrbanRadiation.txt"
     else
          file='cpu_UrbanRadiation.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc lun_pp%wtroad_perv, & 
     !$acc lun_pp%canyon_hwr )
     !$acc update self(& 
     !$acc urbanparams_vars%em_wall )
     !$acc update self(& 
     !$acc solarabs_vars%sabs_shadewall_dir_lun, & 
     !$acc solarabs_vars%sabs_sunwall_dir_lun, & 
     !$acc solarabs_vars%sabs_shadewall_dif_lun, & 
     !$acc solarabs_vars%sabs_improad_dir_lun, & 
     !$acc solarabs_vars%sabs_perroad_dif_lun, & 
     !$acc solarabs_vars%fsa_patch, & 
     !$acc solarabs_vars%sabv_patch, & 
     !$acc solarabs_vars%sabg_patch, & 
     !$acc solarabs_vars%fsa_u_patch, & 
     !$acc solarabs_vars%sabs_sunwall_dif_lun, & 
     !$acc solarabs_vars%sabs_perroad_dir_lun, & 
     !$acc solarabs_vars%sabs_improad_dif_lun, & 
     !$acc solarabs_vars%sabs_roof_dif_lun, & 
     !$acc solarabs_vars%sabs_roof_dir_lun )
     !$acc update self(& 
     !$acc veg_ef%eflx_lwrad_net, & 
     !$acc veg_ef%eflx_lwrad_out, & 
     !$acc veg_ef%eflx_lwrad_net_u )
     end if 
     !! CPU print statements !! 
     write(10,*) 'lun_pp%wtroad_perv' 
     write(10,*) lun_pp%wtroad_perv
     write(10,*) 'lun_pp%canyon_hwr' 
     write(10,*) lun_pp%canyon_hwr
     write(10,*) 'urbanparams_vars%em_wall' 
     write(10,*) urbanparams_vars%em_wall
     write(10,*) 'solarabs_vars%sabs_shadewall_dir_lun' 
     write(10,*) solarabs_vars%sabs_shadewall_dir_lun
     write(10,*) 'solarabs_vars%sabs_sunwall_dir_lun' 
     write(10,*) solarabs_vars%sabs_sunwall_dir_lun
     write(10,*) 'solarabs_vars%sabs_shadewall_dif_lun' 
     write(10,*) solarabs_vars%sabs_shadewall_dif_lun
     write(10,*) 'solarabs_vars%sabs_improad_dir_lun' 
     write(10,*) solarabs_vars%sabs_improad_dir_lun
     write(10,*) 'solarabs_vars%sabs_perroad_dif_lun' 
     write(10,*) solarabs_vars%sabs_perroad_dif_lun
     write(10,*) 'solarabs_vars%fsa_patch' 
     write(10,*) solarabs_vars%fsa_patch
     write(10,*) 'solarabs_vars%sabv_patch' 
     write(10,*) solarabs_vars%sabv_patch
     write(10,*) 'solarabs_vars%sabg_patch' 
     write(10,*) solarabs_vars%sabg_patch
     write(10,*) 'solarabs_vars%fsa_u_patch' 
     write(10,*) solarabs_vars%fsa_u_patch
     write(10,*) 'solarabs_vars%sabs_sunwall_dif_lun' 
     write(10,*) solarabs_vars%sabs_sunwall_dif_lun
     write(10,*) 'solarabs_vars%sabs_perroad_dir_lun' 
     write(10,*) solarabs_vars%sabs_perroad_dir_lun
     write(10,*) 'solarabs_vars%sabs_improad_dif_lun' 
     write(10,*) solarabs_vars%sabs_improad_dif_lun
     write(10,*) 'solarabs_vars%sabs_roof_dif_lun' 
     write(10,*) solarabs_vars%sabs_roof_dif_lun
     write(10,*) 'solarabs_vars%sabs_roof_dir_lun' 
     write(10,*) solarabs_vars%sabs_roof_dir_lun
     write(10,*) 'veg_ef%eflx_lwrad_net' 
     write(10,*) veg_ef%eflx_lwrad_net
     write(10,*) 'veg_ef%eflx_lwrad_out' 
     write(10,*) veg_ef%eflx_lwrad_out
     write(10,*) 'veg_ef%eflx_lwrad_net_u' 
     write(10,*) veg_ef%eflx_lwrad_net_u
     close(10)
end subroutine 
subroutine update_vars_CanopyTemperature(gpu)
     use ColumnType, only : col_pp 
     use TopounitDataType, only : top_as 
     use ColumnDataType, only : col_ws 
     use VegetationDataType, only : veg_wf 
     use ColumnDataType, only : col_ef 
     use VegetationDataType, only : veg_ef 
     use clm_instMod, only : frictionvel_vars 
     use clm_instMod, only : canopystate_vars 
     use clm_instMod, only : soilstate_vars 
     use ColumnDataType, only : col_es 
     use VegetationDataType, only : veg_es 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_CanopyTemperature.txt"
     else
          file='cpu_CanopyTemperature.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pp%snl, & 
     !$acc col_pp%zii )
     !$acc update self(& 
     !$acc top_as%pbot )
     !$acc update self(& 
     !$acc col_ws%qg, & 
     !$acc col_ws%qg_soil, & 
     !$acc col_ws%dqgdT, & 
     !$acc col_ws%qg_h2osfc, & 
     !$acc col_ws%qg_snow )
     !$acc update self(& 
     !$acc veg_wf%qflx_evap_veg, & 
     !$acc veg_wf%qflx_tran_veg, & 
     !$acc veg_wf%qflx_evap_tot )
     !$acc update self(& 
     !$acc col_ef%htvp )
     !$acc update self(& 
     !$acc veg_ef%cgrnds, & 
     !$acc veg_ef%cgrnd, & 
     !$acc veg_ef%eflx_sh_tot, & 
     !$acc veg_ef%eflx_lh_tot, & 
     !$acc veg_ef%cgrndl, & 
     !$acc veg_ef%eflx_lh_tot_u, & 
     !$acc veg_ef%eflx_sh_veg, & 
     !$acc veg_ef%eflx_lh_tot_r, & 
     !$acc veg_ef%eflx_sh_tot_r, & 
     !$acc veg_ef%eflx_sh_tot_u )
     !$acc update self(& 
     !$acc frictionvel_vars%forc_hgt_u_patch, & 
     !$acc frictionvel_vars%z0qg_col, & 
     !$acc frictionvel_vars%z0mv_patch, & 
     !$acc frictionvel_vars%z0qv_patch, & 
     !$acc frictionvel_vars%forc_hgt_q_patch, & 
     !$acc frictionvel_vars%z0mg_col, & 
     !$acc frictionvel_vars%z0hg_col, & 
     !$acc frictionvel_vars%z0hv_patch, & 
     !$acc frictionvel_vars%forc_hgt_t_patch, & 
     !$acc frictionvel_vars%z0m_patch )
     !$acc update self(& 
     !$acc canopystate_vars%displa_patch )
     !$acc update self(& 
     !$acc soilstate_vars%soilalpha_u_col, & 
     !$acc soilstate_vars%rootr_road_perv_col, & 
     !$acc soilstate_vars%soilalpha_col, & 
     !$acc soilstate_vars%soilbeta_col )
     !$acc update self(& 
     !$acc col_es%t_h2osfc_bef, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_h2osfc, & 
     !$acc col_es%thv, & 
     !$acc col_es%emg, & 
     !$acc col_es%t_ssbef, & 
     !$acc col_es%t_grnd )
     !$acc update self(& 
     !$acc veg_es%thm, & 
     !$acc veg_es%emv )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zii' 
     write(10,*) col_pp%zii
     write(10,*) 'top_as%pbot' 
     write(10,*) top_as%pbot
     write(10,*) 'col_ws%qg' 
     write(10,*) col_ws%qg
     write(10,*) 'col_ws%qg_soil' 
     write(10,*) col_ws%qg_soil
     write(10,*) 'col_ws%dqgdT' 
     write(10,*) col_ws%dqgdT
     write(10,*) 'col_ws%qg_h2osfc' 
     write(10,*) col_ws%qg_h2osfc
     write(10,*) 'col_ws%qg_snow' 
     write(10,*) col_ws%qg_snow
     write(10,*) 'veg_wf%qflx_evap_veg' 
     write(10,*) veg_wf%qflx_evap_veg
     write(10,*) 'veg_wf%qflx_tran_veg' 
     write(10,*) veg_wf%qflx_tran_veg
     write(10,*) 'veg_wf%qflx_evap_tot' 
     write(10,*) veg_wf%qflx_evap_tot
     write(10,*) 'col_ef%htvp' 
     write(10,*) col_ef%htvp
     write(10,*) 'veg_ef%cgrnds' 
     write(10,*) veg_ef%cgrnds
     write(10,*) 'veg_ef%cgrnd' 
     write(10,*) veg_ef%cgrnd
     write(10,*) 'veg_ef%eflx_sh_tot' 
     write(10,*) veg_ef%eflx_sh_tot
     write(10,*) 'veg_ef%eflx_lh_tot' 
     write(10,*) veg_ef%eflx_lh_tot
     write(10,*) 'veg_ef%cgrndl' 
     write(10,*) veg_ef%cgrndl
     write(10,*) 'veg_ef%eflx_lh_tot_u' 
     write(10,*) veg_ef%eflx_lh_tot_u
     write(10,*) 'veg_ef%eflx_sh_veg' 
     write(10,*) veg_ef%eflx_sh_veg
     write(10,*) 'veg_ef%eflx_lh_tot_r' 
     write(10,*) veg_ef%eflx_lh_tot_r
     write(10,*) 'veg_ef%eflx_sh_tot_r' 
     write(10,*) veg_ef%eflx_sh_tot_r
     write(10,*) 'veg_ef%eflx_sh_tot_u' 
     write(10,*) veg_ef%eflx_sh_tot_u
     write(10,*) 'frictionvel_vars%forc_hgt_u_patch' 
     write(10,*) frictionvel_vars%forc_hgt_u_patch
     write(10,*) 'frictionvel_vars%z0qg_col' 
     write(10,*) frictionvel_vars%z0qg_col
     write(10,*) 'frictionvel_vars%z0mv_patch' 
     write(10,*) frictionvel_vars%z0mv_patch
     write(10,*) 'frictionvel_vars%z0qv_patch' 
     write(10,*) frictionvel_vars%z0qv_patch
     write(10,*) 'frictionvel_vars%forc_hgt_q_patch' 
     write(10,*) frictionvel_vars%forc_hgt_q_patch
     write(10,*) 'frictionvel_vars%z0mg_col' 
     write(10,*) frictionvel_vars%z0mg_col
     write(10,*) 'frictionvel_vars%z0hg_col' 
     write(10,*) frictionvel_vars%z0hg_col
     write(10,*) 'frictionvel_vars%z0hv_patch' 
     write(10,*) frictionvel_vars%z0hv_patch
     write(10,*) 'frictionvel_vars%forc_hgt_t_patch' 
     write(10,*) frictionvel_vars%forc_hgt_t_patch
     write(10,*) 'frictionvel_vars%z0m_patch' 
     write(10,*) frictionvel_vars%z0m_patch
     write(10,*) 'canopystate_vars%displa_patch' 
     write(10,*) canopystate_vars%displa_patch
     write(10,*) 'soilstate_vars%soilalpha_u_col' 
     write(10,*) soilstate_vars%soilalpha_u_col
     write(10,*) 'soilstate_vars%rootr_road_perv_col' 
     write(10,*) soilstate_vars%rootr_road_perv_col
     write(10,*) 'soilstate_vars%soilalpha_col' 
     write(10,*) soilstate_vars%soilalpha_col
     write(10,*) 'soilstate_vars%soilbeta_col' 
     write(10,*) soilstate_vars%soilbeta_col
     write(10,*) 'col_es%t_h2osfc_bef' 
     write(10,*) col_es%t_h2osfc_bef
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_h2osfc' 
     write(10,*) col_es%t_h2osfc
     write(10,*) 'col_es%thv' 
     write(10,*) col_es%thv
     write(10,*) 'col_es%emg' 
     write(10,*) col_es%emg
     write(10,*) 'col_es%t_ssbef' 
     write(10,*) col_es%t_ssbef
     write(10,*) 'col_es%t_grnd' 
     write(10,*) col_es%t_grnd
     write(10,*) 'veg_es%thm' 
     write(10,*) veg_es%thm
     write(10,*) 'veg_es%emv' 
     write(10,*) veg_es%emv
     close(10)
end subroutine 
subroutine update_vars_CanopyFluxes(gpu)
     use TopounitDataType, only : top_as 
     use clm_instMod, only : canopystate_vars 
     use clm_instMod, only : soilstate_vars 
     use clm_instMod, only : frictionvel_vars 
     use ColumnDataType, only : col_es 
     use VegetationDataType, only : veg_es 
     use ColumnDataType, only : col_ws 
     use VegetationDataType, only : veg_ws 
     use VegetationDataType, only : veg_wf 
     use clm_instMod, only : photosyns_vars 
     use clm_instMod, only : ch4_vars 
     use clm_instMod, only : energyflux_vars 
     use VegetationDataType, only : veg_ef 
     use clm_instMod, only : soilstate_vars 
     use clm_instMod, only : photosyns_vars 
     use clm_instMod, only : solarabs_vars 
     use clm_instMod, only : canopystate_vars 
     use clm_instMod, only : solarabs_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_CanopyFluxes.txt"
     else
          file='cpu_CanopyFluxes.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc top_as%pbot )
     !$acc update self(& 
     !$acc canopystate_vars%displa_patch, & 
     !$acc canopystate_vars%lbl_rsc_h2o_patch, & 
     !$acc canopystate_vars%dleaf_patch )
     !$acc update self(& 
     !$acc soilstate_vars%rootr_patch, & 
     !$acc soilstate_vars%watsat_col, & 
     !$acc soilstate_vars%eff_porosity_col, & 
     !$acc soilstate_vars%rootr_patch )
     !$acc update self(& 
     !$acc frictionvel_vars%ram1_patch, & 
     !$acc frictionvel_vars%z0mv_patch, & 
     !$acc frictionvel_vars%rb1_patch, & 
     !$acc frictionvel_vars%z0qv_patch, & 
     !$acc frictionvel_vars%z0hv_patch, & 
     !$acc frictionvel_vars%vds_patch, & 
     !$acc frictionvel_vars%u10_clm_patch, & 
     !$acc frictionvel_vars%va_patch, & 
     !$acc frictionvel_vars%u10_patch, & 
     !$acc frictionvel_vars%fv_patch )
     !$acc update self(& 
     !$acc col_es%thv )
     !$acc update self(& 
     !$acc veg_es%t_ref2m_r, & 
     !$acc veg_es%t_veg, & 
     !$acc veg_es%t_ref2m )
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osoi_liqvol )
     !$acc update self(& 
     !$acc veg_ws%h2ocan, & 
     !$acc veg_ws%rh_af, & 
     !$acc veg_ws%rh_ref2m, & 
     !$acc veg_ws%q_ref2m, & 
     !$acc veg_ws%rh_ref2m_r )
     !$acc update self(& 
     !$acc veg_wf%qflx_ev_snow, & 
     !$acc veg_wf%qflx_evap_soi, & 
     !$acc veg_wf%n_irrig_steps_left, & 
     !$acc veg_wf%irrig_rate, & 
     !$acc veg_wf%qflx_tran_veg, & 
     !$acc veg_wf%qflx_ev_soil, & 
     !$acc veg_wf%qflx_evap_veg, & 
     !$acc veg_wf%qflx_ev_h2osfc, & 
     !$acc veg_wf%qflx_tran_veg, & 
     !$acc veg_wf%qflx_tran_veg, & 
     !$acc veg_wf%qflx_tran_veg, & 
     !$acc veg_wf%qflx_tran_veg )
     !$acc update self(& 
     !$acc photosyns_vars%rssun_patch, & 
     !$acc photosyns_vars%rssha_patch, & 
     !$acc photosyns_vars%fpsn_wc_patch, & 
     !$acc photosyns_vars%psnsha_wp_patch, & 
     !$acc photosyns_vars%psnsun_wp_patch, & 
     !$acc photosyns_vars%psnsun_wj_patch, & 
     !$acc photosyns_vars%c13_psnsun_patch, & 
     !$acc photosyns_vars%rc13_psnsun_patch, & 
     !$acc photosyns_vars%alphapsnsha_patch, & 
     !$acc photosyns_vars%psnsun_wc_patch, & 
     !$acc photosyns_vars%rc13_psnsha_patch, & 
     !$acc photosyns_vars%psnsha_wc_patch, & 
     !$acc photosyns_vars%fpsn_wj_patch, & 
     !$acc photosyns_vars%alphapsnsun_patch, & 
     !$acc photosyns_vars%psnsha_wj_patch, & 
     !$acc photosyns_vars%rc13_canair_patch, & 
     !$acc photosyns_vars%fpsn_patch, & 
     !$acc photosyns_vars%psnsun_patch, & 
     !$acc photosyns_vars%fpsn_wp_patch, & 
     !$acc photosyns_vars%c13_psnsha_patch, & 
     !$acc photosyns_vars%c14_psnsha_patch, & 
     !$acc photosyns_vars%psnsha_patch, & 
     !$acc photosyns_vars%c14_psnsun_patch, & 
     !$acc photosyns_vars%ko_patch, & 
     !$acc photosyns_vars%gb_mol_patch, & 
     !$acc photosyns_vars%ac_patch, & 
     !$acc photosyns_vars%qe_patch, & 
     !$acc photosyns_vars%lmrsha_z_patch, & 
     !$acc photosyns_vars%ag_patch, & 
     !$acc photosyns_vars%psnsha_wp_patch, & 
     !$acc photosyns_vars%gs_mol_patch, & 
     !$acc photosyns_vars%rssha_patch, & 
     !$acc photosyns_vars%psnsun_wp_patch, & 
     !$acc photosyns_vars%psnsun_wj_patch, & 
     !$acc photosyns_vars%an_patch, & 
     !$acc photosyns_vars%alphapsnsha_patch, & 
     !$acc photosyns_vars%psnsun_wc_patch, & 
     !$acc photosyns_vars%lmrsun_z_patch, & 
     !$acc photosyns_vars%mbb_patch, & 
     !$acc photosyns_vars%psnsha_wc_patch, & 
     !$acc photosyns_vars%rh_leaf_patch, & 
     !$acc photosyns_vars%tpu_z_patch, & 
     !$acc photosyns_vars%kc_patch, & 
     !$acc photosyns_vars%bbb_patch, & 
     !$acc photosyns_vars%theta_cj_patch, & 
     !$acc photosyns_vars%alphapsnsun_patch, & 
     !$acc photosyns_vars%aj_patch, & 
     !$acc photosyns_vars%vcmax_z_patch, & 
     !$acc photosyns_vars%ap_patch, & 
     !$acc photosyns_vars%psnsun_z_patch, & 
     !$acc photosyns_vars%psnsha_wj_patch, & 
     !$acc photosyns_vars%rssun_z_patch, & 
     !$acc photosyns_vars%lmrsun_patch, & 
     !$acc photosyns_vars%cp_patch, & 
     !$acc photosyns_vars%c3flag_patch, & 
     !$acc photosyns_vars%cisun_z_patch, & 
     !$acc photosyns_vars%psnsun_patch, & 
     !$acc photosyns_vars%psnsha_z_patch, & 
     !$acc photosyns_vars%kp_z_patch, & 
     !$acc photosyns_vars%rssha_z_patch, & 
     !$acc photosyns_vars%rssun_patch, & 
     !$acc photosyns_vars%lmrsha_patch, & 
     !$acc photosyns_vars%psnsha_patch, & 
     !$acc photosyns_vars%cisha_z_patch, & 
     !$acc photosyns_vars%an_patch, & 
     !$acc photosyns_vars%aj_patch, & 
     !$acc photosyns_vars%ap_patch, & 
     !$acc photosyns_vars%ac_patch, & 
     !$acc photosyns_vars%ag_patch, & 
     !$acc photosyns_vars%an_patch, & 
     !$acc photosyns_vars%aj_patch, & 
     !$acc photosyns_vars%ap_patch, & 
     !$acc photosyns_vars%ac_patch, & 
     !$acc photosyns_vars%ag_patch, & 
     !$acc photosyns_vars%alphapsnsun_patch, & 
     !$acc photosyns_vars%alphapsnsha_patch, & 
     !$acc photosyns_vars%c13_psnsun_patch, & 
     !$acc photosyns_vars%rc13_psnsun_patch, & 
     !$acc photosyns_vars%fpsn_wp_patch, & 
     !$acc photosyns_vars%c13_psnsha_patch, & 
     !$acc photosyns_vars%c14_psnsha_patch, & 
     !$acc photosyns_vars%rc13_canair_patch, & 
     !$acc photosyns_vars%rc13_psnsha_patch, & 
     !$acc photosyns_vars%fpsn_patch, & 
     !$acc photosyns_vars%fpsn_wc_patch, & 
     !$acc photosyns_vars%c14_psnsun_patch, & 
     !$acc photosyns_vars%fpsn_wj_patch )
     !$acc update self(& 
     !$acc ch4_vars%grnd_ch4_cond_patch )
     !$acc update self(& 
     !$acc energyflux_vars%canopy_cond_patch, & 
     !$acc energyflux_vars%btran_patch, & 
     !$acc energyflux_vars%bsun_patch, & 
     !$acc energyflux_vars%bsha_patch, & 
     !$acc energyflux_vars%btran2_patch, & 
     !$acc energyflux_vars%rresis_patch, & 
     !$acc energyflux_vars%btran_patch, & 
     !$acc energyflux_vars%btran2_patch, & 
     !$acc energyflux_vars%rresis_patch )
     !$acc update self(& 
     !$acc veg_ef%cgrnds, & 
     !$acc veg_ef%cgrnd, & 
     !$acc veg_ef%eflx_sh_grnd, & 
     !$acc veg_ef%eflx_sh_soil, & 
     !$acc veg_ef%dlrad, & 
     !$acc veg_ef%taux, & 
     !$acc veg_ef%cgrndl, & 
     !$acc veg_ef%eflx_sh_snow, & 
     !$acc veg_ef%eflx_sh_h2osfc, & 
     !$acc veg_ef%eflx_sh_veg, & 
     !$acc veg_ef%tauy, & 
     !$acc veg_ef%ulrad )
     !$acc update self(& 
     !$acc soilstate_vars%soil_conductance_patch, & 
     !$acc soilstate_vars%root_conductance_patch, & 
     !$acc soilstate_vars%k_soil_root_patch )
     !$acc update self(& 
     !$acc photosyns_vars%ko_patch, & 
     !$acc photosyns_vars%vcmax_z_phs_patch, & 
     !$acc photosyns_vars%gb_mol_patch, & 
     !$acc photosyns_vars%qe_patch, & 
     !$acc photosyns_vars%lmrsha_z_patch, & 
     !$acc photosyns_vars%aj_phs_patch, & 
     !$acc photosyns_vars%psnsha_wp_patch, & 
     !$acc photosyns_vars%gs_mol_sun_patch, & 
     !$acc photosyns_vars%rssha_patch, & 
     !$acc photosyns_vars%psnsun_wp_patch, & 
     !$acc photosyns_vars%an_sun_patch, & 
     !$acc photosyns_vars%psnsun_wj_patch, & 
     !$acc photosyns_vars%alphapsnsha_patch, & 
     !$acc photosyns_vars%psnsun_wc_patch, & 
     !$acc photosyns_vars%lmrsun_z_patch, & 
     !$acc photosyns_vars%mbb_patch, & 
     !$acc photosyns_vars%an_sha_patch, & 
     !$acc photosyns_vars%psnsha_wc_patch, & 
     !$acc photosyns_vars%kc_patch, & 
     !$acc photosyns_vars%bbb_patch, & 
     !$acc photosyns_vars%theta_cj_patch, & 
     !$acc photosyns_vars%alphapsnsun_patch, & 
     !$acc photosyns_vars%tpu_z_phs_patch, & 
     !$acc photosyns_vars%psnsun_z_patch, & 
     !$acc photosyns_vars%psnsha_wj_patch, & 
     !$acc photosyns_vars%rssun_z_patch, & 
     !$acc photosyns_vars%kp_z_phs_patch, & 
     !$acc photosyns_vars%lmrsun_patch, & 
     !$acc photosyns_vars%cp_patch, & 
     !$acc photosyns_vars%ac_phs_patch, & 
     !$acc photosyns_vars%c3flag_patch, & 
     !$acc photosyns_vars%cisun_z_patch, & 
     !$acc photosyns_vars%gs_mol_sha_patch, & 
     !$acc photosyns_vars%ag_phs_patch, & 
     !$acc photosyns_vars%psnsun_patch, & 
     !$acc photosyns_vars%psnsha_z_patch, & 
     !$acc photosyns_vars%rssha_z_patch, & 
     !$acc photosyns_vars%rssun_patch, & 
     !$acc photosyns_vars%ap_phs_patch, & 
     !$acc photosyns_vars%lmrsha_patch, & 
     !$acc photosyns_vars%psnsha_patch, & 
     !$acc photosyns_vars%cisha_z_patch, & 
     !$acc photosyns_vars%an_sun_patch, & 
     !$acc photosyns_vars%ag_phs_patch, & 
     !$acc photosyns_vars%an_sha_patch, & 
     !$acc photosyns_vars%ap_phs_patch, & 
     !$acc photosyns_vars%aj_phs_patch, & 
     !$acc photosyns_vars%ac_phs_patch, & 
     !$acc photosyns_vars%an_sun_patch, & 
     !$acc photosyns_vars%ag_phs_patch, & 
     !$acc photosyns_vars%an_sha_patch, & 
     !$acc photosyns_vars%ap_phs_patch, & 
     !$acc photosyns_vars%aj_phs_patch, & 
     !$acc photosyns_vars%ac_phs_patch )
     !$acc update self(& 
     !$acc solarabs_vars%parsha_z_patch, & 
     !$acc solarabs_vars%parsun_z_patch )
     !$acc update self(& 
     !$acc canopystate_vars%vegwp_patch, & 
     !$acc canopystate_vars%vegwp_patch )
     !$acc update self(& 
     !$acc solarabs_vars%parsha_z_patch, & 
     !$acc solarabs_vars%parsun_z_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'top_as%pbot' 
     write(10,*) top_as%pbot
     write(10,*) 'canopystate_vars%displa_patch' 
     write(10,*) canopystate_vars%displa_patch
     write(10,*) 'canopystate_vars%lbl_rsc_h2o_patch' 
     write(10,*) canopystate_vars%lbl_rsc_h2o_patch
     write(10,*) 'canopystate_vars%dleaf_patch' 
     write(10,*) canopystate_vars%dleaf_patch
     write(10,*) 'soilstate_vars%rootr_patch' 
     write(10,*) soilstate_vars%rootr_patch
     write(10,*) 'soilstate_vars%watsat_col' 
     write(10,*) soilstate_vars%watsat_col
     write(10,*) 'soilstate_vars%eff_porosity_col' 
     write(10,*) soilstate_vars%eff_porosity_col
     write(10,*) 'soilstate_vars%rootr_patch' 
     write(10,*) soilstate_vars%rootr_patch
     write(10,*) 'frictionvel_vars%ram1_patch' 
     write(10,*) frictionvel_vars%ram1_patch
     write(10,*) 'frictionvel_vars%z0mv_patch' 
     write(10,*) frictionvel_vars%z0mv_patch
     write(10,*) 'frictionvel_vars%rb1_patch' 
     write(10,*) frictionvel_vars%rb1_patch
     write(10,*) 'frictionvel_vars%z0qv_patch' 
     write(10,*) frictionvel_vars%z0qv_patch
     write(10,*) 'frictionvel_vars%z0hv_patch' 
     write(10,*) frictionvel_vars%z0hv_patch
     write(10,*) 'frictionvel_vars%vds_patch' 
     write(10,*) frictionvel_vars%vds_patch
     write(10,*) 'frictionvel_vars%u10_clm_patch' 
     write(10,*) frictionvel_vars%u10_clm_patch
     write(10,*) 'frictionvel_vars%va_patch' 
     write(10,*) frictionvel_vars%va_patch
     write(10,*) 'frictionvel_vars%u10_patch' 
     write(10,*) frictionvel_vars%u10_patch
     write(10,*) 'frictionvel_vars%fv_patch' 
     write(10,*) frictionvel_vars%fv_patch
     write(10,*) 'col_es%thv' 
     write(10,*) col_es%thv
     write(10,*) 'veg_es%t_ref2m_r' 
     write(10,*) veg_es%t_ref2m_r
     write(10,*) 'veg_es%t_veg' 
     write(10,*) veg_es%t_veg
     write(10,*) 'veg_es%t_ref2m' 
     write(10,*) veg_es%t_ref2m
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osoi_liqvol' 
     write(10,*) col_ws%h2osoi_liqvol
     write(10,*) 'veg_ws%h2ocan' 
     write(10,*) veg_ws%h2ocan
     write(10,*) 'veg_ws%rh_af' 
     write(10,*) veg_ws%rh_af
     write(10,*) 'veg_ws%rh_ref2m' 
     write(10,*) veg_ws%rh_ref2m
     write(10,*) 'veg_ws%q_ref2m' 
     write(10,*) veg_ws%q_ref2m
     write(10,*) 'veg_ws%rh_ref2m_r' 
     write(10,*) veg_ws%rh_ref2m_r
     write(10,*) 'veg_wf%qflx_ev_snow' 
     write(10,*) veg_wf%qflx_ev_snow
     write(10,*) 'veg_wf%qflx_evap_soi' 
     write(10,*) veg_wf%qflx_evap_soi
     write(10,*) 'veg_wf%n_irrig_steps_left' 
     write(10,*) veg_wf%n_irrig_steps_left
     write(10,*) 'veg_wf%irrig_rate' 
     write(10,*) veg_wf%irrig_rate
     write(10,*) 'veg_wf%qflx_tran_veg' 
     write(10,*) veg_wf%qflx_tran_veg
     write(10,*) 'veg_wf%qflx_ev_soil' 
     write(10,*) veg_wf%qflx_ev_soil
     write(10,*) 'veg_wf%qflx_evap_veg' 
     write(10,*) veg_wf%qflx_evap_veg
     write(10,*) 'veg_wf%qflx_ev_h2osfc' 
     write(10,*) veg_wf%qflx_ev_h2osfc
     write(10,*) 'veg_wf%qflx_tran_veg' 
     write(10,*) veg_wf%qflx_tran_veg
     write(10,*) 'veg_wf%qflx_tran_veg' 
     write(10,*) veg_wf%qflx_tran_veg
     write(10,*) 'veg_wf%qflx_tran_veg' 
     write(10,*) veg_wf%qflx_tran_veg
     write(10,*) 'veg_wf%qflx_tran_veg' 
     write(10,*) veg_wf%qflx_tran_veg
     write(10,*) 'photosyns_vars%rssun_patch' 
     write(10,*) photosyns_vars%rssun_patch
     write(10,*) 'photosyns_vars%rssha_patch' 
     write(10,*) photosyns_vars%rssha_patch
     write(10,*) 'photosyns_vars%fpsn_wc_patch' 
     write(10,*) photosyns_vars%fpsn_wc_patch
     write(10,*) 'photosyns_vars%psnsha_wp_patch' 
     write(10,*) photosyns_vars%psnsha_wp_patch
     write(10,*) 'photosyns_vars%psnsun_wp_patch' 
     write(10,*) photosyns_vars%psnsun_wp_patch
     write(10,*) 'photosyns_vars%psnsun_wj_patch' 
     write(10,*) photosyns_vars%psnsun_wj_patch
     write(10,*) 'photosyns_vars%c13_psnsun_patch' 
     write(10,*) photosyns_vars%c13_psnsun_patch
     write(10,*) 'photosyns_vars%rc13_psnsun_patch' 
     write(10,*) photosyns_vars%rc13_psnsun_patch
     write(10,*) 'photosyns_vars%alphapsnsha_patch' 
     write(10,*) photosyns_vars%alphapsnsha_patch
     write(10,*) 'photosyns_vars%psnsun_wc_patch' 
     write(10,*) photosyns_vars%psnsun_wc_patch
     write(10,*) 'photosyns_vars%rc13_psnsha_patch' 
     write(10,*) photosyns_vars%rc13_psnsha_patch
     write(10,*) 'photosyns_vars%psnsha_wc_patch' 
     write(10,*) photosyns_vars%psnsha_wc_patch
     write(10,*) 'photosyns_vars%fpsn_wj_patch' 
     write(10,*) photosyns_vars%fpsn_wj_patch
     write(10,*) 'photosyns_vars%alphapsnsun_patch' 
     write(10,*) photosyns_vars%alphapsnsun_patch
     write(10,*) 'photosyns_vars%psnsha_wj_patch' 
     write(10,*) photosyns_vars%psnsha_wj_patch
     write(10,*) 'photosyns_vars%rc13_canair_patch' 
     write(10,*) photosyns_vars%rc13_canair_patch
     write(10,*) 'photosyns_vars%fpsn_patch' 
     write(10,*) photosyns_vars%fpsn_patch
     write(10,*) 'photosyns_vars%psnsun_patch' 
     write(10,*) photosyns_vars%psnsun_patch
     write(10,*) 'photosyns_vars%fpsn_wp_patch' 
     write(10,*) photosyns_vars%fpsn_wp_patch
     write(10,*) 'photosyns_vars%c13_psnsha_patch' 
     write(10,*) photosyns_vars%c13_psnsha_patch
     write(10,*) 'photosyns_vars%c14_psnsha_patch' 
     write(10,*) photosyns_vars%c14_psnsha_patch
     write(10,*) 'photosyns_vars%psnsha_patch' 
     write(10,*) photosyns_vars%psnsha_patch
     write(10,*) 'photosyns_vars%c14_psnsun_patch' 
     write(10,*) photosyns_vars%c14_psnsun_patch
     write(10,*) 'photosyns_vars%ko_patch' 
     write(10,*) photosyns_vars%ko_patch
     write(10,*) 'photosyns_vars%gb_mol_patch' 
     write(10,*) photosyns_vars%gb_mol_patch
     write(10,*) 'photosyns_vars%ac_patch' 
     write(10,*) photosyns_vars%ac_patch
     write(10,*) 'photosyns_vars%qe_patch' 
     write(10,*) photosyns_vars%qe_patch
     write(10,*) 'photosyns_vars%lmrsha_z_patch' 
     write(10,*) photosyns_vars%lmrsha_z_patch
     write(10,*) 'photosyns_vars%ag_patch' 
     write(10,*) photosyns_vars%ag_patch
     write(10,*) 'photosyns_vars%psnsha_wp_patch' 
     write(10,*) photosyns_vars%psnsha_wp_patch
     write(10,*) 'photosyns_vars%gs_mol_patch' 
     write(10,*) photosyns_vars%gs_mol_patch
     write(10,*) 'photosyns_vars%rssha_patch' 
     write(10,*) photosyns_vars%rssha_patch
     write(10,*) 'photosyns_vars%psnsun_wp_patch' 
     write(10,*) photosyns_vars%psnsun_wp_patch
     write(10,*) 'photosyns_vars%psnsun_wj_patch' 
     write(10,*) photosyns_vars%psnsun_wj_patch
     write(10,*) 'photosyns_vars%an_patch' 
     write(10,*) photosyns_vars%an_patch
     write(10,*) 'photosyns_vars%alphapsnsha_patch' 
     write(10,*) photosyns_vars%alphapsnsha_patch
     write(10,*) 'photosyns_vars%psnsun_wc_patch' 
     write(10,*) photosyns_vars%psnsun_wc_patch
     write(10,*) 'photosyns_vars%lmrsun_z_patch' 
     write(10,*) photosyns_vars%lmrsun_z_patch
     write(10,*) 'photosyns_vars%mbb_patch' 
     write(10,*) photosyns_vars%mbb_patch
     write(10,*) 'photosyns_vars%psnsha_wc_patch' 
     write(10,*) photosyns_vars%psnsha_wc_patch
     write(10,*) 'photosyns_vars%rh_leaf_patch' 
     write(10,*) photosyns_vars%rh_leaf_patch
     write(10,*) 'photosyns_vars%tpu_z_patch' 
     write(10,*) photosyns_vars%tpu_z_patch
     write(10,*) 'photosyns_vars%kc_patch' 
     write(10,*) photosyns_vars%kc_patch
     write(10,*) 'photosyns_vars%bbb_patch' 
     write(10,*) photosyns_vars%bbb_patch
     write(10,*) 'photosyns_vars%theta_cj_patch' 
     write(10,*) photosyns_vars%theta_cj_patch
     write(10,*) 'photosyns_vars%alphapsnsun_patch' 
     write(10,*) photosyns_vars%alphapsnsun_patch
     write(10,*) 'photosyns_vars%aj_patch' 
     write(10,*) photosyns_vars%aj_patch
     write(10,*) 'photosyns_vars%vcmax_z_patch' 
     write(10,*) photosyns_vars%vcmax_z_patch
     write(10,*) 'photosyns_vars%ap_patch' 
     write(10,*) photosyns_vars%ap_patch
     write(10,*) 'photosyns_vars%psnsun_z_patch' 
     write(10,*) photosyns_vars%psnsun_z_patch
     write(10,*) 'photosyns_vars%psnsha_wj_patch' 
     write(10,*) photosyns_vars%psnsha_wj_patch
     write(10,*) 'photosyns_vars%rssun_z_patch' 
     write(10,*) photosyns_vars%rssun_z_patch
     write(10,*) 'photosyns_vars%lmrsun_patch' 
     write(10,*) photosyns_vars%lmrsun_patch
     write(10,*) 'photosyns_vars%cp_patch' 
     write(10,*) photosyns_vars%cp_patch
     write(10,*) 'photosyns_vars%c3flag_patch' 
     write(10,*) photosyns_vars%c3flag_patch
     write(10,*) 'photosyns_vars%cisun_z_patch' 
     write(10,*) photosyns_vars%cisun_z_patch
     write(10,*) 'photosyns_vars%psnsun_patch' 
     write(10,*) photosyns_vars%psnsun_patch
     write(10,*) 'photosyns_vars%psnsha_z_patch' 
     write(10,*) photosyns_vars%psnsha_z_patch
     write(10,*) 'photosyns_vars%kp_z_patch' 
     write(10,*) photosyns_vars%kp_z_patch
     write(10,*) 'photosyns_vars%rssha_z_patch' 
     write(10,*) photosyns_vars%rssha_z_patch
     write(10,*) 'photosyns_vars%rssun_patch' 
     write(10,*) photosyns_vars%rssun_patch
     write(10,*) 'photosyns_vars%lmrsha_patch' 
     write(10,*) photosyns_vars%lmrsha_patch
     write(10,*) 'photosyns_vars%psnsha_patch' 
     write(10,*) photosyns_vars%psnsha_patch
     write(10,*) 'photosyns_vars%cisha_z_patch' 
     write(10,*) photosyns_vars%cisha_z_patch
     write(10,*) 'photosyns_vars%an_patch' 
     write(10,*) photosyns_vars%an_patch
     write(10,*) 'photosyns_vars%aj_patch' 
     write(10,*) photosyns_vars%aj_patch
     write(10,*) 'photosyns_vars%ap_patch' 
     write(10,*) photosyns_vars%ap_patch
     write(10,*) 'photosyns_vars%ac_patch' 
     write(10,*) photosyns_vars%ac_patch
     write(10,*) 'photosyns_vars%ag_patch' 
     write(10,*) photosyns_vars%ag_patch
     write(10,*) 'photosyns_vars%an_patch' 
     write(10,*) photosyns_vars%an_patch
     write(10,*) 'photosyns_vars%aj_patch' 
     write(10,*) photosyns_vars%aj_patch
     write(10,*) 'photosyns_vars%ap_patch' 
     write(10,*) photosyns_vars%ap_patch
     write(10,*) 'photosyns_vars%ac_patch' 
     write(10,*) photosyns_vars%ac_patch
     write(10,*) 'photosyns_vars%ag_patch' 
     write(10,*) photosyns_vars%ag_patch
     write(10,*) 'photosyns_vars%alphapsnsun_patch' 
     write(10,*) photosyns_vars%alphapsnsun_patch
     write(10,*) 'photosyns_vars%alphapsnsha_patch' 
     write(10,*) photosyns_vars%alphapsnsha_patch
     write(10,*) 'photosyns_vars%c13_psnsun_patch' 
     write(10,*) photosyns_vars%c13_psnsun_patch
     write(10,*) 'photosyns_vars%rc13_psnsun_patch' 
     write(10,*) photosyns_vars%rc13_psnsun_patch
     write(10,*) 'photosyns_vars%fpsn_wp_patch' 
     write(10,*) photosyns_vars%fpsn_wp_patch
     write(10,*) 'photosyns_vars%c13_psnsha_patch' 
     write(10,*) photosyns_vars%c13_psnsha_patch
     write(10,*) 'photosyns_vars%c14_psnsha_patch' 
     write(10,*) photosyns_vars%c14_psnsha_patch
     write(10,*) 'photosyns_vars%rc13_canair_patch' 
     write(10,*) photosyns_vars%rc13_canair_patch
     write(10,*) 'photosyns_vars%rc13_psnsha_patch' 
     write(10,*) photosyns_vars%rc13_psnsha_patch
     write(10,*) 'photosyns_vars%fpsn_patch' 
     write(10,*) photosyns_vars%fpsn_patch
     write(10,*) 'photosyns_vars%fpsn_wc_patch' 
     write(10,*) photosyns_vars%fpsn_wc_patch
     write(10,*) 'photosyns_vars%c14_psnsun_patch' 
     write(10,*) photosyns_vars%c14_psnsun_patch
     write(10,*) 'photosyns_vars%fpsn_wj_patch' 
     write(10,*) photosyns_vars%fpsn_wj_patch
     write(10,*) 'ch4_vars%grnd_ch4_cond_patch' 
     write(10,*) ch4_vars%grnd_ch4_cond_patch
     write(10,*) 'energyflux_vars%canopy_cond_patch' 
     write(10,*) energyflux_vars%canopy_cond_patch
     write(10,*) 'energyflux_vars%btran_patch' 
     write(10,*) energyflux_vars%btran_patch
     write(10,*) 'energyflux_vars%bsun_patch' 
     write(10,*) energyflux_vars%bsun_patch
     write(10,*) 'energyflux_vars%bsha_patch' 
     write(10,*) energyflux_vars%bsha_patch
     write(10,*) 'energyflux_vars%btran2_patch' 
     write(10,*) energyflux_vars%btran2_patch
     write(10,*) 'energyflux_vars%rresis_patch' 
     write(10,*) energyflux_vars%rresis_patch
     write(10,*) 'energyflux_vars%btran_patch' 
     write(10,*) energyflux_vars%btran_patch
     write(10,*) 'energyflux_vars%btran2_patch' 
     write(10,*) energyflux_vars%btran2_patch
     write(10,*) 'energyflux_vars%rresis_patch' 
     write(10,*) energyflux_vars%rresis_patch
     write(10,*) 'veg_ef%cgrnds' 
     write(10,*) veg_ef%cgrnds
     write(10,*) 'veg_ef%cgrnd' 
     write(10,*) veg_ef%cgrnd
     write(10,*) 'veg_ef%eflx_sh_grnd' 
     write(10,*) veg_ef%eflx_sh_grnd
     write(10,*) 'veg_ef%eflx_sh_soil' 
     write(10,*) veg_ef%eflx_sh_soil
     write(10,*) 'veg_ef%dlrad' 
     write(10,*) veg_ef%dlrad
     write(10,*) 'veg_ef%taux' 
     write(10,*) veg_ef%taux
     write(10,*) 'veg_ef%cgrndl' 
     write(10,*) veg_ef%cgrndl
     write(10,*) 'veg_ef%eflx_sh_snow' 
     write(10,*) veg_ef%eflx_sh_snow
     write(10,*) 'veg_ef%eflx_sh_h2osfc' 
     write(10,*) veg_ef%eflx_sh_h2osfc
     write(10,*) 'veg_ef%eflx_sh_veg' 
     write(10,*) veg_ef%eflx_sh_veg
     write(10,*) 'veg_ef%tauy' 
     write(10,*) veg_ef%tauy
     write(10,*) 'veg_ef%ulrad' 
     write(10,*) veg_ef%ulrad
     write(10,*) 'soilstate_vars%soil_conductance_patch' 
     write(10,*) soilstate_vars%soil_conductance_patch
     write(10,*) 'soilstate_vars%root_conductance_patch' 
     write(10,*) soilstate_vars%root_conductance_patch
     write(10,*) 'soilstate_vars%k_soil_root_patch' 
     write(10,*) soilstate_vars%k_soil_root_patch
     write(10,*) 'photosyns_vars%ko_patch' 
     write(10,*) photosyns_vars%ko_patch
     write(10,*) 'photosyns_vars%vcmax_z_phs_patch' 
     write(10,*) photosyns_vars%vcmax_z_phs_patch
     write(10,*) 'photosyns_vars%gb_mol_patch' 
     write(10,*) photosyns_vars%gb_mol_patch
     write(10,*) 'photosyns_vars%qe_patch' 
     write(10,*) photosyns_vars%qe_patch
     write(10,*) 'photosyns_vars%lmrsha_z_patch' 
     write(10,*) photosyns_vars%lmrsha_z_patch
     write(10,*) 'photosyns_vars%aj_phs_patch' 
     write(10,*) photosyns_vars%aj_phs_patch
     write(10,*) 'photosyns_vars%psnsha_wp_patch' 
     write(10,*) photosyns_vars%psnsha_wp_patch
     write(10,*) 'photosyns_vars%gs_mol_sun_patch' 
     write(10,*) photosyns_vars%gs_mol_sun_patch
     write(10,*) 'photosyns_vars%rssha_patch' 
     write(10,*) photosyns_vars%rssha_patch
     write(10,*) 'photosyns_vars%psnsun_wp_patch' 
     write(10,*) photosyns_vars%psnsun_wp_patch
     write(10,*) 'photosyns_vars%an_sun_patch' 
     write(10,*) photosyns_vars%an_sun_patch
     write(10,*) 'photosyns_vars%psnsun_wj_patch' 
     write(10,*) photosyns_vars%psnsun_wj_patch
     write(10,*) 'photosyns_vars%alphapsnsha_patch' 
     write(10,*) photosyns_vars%alphapsnsha_patch
     write(10,*) 'photosyns_vars%psnsun_wc_patch' 
     write(10,*) photosyns_vars%psnsun_wc_patch
     write(10,*) 'photosyns_vars%lmrsun_z_patch' 
     write(10,*) photosyns_vars%lmrsun_z_patch
     write(10,*) 'photosyns_vars%mbb_patch' 
     write(10,*) photosyns_vars%mbb_patch
     write(10,*) 'photosyns_vars%an_sha_patch' 
     write(10,*) photosyns_vars%an_sha_patch
     write(10,*) 'photosyns_vars%psnsha_wc_patch' 
     write(10,*) photosyns_vars%psnsha_wc_patch
     write(10,*) 'photosyns_vars%kc_patch' 
     write(10,*) photosyns_vars%kc_patch
     write(10,*) 'photosyns_vars%bbb_patch' 
     write(10,*) photosyns_vars%bbb_patch
     write(10,*) 'photosyns_vars%theta_cj_patch' 
     write(10,*) photosyns_vars%theta_cj_patch
     write(10,*) 'photosyns_vars%alphapsnsun_patch' 
     write(10,*) photosyns_vars%alphapsnsun_patch
     write(10,*) 'photosyns_vars%tpu_z_phs_patch' 
     write(10,*) photosyns_vars%tpu_z_phs_patch
     write(10,*) 'photosyns_vars%psnsun_z_patch' 
     write(10,*) photosyns_vars%psnsun_z_patch
     write(10,*) 'photosyns_vars%psnsha_wj_patch' 
     write(10,*) photosyns_vars%psnsha_wj_patch
     write(10,*) 'photosyns_vars%rssun_z_patch' 
     write(10,*) photosyns_vars%rssun_z_patch
     write(10,*) 'photosyns_vars%kp_z_phs_patch' 
     write(10,*) photosyns_vars%kp_z_phs_patch
     write(10,*) 'photosyns_vars%lmrsun_patch' 
     write(10,*) photosyns_vars%lmrsun_patch
     write(10,*) 'photosyns_vars%cp_patch' 
     write(10,*) photosyns_vars%cp_patch
     write(10,*) 'photosyns_vars%ac_phs_patch' 
     write(10,*) photosyns_vars%ac_phs_patch
     write(10,*) 'photosyns_vars%c3flag_patch' 
     write(10,*) photosyns_vars%c3flag_patch
     write(10,*) 'photosyns_vars%cisun_z_patch' 
     write(10,*) photosyns_vars%cisun_z_patch
     write(10,*) 'photosyns_vars%gs_mol_sha_patch' 
     write(10,*) photosyns_vars%gs_mol_sha_patch
     write(10,*) 'photosyns_vars%ag_phs_patch' 
     write(10,*) photosyns_vars%ag_phs_patch
     write(10,*) 'photosyns_vars%psnsun_patch' 
     write(10,*) photosyns_vars%psnsun_patch
     write(10,*) 'photosyns_vars%psnsha_z_patch' 
     write(10,*) photosyns_vars%psnsha_z_patch
     write(10,*) 'photosyns_vars%rssha_z_patch' 
     write(10,*) photosyns_vars%rssha_z_patch
     write(10,*) 'photosyns_vars%rssun_patch' 
     write(10,*) photosyns_vars%rssun_patch
     write(10,*) 'photosyns_vars%ap_phs_patch' 
     write(10,*) photosyns_vars%ap_phs_patch
     write(10,*) 'photosyns_vars%lmrsha_patch' 
     write(10,*) photosyns_vars%lmrsha_patch
     write(10,*) 'photosyns_vars%psnsha_patch' 
     write(10,*) photosyns_vars%psnsha_patch
     write(10,*) 'photosyns_vars%cisha_z_patch' 
     write(10,*) photosyns_vars%cisha_z_patch
     write(10,*) 'photosyns_vars%an_sun_patch' 
     write(10,*) photosyns_vars%an_sun_patch
     write(10,*) 'photosyns_vars%ag_phs_patch' 
     write(10,*) photosyns_vars%ag_phs_patch
     write(10,*) 'photosyns_vars%an_sha_patch' 
     write(10,*) photosyns_vars%an_sha_patch
     write(10,*) 'photosyns_vars%ap_phs_patch' 
     write(10,*) photosyns_vars%ap_phs_patch
     write(10,*) 'photosyns_vars%aj_phs_patch' 
     write(10,*) photosyns_vars%aj_phs_patch
     write(10,*) 'photosyns_vars%ac_phs_patch' 
     write(10,*) photosyns_vars%ac_phs_patch
     write(10,*) 'photosyns_vars%an_sun_patch' 
     write(10,*) photosyns_vars%an_sun_patch
     write(10,*) 'photosyns_vars%ag_phs_patch' 
     write(10,*) photosyns_vars%ag_phs_patch
     write(10,*) 'photosyns_vars%an_sha_patch' 
     write(10,*) photosyns_vars%an_sha_patch
     write(10,*) 'photosyns_vars%ap_phs_patch' 
     write(10,*) photosyns_vars%ap_phs_patch
     write(10,*) 'photosyns_vars%aj_phs_patch' 
     write(10,*) photosyns_vars%aj_phs_patch
     write(10,*) 'photosyns_vars%ac_phs_patch' 
     write(10,*) photosyns_vars%ac_phs_patch
     write(10,*) 'solarabs_vars%parsha_z_patch' 
     write(10,*) solarabs_vars%parsha_z_patch
     write(10,*) 'solarabs_vars%parsun_z_patch' 
     write(10,*) solarabs_vars%parsun_z_patch
     write(10,*) 'canopystate_vars%vegwp_patch' 
     write(10,*) canopystate_vars%vegwp_patch
     write(10,*) 'canopystate_vars%vegwp_patch' 
     write(10,*) canopystate_vars%vegwp_patch
     write(10,*) 'solarabs_vars%parsha_z_patch' 
     write(10,*) solarabs_vars%parsha_z_patch
     write(10,*) 'solarabs_vars%parsun_z_patch' 
     write(10,*) solarabs_vars%parsun_z_patch
     close(10)
end subroutine 
subroutine update_vars_BareGroundFluxes(gpu)
     use TopounitDataType, only : top_as 
     use ColumnDataType, only : col_es 
     use clm_instMod, only : ch4_vars 
     use VegetationDataType, only : veg_ef 
     use VegetationDataType, only : veg_es 
     use VegetationDataType, only : veg_ws 
     use clm_instMod, only : frictionvel_vars 
     use VegetationDataType, only : veg_wf 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_BareGroundFluxes.txt"
     else
          file='cpu_BareGroundFluxes.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc top_as%pbot )
     !$acc update self(& 
     !$acc col_es%thv )
     !$acc update self(& 
     !$acc ch4_vars%grnd_ch4_cond_patch )
     !$acc update self(& 
     !$acc veg_ef%eflx_sh_grnd, & 
     !$acc veg_ef%cgrnds, & 
     !$acc veg_ef%cgrnd, & 
     !$acc veg_ef%eflx_sh_tot, & 
     !$acc veg_ef%eflx_sh_soil, & 
     !$acc veg_ef%dlrad, & 
     !$acc veg_ef%taux, & 
     !$acc veg_ef%cgrndl, & 
     !$acc veg_ef%eflx_sh_snow, & 
     !$acc veg_ef%eflx_sh_h2osfc, & 
     !$acc veg_ef%tauy, & 
     !$acc veg_ef%ulrad )
     !$acc update self(& 
     !$acc veg_es%t_ref2m_r, & 
     !$acc veg_es%t_ref2m )
     !$acc update self(& 
     !$acc veg_ws%rh_ref2m_r, & 
     !$acc veg_ws%rh_ref2m, & 
     !$acc veg_ws%q_ref2m )
     !$acc update self(& 
     !$acc frictionvel_vars%ram1_patch, & 
     !$acc frictionvel_vars%vds_patch, & 
     !$acc frictionvel_vars%u10_clm_patch, & 
     !$acc frictionvel_vars%va_patch, & 
     !$acc frictionvel_vars%u10_patch, & 
     !$acc frictionvel_vars%fv_patch )
     !$acc update self(& 
     !$acc veg_wf%qflx_ev_snow, & 
     !$acc veg_wf%qflx_evap_soi, & 
     !$acc veg_wf%qflx_ev_soil, & 
     !$acc veg_wf%qflx_evap_tot, & 
     !$acc veg_wf%qflx_ev_h2osfc )
     end if 
     !! CPU print statements !! 
     write(10,*) 'top_as%pbot' 
     write(10,*) top_as%pbot
     write(10,*) 'col_es%thv' 
     write(10,*) col_es%thv
     write(10,*) 'ch4_vars%grnd_ch4_cond_patch' 
     write(10,*) ch4_vars%grnd_ch4_cond_patch
     write(10,*) 'veg_ef%eflx_sh_grnd' 
     write(10,*) veg_ef%eflx_sh_grnd
     write(10,*) 'veg_ef%cgrnds' 
     write(10,*) veg_ef%cgrnds
     write(10,*) 'veg_ef%cgrnd' 
     write(10,*) veg_ef%cgrnd
     write(10,*) 'veg_ef%eflx_sh_tot' 
     write(10,*) veg_ef%eflx_sh_tot
     write(10,*) 'veg_ef%eflx_sh_soil' 
     write(10,*) veg_ef%eflx_sh_soil
     write(10,*) 'veg_ef%dlrad' 
     write(10,*) veg_ef%dlrad
     write(10,*) 'veg_ef%taux' 
     write(10,*) veg_ef%taux
     write(10,*) 'veg_ef%cgrndl' 
     write(10,*) veg_ef%cgrndl
     write(10,*) 'veg_ef%eflx_sh_snow' 
     write(10,*) veg_ef%eflx_sh_snow
     write(10,*) 'veg_ef%eflx_sh_h2osfc' 
     write(10,*) veg_ef%eflx_sh_h2osfc
     write(10,*) 'veg_ef%tauy' 
     write(10,*) veg_ef%tauy
     write(10,*) 'veg_ef%ulrad' 
     write(10,*) veg_ef%ulrad
     write(10,*) 'veg_es%t_ref2m_r' 
     write(10,*) veg_es%t_ref2m_r
     write(10,*) 'veg_es%t_ref2m' 
     write(10,*) veg_es%t_ref2m
     write(10,*) 'veg_ws%rh_ref2m_r' 
     write(10,*) veg_ws%rh_ref2m_r
     write(10,*) 'veg_ws%rh_ref2m' 
     write(10,*) veg_ws%rh_ref2m
     write(10,*) 'veg_ws%q_ref2m' 
     write(10,*) veg_ws%q_ref2m
     write(10,*) 'frictionvel_vars%ram1_patch' 
     write(10,*) frictionvel_vars%ram1_patch
     write(10,*) 'frictionvel_vars%vds_patch' 
     write(10,*) frictionvel_vars%vds_patch
     write(10,*) 'frictionvel_vars%u10_clm_patch' 
     write(10,*) frictionvel_vars%u10_clm_patch
     write(10,*) 'frictionvel_vars%va_patch' 
     write(10,*) frictionvel_vars%va_patch
     write(10,*) 'frictionvel_vars%u10_patch' 
     write(10,*) frictionvel_vars%u10_patch
     write(10,*) 'frictionvel_vars%fv_patch' 
     write(10,*) frictionvel_vars%fv_patch
     write(10,*) 'veg_wf%qflx_ev_snow' 
     write(10,*) veg_wf%qflx_ev_snow
     write(10,*) 'veg_wf%qflx_evap_soi' 
     write(10,*) veg_wf%qflx_evap_soi
     write(10,*) 'veg_wf%qflx_ev_soil' 
     write(10,*) veg_wf%qflx_ev_soil
     write(10,*) 'veg_wf%qflx_evap_tot' 
     write(10,*) veg_wf%qflx_evap_tot
     write(10,*) 'veg_wf%qflx_ev_h2osfc' 
     write(10,*) veg_wf%qflx_ev_h2osfc
     close(10)
end subroutine 
subroutine update_vars_LakeFluxes(gpu)
     use TopounitDataType, only : top_as 
     use clm_instMod, only : solarabs_vars 
     use clm_instMod, only : frictionvel_vars 
     use VegetationDataType, only : veg_ws 
     use VegetationDataType, only : veg_wf 
     use VegetationDataType, only : veg_es 
     use ColumnDataType, only : col_es 
     use VegetationDataType, only : veg_ef 
     use clm_instMod, only : lakestate_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_LakeFluxes.txt"
     else
          file='cpu_LakeFluxes.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc top_as%pbot )
     !$acc update self(& 
     !$acc solarabs_vars%sabg_chk_patch )
     !$acc update self(& 
     !$acc frictionvel_vars%ram1_patch, & 
     !$acc frictionvel_vars%forc_hgt_u_patch, & 
     !$acc frictionvel_vars%z0qg_col, & 
     !$acc frictionvel_vars%z0mg_col, & 
     !$acc frictionvel_vars%forc_hgt_q_patch, & 
     !$acc frictionvel_vars%z0hg_col, & 
     !$acc frictionvel_vars%forc_hgt_t_patch, & 
     !$acc frictionvel_vars%vds_patch, & 
     !$acc frictionvel_vars%u10_clm_patch, & 
     !$acc frictionvel_vars%va_patch, & 
     !$acc frictionvel_vars%u10_patch, & 
     !$acc frictionvel_vars%fv_patch )
     !$acc update self(& 
     !$acc veg_ws%rh_ref2m, & 
     !$acc veg_ws%q_ref2m )
     !$acc update self(& 
     !$acc veg_wf%qflx_evap_soi, & 
     !$acc veg_wf%qflx_snwcp_ice, & 
     !$acc veg_wf%qflx_snwcp_liq, & 
     !$acc veg_wf%qflx_evap_tot, & 
     !$acc veg_wf%qflx_prec_grnd )
     !$acc update self(& 
     !$acc veg_es%t_veg, & 
     !$acc veg_es%t_ref2m )
     !$acc update self(& 
     !$acc col_es%t_grnd )
     !$acc update self(& 
     !$acc veg_ef%eflx_lwrad_net, & 
     !$acc veg_ef%eflx_sh_grnd, & 
     !$acc veg_ef%eflx_gnet, & 
     !$acc veg_ef%eflx_sh_tot, & 
     !$acc veg_ef%eflx_lh_tot, & 
     !$acc veg_ef%eflx_soil_grnd, & 
     !$acc veg_ef%taux, & 
     !$acc veg_ef%eflx_lwrad_out, & 
     !$acc veg_ef%tauy, & 
     !$acc veg_ef%eflx_lh_grnd )
     !$acc update self(& 
     !$acc lakestate_vars%ks_col, & 
     !$acc lakestate_vars%lake_raw_col, & 
     !$acc lakestate_vars%betaprime_col, & 
     !$acc lakestate_vars%ram1_lake_patch, & 
     !$acc lakestate_vars%ust_lake_col, & 
     !$acc lakestate_vars%ws_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'top_as%pbot' 
     write(10,*) top_as%pbot
     write(10,*) 'solarabs_vars%sabg_chk_patch' 
     write(10,*) solarabs_vars%sabg_chk_patch
     write(10,*) 'frictionvel_vars%ram1_patch' 
     write(10,*) frictionvel_vars%ram1_patch
     write(10,*) 'frictionvel_vars%forc_hgt_u_patch' 
     write(10,*) frictionvel_vars%forc_hgt_u_patch
     write(10,*) 'frictionvel_vars%z0qg_col' 
     write(10,*) frictionvel_vars%z0qg_col
     write(10,*) 'frictionvel_vars%z0mg_col' 
     write(10,*) frictionvel_vars%z0mg_col
     write(10,*) 'frictionvel_vars%forc_hgt_q_patch' 
     write(10,*) frictionvel_vars%forc_hgt_q_patch
     write(10,*) 'frictionvel_vars%z0hg_col' 
     write(10,*) frictionvel_vars%z0hg_col
     write(10,*) 'frictionvel_vars%forc_hgt_t_patch' 
     write(10,*) frictionvel_vars%forc_hgt_t_patch
     write(10,*) 'frictionvel_vars%vds_patch' 
     write(10,*) frictionvel_vars%vds_patch
     write(10,*) 'frictionvel_vars%u10_clm_patch' 
     write(10,*) frictionvel_vars%u10_clm_patch
     write(10,*) 'frictionvel_vars%va_patch' 
     write(10,*) frictionvel_vars%va_patch
     write(10,*) 'frictionvel_vars%u10_patch' 
     write(10,*) frictionvel_vars%u10_patch
     write(10,*) 'frictionvel_vars%fv_patch' 
     write(10,*) frictionvel_vars%fv_patch
     write(10,*) 'veg_ws%rh_ref2m' 
     write(10,*) veg_ws%rh_ref2m
     write(10,*) 'veg_ws%q_ref2m' 
     write(10,*) veg_ws%q_ref2m
     write(10,*) 'veg_wf%qflx_evap_soi' 
     write(10,*) veg_wf%qflx_evap_soi
     write(10,*) 'veg_wf%qflx_snwcp_ice' 
     write(10,*) veg_wf%qflx_snwcp_ice
     write(10,*) 'veg_wf%qflx_snwcp_liq' 
     write(10,*) veg_wf%qflx_snwcp_liq
     write(10,*) 'veg_wf%qflx_evap_tot' 
     write(10,*) veg_wf%qflx_evap_tot
     write(10,*) 'veg_wf%qflx_prec_grnd' 
     write(10,*) veg_wf%qflx_prec_grnd
     write(10,*) 'veg_es%t_veg' 
     write(10,*) veg_es%t_veg
     write(10,*) 'veg_es%t_ref2m' 
     write(10,*) veg_es%t_ref2m
     write(10,*) 'col_es%t_grnd' 
     write(10,*) col_es%t_grnd
     write(10,*) 'veg_ef%eflx_lwrad_net' 
     write(10,*) veg_ef%eflx_lwrad_net
     write(10,*) 'veg_ef%eflx_sh_grnd' 
     write(10,*) veg_ef%eflx_sh_grnd
     write(10,*) 'veg_ef%eflx_gnet' 
     write(10,*) veg_ef%eflx_gnet
     write(10,*) 'veg_ef%eflx_sh_tot' 
     write(10,*) veg_ef%eflx_sh_tot
     write(10,*) 'veg_ef%eflx_lh_tot' 
     write(10,*) veg_ef%eflx_lh_tot
     write(10,*) 'veg_ef%eflx_soil_grnd' 
     write(10,*) veg_ef%eflx_soil_grnd
     write(10,*) 'veg_ef%taux' 
     write(10,*) veg_ef%taux
     write(10,*) 'veg_ef%eflx_lwrad_out' 
     write(10,*) veg_ef%eflx_lwrad_out
     write(10,*) 'veg_ef%tauy' 
     write(10,*) veg_ef%tauy
     write(10,*) 'veg_ef%eflx_lh_grnd' 
     write(10,*) veg_ef%eflx_lh_grnd
     write(10,*) 'lakestate_vars%ks_col' 
     write(10,*) lakestate_vars%ks_col
     write(10,*) 'lakestate_vars%lake_raw_col' 
     write(10,*) lakestate_vars%lake_raw_col
     write(10,*) 'lakestate_vars%betaprime_col' 
     write(10,*) lakestate_vars%betaprime_col
     write(10,*) 'lakestate_vars%ram1_lake_patch' 
     write(10,*) lakestate_vars%ram1_lake_patch
     write(10,*) 'lakestate_vars%ust_lake_col' 
     write(10,*) lakestate_vars%ust_lake_col
     write(10,*) 'lakestate_vars%ws_col' 
     write(10,*) lakestate_vars%ws_col
     close(10)
end subroutine 
subroutine update_vars_DustEmission(gpu)
     use clm_instMod, only : dust_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_DustEmission.txt"
     else
          file='cpu_DustEmission.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc dust_vars%flx_mss_vrt_dst_patch, & 
     !$acc dust_vars%flx_mss_vrt_dst_tot_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'dust_vars%flx_mss_vrt_dst_patch' 
     write(10,*) dust_vars%flx_mss_vrt_dst_patch
     write(10,*) 'dust_vars%flx_mss_vrt_dst_tot_patch' 
     write(10,*) dust_vars%flx_mss_vrt_dst_tot_patch
     close(10)
end subroutine 
subroutine update_vars_DustDryDep(gpu)
     use clm_instMod, only : dust_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_DustDryDep.txt"
     else
          file='cpu_DustDryDep.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc dust_vars%vlc_trb_1_patch, & 
     !$acc dust_vars%vlc_trb_patch, & 
     !$acc dust_vars%vlc_trb_4_patch, & 
     !$acc dust_vars%vlc_trb_2_patch, & 
     !$acc dust_vars%vlc_trb_3_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'dust_vars%vlc_trb_1_patch' 
     write(10,*) dust_vars%vlc_trb_1_patch
     write(10,*) 'dust_vars%vlc_trb_patch' 
     write(10,*) dust_vars%vlc_trb_patch
     write(10,*) 'dust_vars%vlc_trb_4_patch' 
     write(10,*) dust_vars%vlc_trb_4_patch
     write(10,*) 'dust_vars%vlc_trb_2_patch' 
     write(10,*) dust_vars%vlc_trb_2_patch
     write(10,*) 'dust_vars%vlc_trb_3_patch' 
     write(10,*) dust_vars%vlc_trb_3_patch
     close(10)
end subroutine 
subroutine update_vars_LakeTemperature(gpu)
     use ColumnDataType, only : col_ws 
     use ColumnDataType, only : col_wf 
     use ColumnDataType, only : col_es 
     use clm_instMod, only : lakestate_vars 
     use clm_instMod, only : ch4_vars 
     use VegetationDataType, only : veg_ef 
     use ColumnDataType, only : col_ef 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_LakeTemperature.txt"
     else
          file='cpu_LakeTemperature.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%frac_iceold, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%snow_depth )
     !$acc update self(& 
     !$acc col_wf%qflx_snofrz, & 
     !$acc col_wf%qflx_snofrz_lyr, & 
     !$acc col_wf%qflx_snow_melt, & 
     !$acc col_wf%qflx_snofrz, & 
     !$acc col_wf%qflx_snomelt )
     !$acc update self(& 
     !$acc col_es%hc_soi, & 
     !$acc col_es%t_lake, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%hc_soisno, & 
     !$acc col_es%t_lake, & 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc lakestate_vars%lakeresist_col, & 
     !$acc lakestate_vars%betaprime_col, & 
     !$acc lakestate_vars%lake_icefrac_col, & 
     !$acc lakestate_vars%lake_icethick_col, & 
     !$acc lakestate_vars%savedtke1_col, & 
     !$acc lakestate_vars%lake_icefrac_col )
     !$acc update self(& 
     !$acc ch4_vars%grnd_ch4_cond_col )
     !$acc update self(& 
     !$acc veg_ef%eflx_soil_grnd, & 
     !$acc veg_ef%eflx_sh_grnd, & 
     !$acc veg_ef%eflx_sh_tot, & 
     !$acc veg_ef%eflx_gnet )
     !$acc update self(& 
     !$acc col_ef%errsoi, & 
     !$acc col_ef%eflx_snomelt, & 
     !$acc col_ef%imelt )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%frac_iceold' 
     write(10,*) col_ws%frac_iceold
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'col_wf%qflx_snofrz' 
     write(10,*) col_wf%qflx_snofrz
     write(10,*) 'col_wf%qflx_snofrz_lyr' 
     write(10,*) col_wf%qflx_snofrz_lyr
     write(10,*) 'col_wf%qflx_snow_melt' 
     write(10,*) col_wf%qflx_snow_melt
     write(10,*) 'col_wf%qflx_snofrz' 
     write(10,*) col_wf%qflx_snofrz
     write(10,*) 'col_wf%qflx_snomelt' 
     write(10,*) col_wf%qflx_snomelt
     write(10,*) 'col_es%hc_soi' 
     write(10,*) col_es%hc_soi
     write(10,*) 'col_es%t_lake' 
     write(10,*) col_es%t_lake
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%hc_soisno' 
     write(10,*) col_es%hc_soisno
     write(10,*) 'col_es%t_lake' 
     write(10,*) col_es%t_lake
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'lakestate_vars%lakeresist_col' 
     write(10,*) lakestate_vars%lakeresist_col
     write(10,*) 'lakestate_vars%betaprime_col' 
     write(10,*) lakestate_vars%betaprime_col
     write(10,*) 'lakestate_vars%lake_icefrac_col' 
     write(10,*) lakestate_vars%lake_icefrac_col
     write(10,*) 'lakestate_vars%lake_icethick_col' 
     write(10,*) lakestate_vars%lake_icethick_col
     write(10,*) 'lakestate_vars%savedtke1_col' 
     write(10,*) lakestate_vars%savedtke1_col
     write(10,*) 'lakestate_vars%lake_icefrac_col' 
     write(10,*) lakestate_vars%lake_icefrac_col
     write(10,*) 'ch4_vars%grnd_ch4_cond_col' 
     write(10,*) ch4_vars%grnd_ch4_cond_col
     write(10,*) 'veg_ef%eflx_soil_grnd' 
     write(10,*) veg_ef%eflx_soil_grnd
     write(10,*) 'veg_ef%eflx_sh_grnd' 
     write(10,*) veg_ef%eflx_sh_grnd
     write(10,*) 'veg_ef%eflx_sh_tot' 
     write(10,*) veg_ef%eflx_sh_tot
     write(10,*) 'veg_ef%eflx_gnet' 
     write(10,*) veg_ef%eflx_gnet
     write(10,*) 'col_ef%errsoi' 
     write(10,*) col_ef%errsoi
     write(10,*) 'col_ef%eflx_snomelt' 
     write(10,*) col_ef%eflx_snomelt
     write(10,*) 'col_ef%imelt' 
     write(10,*) col_ef%imelt
     close(10)
end subroutine 
subroutine update_vars_SoilTemperature(gpu)
     use ColumnDataType, only : col_ef 
     use ColumnDataType, only : col_es 
     use LandunitDataType, only : lun_es 
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : soilstate_vars 
     use VegetationDataType, only : veg_ef 
     use clm_instMod, only : solarabs_vars 
     use ColumnDataType, only : col_wf 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SoilTemperature.txt"
     else
          file='cpu_SoilTemperature.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ef%eflx_urban_heat, & 
     !$acc col_ef%eflx_urban_ac, & 
     !$acc col_ef%eflx_building_heat, & 
     !$acc col_ef%eflx_fgr12, & 
     !$acc col_ef%eflx_fgr, & 
     !$acc col_ef%xmf_h2osfc, & 
     !$acc col_ef%eflx_h2osfc_to_snow, & 
     !$acc col_ef%xmf_h2osfc, & 
     !$acc col_ef%imelt, & 
     !$acc col_ef%xmf, & 
     !$acc col_ef%eflx_snomelt_r, & 
     !$acc col_ef%eflx_snomelt_u, & 
     !$acc col_ef%eflx_snomelt )
     !$acc update self(& 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_h2osfc, & 
     !$acc col_es%hc_soisno, & 
     !$acc col_es%fact, & 
     !$acc col_es%hc_soi, & 
     !$acc col_es%c_h2osfc, & 
     !$acc col_es%t_grnd, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_h2osfc, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_h2osfc, & 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc lun_es%t_building )
     !$acc update self(& 
     !$acc col_ws%bw, & 
     !$acc col_ws%frac_h2osfc, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%frac_h2osfc, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%int_snow, & 
     !$acc col_ws%snow_depth, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%h2osfc, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%snow_depth )
     !$acc update self(& 
     !$acc soilstate_vars%thk_col )
     !$acc update self(& 
     !$acc veg_ef%eflx_gnet, & 
     !$acc veg_ef%dgnetdT, & 
     !$acc veg_ef%eflx_heat_from_ac, & 
     !$acc veg_ef%eflx_traffic, & 
     !$acc veg_ef%eflx_anthro, & 
     !$acc veg_ef%eflx_wasteheat )
     !$acc update self(& 
     !$acc solarabs_vars%sabg_chk_patch )
     !$acc update self(& 
     !$acc col_wf%qflx_h2osfc_to_ice, & 
     !$acc col_wf%qflx_snofrz_lyr, & 
     !$acc col_wf%qflx_snomelt, & 
     !$acc col_wf%qflx_snow_melt, & 
     !$acc col_wf%qflx_glcice, & 
     !$acc col_wf%qflx_glcice_melt, & 
     !$acc col_wf%qflx_snofrz )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ef%eflx_urban_heat' 
     write(10,*) col_ef%eflx_urban_heat
     write(10,*) 'col_ef%eflx_urban_ac' 
     write(10,*) col_ef%eflx_urban_ac
     write(10,*) 'col_ef%eflx_building_heat' 
     write(10,*) col_ef%eflx_building_heat
     write(10,*) 'col_ef%eflx_fgr12' 
     write(10,*) col_ef%eflx_fgr12
     write(10,*) 'col_ef%eflx_fgr' 
     write(10,*) col_ef%eflx_fgr
     write(10,*) 'col_ef%xmf_h2osfc' 
     write(10,*) col_ef%xmf_h2osfc
     write(10,*) 'col_ef%eflx_h2osfc_to_snow' 
     write(10,*) col_ef%eflx_h2osfc_to_snow
     write(10,*) 'col_ef%xmf_h2osfc' 
     write(10,*) col_ef%xmf_h2osfc
     write(10,*) 'col_ef%imelt' 
     write(10,*) col_ef%imelt
     write(10,*) 'col_ef%xmf' 
     write(10,*) col_ef%xmf
     write(10,*) 'col_ef%eflx_snomelt_r' 
     write(10,*) col_ef%eflx_snomelt_r
     write(10,*) 'col_ef%eflx_snomelt_u' 
     write(10,*) col_ef%eflx_snomelt_u
     write(10,*) 'col_ef%eflx_snomelt' 
     write(10,*) col_ef%eflx_snomelt
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_h2osfc' 
     write(10,*) col_es%t_h2osfc
     write(10,*) 'col_es%hc_soisno' 
     write(10,*) col_es%hc_soisno
     write(10,*) 'col_es%fact' 
     write(10,*) col_es%fact
     write(10,*) 'col_es%hc_soi' 
     write(10,*) col_es%hc_soi
     write(10,*) 'col_es%c_h2osfc' 
     write(10,*) col_es%c_h2osfc
     write(10,*) 'col_es%t_grnd' 
     write(10,*) col_es%t_grnd
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_h2osfc' 
     write(10,*) col_es%t_h2osfc
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_h2osfc' 
     write(10,*) col_es%t_h2osfc
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'lun_es%t_building' 
     write(10,*) lun_es%t_building
     write(10,*) 'col_ws%bw' 
     write(10,*) col_ws%bw
     write(10,*) 'col_ws%frac_h2osfc' 
     write(10,*) col_ws%frac_h2osfc
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%frac_h2osfc' 
     write(10,*) col_ws%frac_h2osfc
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%int_snow' 
     write(10,*) col_ws%int_snow
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'soilstate_vars%thk_col' 
     write(10,*) soilstate_vars%thk_col
     write(10,*) 'veg_ef%eflx_gnet' 
     write(10,*) veg_ef%eflx_gnet
     write(10,*) 'veg_ef%dgnetdT' 
     write(10,*) veg_ef%dgnetdT
     write(10,*) 'veg_ef%eflx_heat_from_ac' 
     write(10,*) veg_ef%eflx_heat_from_ac
     write(10,*) 'veg_ef%eflx_traffic' 
     write(10,*) veg_ef%eflx_traffic
     write(10,*) 'veg_ef%eflx_anthro' 
     write(10,*) veg_ef%eflx_anthro
     write(10,*) 'veg_ef%eflx_wasteheat' 
     write(10,*) veg_ef%eflx_wasteheat
     write(10,*) 'solarabs_vars%sabg_chk_patch' 
     write(10,*) solarabs_vars%sabg_chk_patch
     write(10,*) 'col_wf%qflx_h2osfc_to_ice' 
     write(10,*) col_wf%qflx_h2osfc_to_ice
     write(10,*) 'col_wf%qflx_snofrz_lyr' 
     write(10,*) col_wf%qflx_snofrz_lyr
     write(10,*) 'col_wf%qflx_snomelt' 
     write(10,*) col_wf%qflx_snomelt
     write(10,*) 'col_wf%qflx_snow_melt' 
     write(10,*) col_wf%qflx_snow_melt
     write(10,*) 'col_wf%qflx_glcice' 
     write(10,*) col_wf%qflx_glcice
     write(10,*) 'col_wf%qflx_glcice_melt' 
     write(10,*) col_wf%qflx_glcice_melt
     write(10,*) 'col_wf%qflx_snofrz' 
     write(10,*) col_wf%qflx_snofrz
     close(10)
end subroutine 
subroutine update_vars_SoilFluxes(gpu)
     use VegetationDataType, only : veg_wf 
     use VegetationDataType, only : veg_ef 
     use ColumnDataType, only : col_ef 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SoilFluxes.txt"
     else
          file='cpu_SoilFluxes.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc veg_wf%qflx_ev_snow, & 
     !$acc veg_wf%qflx_evap_soi, & 
     !$acc veg_wf%qflx_evap_can, & 
     !$acc veg_wf%qflx_snwcp_ice, & 
     !$acc veg_wf%qflx_sub_snow, & 
     !$acc veg_wf%qflx_snwcp_liq, & 
     !$acc veg_wf%qflx_evap_grnd, & 
     !$acc veg_wf%qflx_ev_soil, & 
     !$acc veg_wf%qflx_dew_grnd, & 
     !$acc veg_wf%qflx_evap_tot, & 
     !$acc veg_wf%qflx_dew_snow, & 
     !$acc veg_wf%qflx_ev_h2osfc )
     !$acc update self(& 
     !$acc veg_ef%eflx_lwrad_net, & 
     !$acc veg_ef%eflx_lwrad_net_u, & 
     !$acc veg_ef%eflx_sh_grnd, & 
     !$acc veg_ef%eflx_lh_vegt, & 
     !$acc veg_ef%eflx_sh_tot, & 
     !$acc veg_ef%eflx_lh_tot, & 
     !$acc veg_ef%eflx_soil_grnd, & 
     !$acc veg_ef%eflx_lwrad_out_r, & 
     !$acc veg_ef%eflx_lh_grnd, & 
     !$acc veg_ef%eflx_soil_grnd_u, & 
     !$acc veg_ef%eflx_lh_tot_u, & 
     !$acc veg_ef%eflx_lwrad_net_r, & 
     !$acc veg_ef%errsoi, & 
     !$acc veg_ef%eflx_lwrad_out_u, & 
     !$acc veg_ef%eflx_soil_grnd_r, & 
     !$acc veg_ef%eflx_lh_tot_r, & 
     !$acc veg_ef%eflx_lwrad_out, & 
     !$acc veg_ef%eflx_lh_vege, & 
     !$acc veg_ef%eflx_sh_tot_r, & 
     !$acc veg_ef%eflx_sh_tot_u )
     !$acc update self(& 
     !$acc col_ef%errsoi )
     end if 
     !! CPU print statements !! 
     write(10,*) 'veg_wf%qflx_ev_snow' 
     write(10,*) veg_wf%qflx_ev_snow
     write(10,*) 'veg_wf%qflx_evap_soi' 
     write(10,*) veg_wf%qflx_evap_soi
     write(10,*) 'veg_wf%qflx_evap_can' 
     write(10,*) veg_wf%qflx_evap_can
     write(10,*) 'veg_wf%qflx_snwcp_ice' 
     write(10,*) veg_wf%qflx_snwcp_ice
     write(10,*) 'veg_wf%qflx_sub_snow' 
     write(10,*) veg_wf%qflx_sub_snow
     write(10,*) 'veg_wf%qflx_snwcp_liq' 
     write(10,*) veg_wf%qflx_snwcp_liq
     write(10,*) 'veg_wf%qflx_evap_grnd' 
     write(10,*) veg_wf%qflx_evap_grnd
     write(10,*) 'veg_wf%qflx_ev_soil' 
     write(10,*) veg_wf%qflx_ev_soil
     write(10,*) 'veg_wf%qflx_dew_grnd' 
     write(10,*) veg_wf%qflx_dew_grnd
     write(10,*) 'veg_wf%qflx_evap_tot' 
     write(10,*) veg_wf%qflx_evap_tot
     write(10,*) 'veg_wf%qflx_dew_snow' 
     write(10,*) veg_wf%qflx_dew_snow
     write(10,*) 'veg_wf%qflx_ev_h2osfc' 
     write(10,*) veg_wf%qflx_ev_h2osfc
     write(10,*) 'veg_ef%eflx_lwrad_net' 
     write(10,*) veg_ef%eflx_lwrad_net
     write(10,*) 'veg_ef%eflx_lwrad_net_u' 
     write(10,*) veg_ef%eflx_lwrad_net_u
     write(10,*) 'veg_ef%eflx_sh_grnd' 
     write(10,*) veg_ef%eflx_sh_grnd
     write(10,*) 'veg_ef%eflx_lh_vegt' 
     write(10,*) veg_ef%eflx_lh_vegt
     write(10,*) 'veg_ef%eflx_sh_tot' 
     write(10,*) veg_ef%eflx_sh_tot
     write(10,*) 'veg_ef%eflx_lh_tot' 
     write(10,*) veg_ef%eflx_lh_tot
     write(10,*) 'veg_ef%eflx_soil_grnd' 
     write(10,*) veg_ef%eflx_soil_grnd
     write(10,*) 'veg_ef%eflx_lwrad_out_r' 
     write(10,*) veg_ef%eflx_lwrad_out_r
     write(10,*) 'veg_ef%eflx_lh_grnd' 
     write(10,*) veg_ef%eflx_lh_grnd
     write(10,*) 'veg_ef%eflx_soil_grnd_u' 
     write(10,*) veg_ef%eflx_soil_grnd_u
     write(10,*) 'veg_ef%eflx_lh_tot_u' 
     write(10,*) veg_ef%eflx_lh_tot_u
     write(10,*) 'veg_ef%eflx_lwrad_net_r' 
     write(10,*) veg_ef%eflx_lwrad_net_r
     write(10,*) 'veg_ef%errsoi' 
     write(10,*) veg_ef%errsoi
     write(10,*) 'veg_ef%eflx_lwrad_out_u' 
     write(10,*) veg_ef%eflx_lwrad_out_u
     write(10,*) 'veg_ef%eflx_soil_grnd_r' 
     write(10,*) veg_ef%eflx_soil_grnd_r
     write(10,*) 'veg_ef%eflx_lh_tot_r' 
     write(10,*) veg_ef%eflx_lh_tot_r
     write(10,*) 'veg_ef%eflx_lwrad_out' 
     write(10,*) veg_ef%eflx_lwrad_out
     write(10,*) 'veg_ef%eflx_lh_vege' 
     write(10,*) veg_ef%eflx_lh_vege
     write(10,*) 'veg_ef%eflx_sh_tot_r' 
     write(10,*) veg_ef%eflx_sh_tot_r
     write(10,*) 'veg_ef%eflx_sh_tot_u' 
     write(10,*) veg_ef%eflx_sh_tot_u
     write(10,*) 'col_ef%errsoi' 
     write(10,*) col_ef%errsoi
     close(10)
end subroutine 
subroutine update_vars_clm_drv_patch2col(gpu)
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_clm_drv_patch2col.txt"
     else
          file='cpu_clm_drv_patch2col.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     end if 
     !! CPU print statements !! 
     close(10)
end subroutine 
subroutine update_vars_HydrologyNoDrainage(gpu)
     use ColumnType, only : col_pp 
     use ColumnDataType, only : col_es 
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : soilstate_vars 
     use ColumnDataType, only : col_wf 
     use clm_instMod, only : aerosol_vars 
     use clm_instMod, only : soilhydrology_vars 
     use VegetationDataType, only : veg_wf 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_HydrologyNoDrainage.txt"
     else
          file='cpu_HydrologyNoDrainage.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pp%dz, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z )
     !$acc update self(& 
     !$acc col_es%t_soi17cm, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%dTdz_top, & 
     !$acc col_es%snot_top, & 
     !$acc col_es%t_grnd_r, & 
     !$acc col_es%t_soi10cm, & 
     !$acc col_es%t_grnd_u, & 
     !$acc col_es%t_grnd, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%wf, & 
     !$acc col_ws%snow_persistence, & 
     !$acc col_ws%h2osoi_vol, & 
     !$acc col_ws%h2osno_top, & 
     !$acc col_ws%snw_rds_top, & 
     !$acc col_ws%snowliq, & 
     !$acc col_ws%h2osoi_liqice_10cm, & 
     !$acc col_ws%wf2, & 
     !$acc col_ws%air_vol, & 
     !$acc col_ws%h2osoi_liqvol, & 
     !$acc col_ws%snowice, & 
     !$acc col_ws%sno_liq_top, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%snowdp, & 
     !$acc col_ws%h2osoi_icevol, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%int_snow, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osfc, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_vol, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osoi_vol, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%snow_depth, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%int_snow, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%frac_sno, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%h2osoi_ice )
     !$acc update self(& 
     !$acc soilstate_vars%soilpsi_col, & 
     !$acc soilstate_vars%eff_porosity_col, & 
     !$acc soilstate_vars%eff_porosity_col, & 
     !$acc soilstate_vars%rootr_col, & 
     !$acc soilstate_vars%rootr_col, & 
     !$acc soilstate_vars%hk_l_col, & 
     !$acc soilstate_vars%smp_l_col )
     !$acc update self(& 
     !$acc col_wf%qflx_snow_melt, & 
     !$acc col_wf%qflx_top_soil, & 
     !$acc col_wf%mflx_neg_snow_1d, & 
     !$acc col_wf%qflx_top_soil, & 
     !$acc col_wf%qflx_surf, & 
     !$acc col_wf%qflx_infl, & 
     !$acc col_wf%qflx_h2osfc_surf, & 
     !$acc col_wf%qflx_gross_evap_soil, & 
     !$acc col_wf%qflx_surf, & 
     !$acc col_wf%qflx_gross_infl_soil, & 
     !$acc col_wf%qflx_drain_perched, & 
     !$acc col_wf%qflx_qrgwl, & 
     !$acc col_wf%mflx_drain_perched_1d, & 
     !$acc col_wf%qflx_drain, & 
     !$acc col_wf%qflx_rsub_sat, & 
     !$acc col_wf%qflx_rootsoi, & 
     !$acc col_wf%qflx_rootsoi, & 
     !$acc col_wf%qflx_deficit, & 
     !$acc col_wf%qflx_drain, & 
     !$acc col_wf%qflx_rsub_sat, & 
     !$acc col_wf%qflx_drain_perched, & 
     !$acc col_wf%qflx_sub_snow, & 
     !$acc col_wf%mflx_snowlyr, & 
     !$acc col_wf%qflx_snow2topsoi, & 
     !$acc col_wf%qflx_sl_top_soil )
     !$acc update self(& 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%flx_bc_dep_wet_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry1_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry2_col, & 
     !$acc aerosol_vars%flx_bc_dep_dry_col, & 
     !$acc aerosol_vars%flx_oc_dep_wet_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet1_col, & 
     !$acc aerosol_vars%flx_bc_dep_pho_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet4_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry4_col, & 
     !$acc aerosol_vars%flx_oc_dep_dry_col, & 
     !$acc aerosol_vars%flx_bc_dep_phi_col, & 
     !$acc aerosol_vars%flx_oc_dep_phi_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%flx_dst_dep_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%flx_bc_dep_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet2_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet3_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%flx_oc_dep_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%flx_oc_dep_pho_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry3_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst2_col )
     !$acc update self(& 
     !$acc soilhydrology_vars%moist_col, & 
     !$acc soilhydrology_vars%moist_vol_col, & 
     !$acc soilhydrology_vars%ice_col, & 
     !$acc soilhydrology_vars%max_infil_col, & 
     !$acc soilhydrology_vars%fracice_col, & 
     !$acc soilhydrology_vars%icefrac_col, & 
     !$acc soilhydrology_vars%i_0_col, & 
     !$acc soilhydrology_vars%fcov_col, & 
     !$acc soilhydrology_vars%fsat_col, & 
     !$acc soilhydrology_vars%icefrac_col, & 
     !$acc soilhydrology_vars%frost_table_col, & 
     !$acc soilhydrology_vars%icefrac_col, & 
     !$acc soilhydrology_vars%zwt_col, & 
     !$acc soilhydrology_vars%zwt_perched_col, & 
     !$acc soilhydrology_vars%wa_col, & 
     !$acc soilhydrology_vars%icefrac_col, & 
     !$acc soilhydrology_vars%qcharge_col, & 
     !$acc soilhydrology_vars%frost_table_col, & 
     !$acc soilhydrology_vars%qcharge_col, & 
     !$acc soilhydrology_vars%zwt_perched_col, & 
     !$acc soilhydrology_vars%wa_col, & 
     !$acc soilhydrology_vars%zwt_col )
     !$acc update self(& 
     !$acc veg_wf%qflx_rootsoi_frac )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_es%t_soi17cm' 
     write(10,*) col_es%t_soi17cm
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%dTdz_top' 
     write(10,*) col_es%dTdz_top
     write(10,*) 'col_es%snot_top' 
     write(10,*) col_es%snot_top
     write(10,*) 'col_es%t_grnd_r' 
     write(10,*) col_es%t_grnd_r
     write(10,*) 'col_es%t_soi10cm' 
     write(10,*) col_es%t_soi10cm
     write(10,*) 'col_es%t_grnd_u' 
     write(10,*) col_es%t_grnd_u
     write(10,*) 'col_es%t_grnd' 
     write(10,*) col_es%t_grnd
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%wf' 
     write(10,*) col_ws%wf
     write(10,*) 'col_ws%snow_persistence' 
     write(10,*) col_ws%snow_persistence
     write(10,*) 'col_ws%h2osoi_vol' 
     write(10,*) col_ws%h2osoi_vol
     write(10,*) 'col_ws%h2osno_top' 
     write(10,*) col_ws%h2osno_top
     write(10,*) 'col_ws%snw_rds_top' 
     write(10,*) col_ws%snw_rds_top
     write(10,*) 'col_ws%snowliq' 
     write(10,*) col_ws%snowliq
     write(10,*) 'col_ws%h2osoi_liqice_10cm' 
     write(10,*) col_ws%h2osoi_liqice_10cm
     write(10,*) 'col_ws%wf2' 
     write(10,*) col_ws%wf2
     write(10,*) 'col_ws%air_vol' 
     write(10,*) col_ws%air_vol
     write(10,*) 'col_ws%h2osoi_liqvol' 
     write(10,*) col_ws%h2osoi_liqvol
     write(10,*) 'col_ws%snowice' 
     write(10,*) col_ws%snowice
     write(10,*) 'col_ws%sno_liq_top' 
     write(10,*) col_ws%sno_liq_top
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%snowdp' 
     write(10,*) col_ws%snowdp
     write(10,*) 'col_ws%h2osoi_icevol' 
     write(10,*) col_ws%h2osoi_icevol
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%int_snow' 
     write(10,*) col_ws%int_snow
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_vol' 
     write(10,*) col_ws%h2osoi_vol
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osoi_vol' 
     write(10,*) col_ws%h2osoi_vol
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%int_snow' 
     write(10,*) col_ws%int_snow
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%frac_sno' 
     write(10,*) col_ws%frac_sno
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'soilstate_vars%soilpsi_col' 
     write(10,*) soilstate_vars%soilpsi_col
     write(10,*) 'soilstate_vars%eff_porosity_col' 
     write(10,*) soilstate_vars%eff_porosity_col
     write(10,*) 'soilstate_vars%eff_porosity_col' 
     write(10,*) soilstate_vars%eff_porosity_col
     write(10,*) 'soilstate_vars%rootr_col' 
     write(10,*) soilstate_vars%rootr_col
     write(10,*) 'soilstate_vars%rootr_col' 
     write(10,*) soilstate_vars%rootr_col
     write(10,*) 'soilstate_vars%hk_l_col' 
     write(10,*) soilstate_vars%hk_l_col
     write(10,*) 'soilstate_vars%smp_l_col' 
     write(10,*) soilstate_vars%smp_l_col
     write(10,*) 'col_wf%qflx_snow_melt' 
     write(10,*) col_wf%qflx_snow_melt
     write(10,*) 'col_wf%qflx_top_soil' 
     write(10,*) col_wf%qflx_top_soil
     write(10,*) 'col_wf%mflx_neg_snow_1d' 
     write(10,*) col_wf%mflx_neg_snow_1d
     write(10,*) 'col_wf%qflx_top_soil' 
     write(10,*) col_wf%qflx_top_soil
     write(10,*) 'col_wf%qflx_surf' 
     write(10,*) col_wf%qflx_surf
     write(10,*) 'col_wf%qflx_infl' 
     write(10,*) col_wf%qflx_infl
     write(10,*) 'col_wf%qflx_h2osfc_surf' 
     write(10,*) col_wf%qflx_h2osfc_surf
     write(10,*) 'col_wf%qflx_gross_evap_soil' 
     write(10,*) col_wf%qflx_gross_evap_soil
     write(10,*) 'col_wf%qflx_surf' 
     write(10,*) col_wf%qflx_surf
     write(10,*) 'col_wf%qflx_gross_infl_soil' 
     write(10,*) col_wf%qflx_gross_infl_soil
     write(10,*) 'col_wf%qflx_drain_perched' 
     write(10,*) col_wf%qflx_drain_perched
     write(10,*) 'col_wf%qflx_qrgwl' 
     write(10,*) col_wf%qflx_qrgwl
     write(10,*) 'col_wf%mflx_drain_perched_1d' 
     write(10,*) col_wf%mflx_drain_perched_1d
     write(10,*) 'col_wf%qflx_drain' 
     write(10,*) col_wf%qflx_drain
     write(10,*) 'col_wf%qflx_rsub_sat' 
     write(10,*) col_wf%qflx_rsub_sat
     write(10,*) 'col_wf%qflx_rootsoi' 
     write(10,*) col_wf%qflx_rootsoi
     write(10,*) 'col_wf%qflx_rootsoi' 
     write(10,*) col_wf%qflx_rootsoi
     write(10,*) 'col_wf%qflx_deficit' 
     write(10,*) col_wf%qflx_deficit
     write(10,*) 'col_wf%qflx_drain' 
     write(10,*) col_wf%qflx_drain
     write(10,*) 'col_wf%qflx_rsub_sat' 
     write(10,*) col_wf%qflx_rsub_sat
     write(10,*) 'col_wf%qflx_drain_perched' 
     write(10,*) col_wf%qflx_drain_perched
     write(10,*) 'col_wf%qflx_sub_snow' 
     write(10,*) col_wf%qflx_sub_snow
     write(10,*) 'col_wf%mflx_snowlyr' 
     write(10,*) col_wf%mflx_snowlyr
     write(10,*) 'col_wf%qflx_snow2topsoi' 
     write(10,*) col_wf%qflx_snow2topsoi
     write(10,*) 'col_wf%qflx_sl_top_soil' 
     write(10,*) col_wf%qflx_sl_top_soil
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%flx_bc_dep_wet_col' 
     write(10,*) aerosol_vars%flx_bc_dep_wet_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry1_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry1_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry2_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry2_col
     write(10,*) 'aerosol_vars%flx_bc_dep_dry_col' 
     write(10,*) aerosol_vars%flx_bc_dep_dry_col
     write(10,*) 'aerosol_vars%flx_oc_dep_wet_col' 
     write(10,*) aerosol_vars%flx_oc_dep_wet_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet1_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet1_col
     write(10,*) 'aerosol_vars%flx_bc_dep_pho_col' 
     write(10,*) aerosol_vars%flx_bc_dep_pho_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet4_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet4_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry4_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry4_col
     write(10,*) 'aerosol_vars%flx_oc_dep_dry_col' 
     write(10,*) aerosol_vars%flx_oc_dep_dry_col
     write(10,*) 'aerosol_vars%flx_bc_dep_phi_col' 
     write(10,*) aerosol_vars%flx_bc_dep_phi_col
     write(10,*) 'aerosol_vars%flx_oc_dep_phi_col' 
     write(10,*) aerosol_vars%flx_oc_dep_phi_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%flx_dst_dep_col' 
     write(10,*) aerosol_vars%flx_dst_dep_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%flx_bc_dep_col' 
     write(10,*) aerosol_vars%flx_bc_dep_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet2_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet2_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet3_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet3_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%flx_oc_dep_col' 
     write(10,*) aerosol_vars%flx_oc_dep_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%flx_oc_dep_pho_col' 
     write(10,*) aerosol_vars%flx_oc_dep_pho_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry3_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'soilhydrology_vars%moist_col' 
     write(10,*) soilhydrology_vars%moist_col
     write(10,*) 'soilhydrology_vars%moist_vol_col' 
     write(10,*) soilhydrology_vars%moist_vol_col
     write(10,*) 'soilhydrology_vars%ice_col' 
     write(10,*) soilhydrology_vars%ice_col
     write(10,*) 'soilhydrology_vars%max_infil_col' 
     write(10,*) soilhydrology_vars%max_infil_col
     write(10,*) 'soilhydrology_vars%fracice_col' 
     write(10,*) soilhydrology_vars%fracice_col
     write(10,*) 'soilhydrology_vars%icefrac_col' 
     write(10,*) soilhydrology_vars%icefrac_col
     write(10,*) 'soilhydrology_vars%i_0_col' 
     write(10,*) soilhydrology_vars%i_0_col
     write(10,*) 'soilhydrology_vars%fcov_col' 
     write(10,*) soilhydrology_vars%fcov_col
     write(10,*) 'soilhydrology_vars%fsat_col' 
     write(10,*) soilhydrology_vars%fsat_col
     write(10,*) 'soilhydrology_vars%icefrac_col' 
     write(10,*) soilhydrology_vars%icefrac_col
     write(10,*) 'soilhydrology_vars%frost_table_col' 
     write(10,*) soilhydrology_vars%frost_table_col
     write(10,*) 'soilhydrology_vars%icefrac_col' 
     write(10,*) soilhydrology_vars%icefrac_col
     write(10,*) 'soilhydrology_vars%zwt_col' 
     write(10,*) soilhydrology_vars%zwt_col
     write(10,*) 'soilhydrology_vars%zwt_perched_col' 
     write(10,*) soilhydrology_vars%zwt_perched_col
     write(10,*) 'soilhydrology_vars%wa_col' 
     write(10,*) soilhydrology_vars%wa_col
     write(10,*) 'soilhydrology_vars%icefrac_col' 
     write(10,*) soilhydrology_vars%icefrac_col
     write(10,*) 'soilhydrology_vars%qcharge_col' 
     write(10,*) soilhydrology_vars%qcharge_col
     write(10,*) 'soilhydrology_vars%frost_table_col' 
     write(10,*) soilhydrology_vars%frost_table_col
     write(10,*) 'soilhydrology_vars%qcharge_col' 
     write(10,*) soilhydrology_vars%qcharge_col
     write(10,*) 'soilhydrology_vars%zwt_perched_col' 
     write(10,*) soilhydrology_vars%zwt_perched_col
     write(10,*) 'soilhydrology_vars%wa_col' 
     write(10,*) soilhydrology_vars%wa_col
     write(10,*) 'soilhydrology_vars%zwt_col' 
     write(10,*) soilhydrology_vars%zwt_col
     write(10,*) 'veg_wf%qflx_rootsoi_frac' 
     write(10,*) veg_wf%qflx_rootsoi_frac
     close(10)
end subroutine 
subroutine update_vars_AerosolMasses(gpu)
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : aerosol_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_AerosolMasses.txt"
     else
          file='cpu_AerosolMasses.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%h2osno_top )
     !$acc update self(& 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_cnc_bcpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst_col_col, & 
     !$acc aerosol_vars%mss_oc_top_col, & 
     !$acc aerosol_vars%mss_cnc_dst1_col, & 
     !$acc aerosol_vars%mss_oc_col_col, & 
     !$acc aerosol_vars%mss_dsttot_col, & 
     !$acc aerosol_vars%mss_cnc_dst4_col, & 
     !$acc aerosol_vars%mss_bc_col_col, & 
     !$acc aerosol_vars%mss_cnc_dst3_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_bctot_col, & 
     !$acc aerosol_vars%mss_cnc_ocpho_col, & 
     !$acc aerosol_vars%mss_cnc_ocphi_col, & 
     !$acc aerosol_vars%mss_cnc_dst2_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bc_top_col, & 
     !$acc aerosol_vars%mss_octot_col, & 
     !$acc aerosol_vars%mss_dst_top_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_cnc_bcphi_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%h2osno_top' 
     write(10,*) col_ws%h2osno_top
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_cnc_bcpho_col' 
     write(10,*) aerosol_vars%mss_cnc_bcpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst_col_col' 
     write(10,*) aerosol_vars%mss_dst_col_col
     write(10,*) 'aerosol_vars%mss_oc_top_col' 
     write(10,*) aerosol_vars%mss_oc_top_col
     write(10,*) 'aerosol_vars%mss_cnc_dst1_col' 
     write(10,*) aerosol_vars%mss_cnc_dst1_col
     write(10,*) 'aerosol_vars%mss_oc_col_col' 
     write(10,*) aerosol_vars%mss_oc_col_col
     write(10,*) 'aerosol_vars%mss_dsttot_col' 
     write(10,*) aerosol_vars%mss_dsttot_col
     write(10,*) 'aerosol_vars%mss_cnc_dst4_col' 
     write(10,*) aerosol_vars%mss_cnc_dst4_col
     write(10,*) 'aerosol_vars%mss_bc_col_col' 
     write(10,*) aerosol_vars%mss_bc_col_col
     write(10,*) 'aerosol_vars%mss_cnc_dst3_col' 
     write(10,*) aerosol_vars%mss_cnc_dst3_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_bctot_col' 
     write(10,*) aerosol_vars%mss_bctot_col
     write(10,*) 'aerosol_vars%mss_cnc_ocpho_col' 
     write(10,*) aerosol_vars%mss_cnc_ocpho_col
     write(10,*) 'aerosol_vars%mss_cnc_ocphi_col' 
     write(10,*) aerosol_vars%mss_cnc_ocphi_col
     write(10,*) 'aerosol_vars%mss_cnc_dst2_col' 
     write(10,*) aerosol_vars%mss_cnc_dst2_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bc_top_col' 
     write(10,*) aerosol_vars%mss_bc_top_col
     write(10,*) 'aerosol_vars%mss_octot_col' 
     write(10,*) aerosol_vars%mss_octot_col
     write(10,*) 'aerosol_vars%mss_dst_top_col' 
     write(10,*) aerosol_vars%mss_dst_top_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_cnc_bcphi_col' 
     write(10,*) aerosol_vars%mss_cnc_bcphi_col
     close(10)
end subroutine 
subroutine update_vars_LakeHydrology(gpu)
     use ColumnType, only : col_pp 
     use ColumnDataType, only : col_es 
     use ColumnDataType, only : col_ws 
     use ColumnDataType, only : col_wf 
     use VegetationDataType, only : veg_wf 
     use ColumnDataType, only : col_ef 
     use VegetationDataType, only : veg_ef 
     use clm_instMod, only : lakestate_vars 
     use clm_instMod, only : aerosol_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_LakeHydrology.txt"
     else
          file='cpu_LakeHydrology.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z, & 
     !$acc col_pp%dz, & 
     !$acc col_pp%snl, & 
     !$acc col_pp%zi, & 
     !$acc col_pp%z )
     !$acc update self(& 
     !$acc col_es%t_soisno, & 
     !$acc col_es%snot_top, & 
     !$acc col_es%t_lake, & 
     !$acc col_es%dTdz_top, & 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osno_top, & 
     !$acc col_ws%frac_iceold, & 
     !$acc col_ws%h2osoi_vol, & 
     !$acc col_ws%snw_rds_top, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%snowliq, & 
     !$acc col_ws%begwb, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%snowice, & 
     !$acc col_ws%sno_liq_top, & 
     !$acc col_ws%endwb, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%snow_depth, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%int_snow, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%snow_depth, & 
     !$acc col_ws%frac_sno_eff, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%int_snow, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%frac_sno, & 
     !$acc col_ws%h2osno, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%h2osoi_ice )
     !$acc update self(& 
     !$acc col_wf%qflx_irrig, & 
     !$acc col_wf%qflx_h2osfc_surf, & 
     !$acc col_wf%qflx_drain_perched, & 
     !$acc col_wf%qflx_sub_snow, & 
     !$acc col_wf%qflx_rsub_sat, & 
     !$acc col_wf%qflx_infl, & 
     !$acc col_wf%qflx_floodc, & 
     !$acc col_wf%qflx_runoff, & 
     !$acc col_wf%qflx_snomelt, & 
     !$acc col_wf%qflx_snow_melt, & 
     !$acc col_wf%qflx_dew_grnd, & 
     !$acc col_wf%qflx_evap_tot, & 
     !$acc col_wf%qflx_drain, & 
     !$acc col_wf%qflx_rain_grnd, & 
     !$acc col_wf%qflx_surf, & 
     !$acc col_wf%qflx_leafdrip, & 
     !$acc col_wf%qflx_snow_grnd, & 
     !$acc col_wf%qflx_prec_grnd, & 
     !$acc col_wf%qflx_snwcp_ice, & 
     !$acc col_wf%qflx_top_soil, & 
     !$acc col_wf%qflx_snwcp_liq, & 
     !$acc col_wf%qflx_evap_grnd, & 
     !$acc col_wf%qflx_qrgwl, & 
     !$acc col_wf%qflx_dirct_rain, & 
     !$acc col_wf%qflx_dew_snow, & 
     !$acc col_wf%qflx_sl_top_soil, & 
     !$acc col_wf%qflx_snow_melt, & 
     !$acc col_wf%qflx_top_soil, & 
     !$acc col_wf%mflx_neg_snow_1d, & 
     !$acc col_wf%mflx_snowlyr, & 
     !$acc col_wf%qflx_snow2topsoi, & 
     !$acc col_wf%qflx_sl_top_soil )
     !$acc update self(& 
     !$acc veg_wf%qflx_rain_grnd, & 
     !$acc veg_wf%qflx_irrig_patch, & 
     !$acc veg_wf%qflx_snwcp_ice, & 
     !$acc veg_wf%qflx_sub_snow, & 
     !$acc veg_wf%qflx_snwcp_liq, & 
     !$acc veg_wf%qflx_evap_grnd, & 
     !$acc veg_wf%qflx_dew_grnd, & 
     !$acc veg_wf%qflx_dirct_rain, & 
     !$acc veg_wf%qflx_dew_snow, & 
     !$acc veg_wf%qflx_leafdrip, & 
     !$acc veg_wf%qflx_snow_grnd, & 
     !$acc veg_wf%qflx_prec_grnd )
     !$acc update self(& 
     !$acc col_ef%eflx_snomelt )
     !$acc update self(& 
     !$acc veg_ef%eflx_sh_grnd, & 
     !$acc veg_ef%eflx_gnet, & 
     !$acc veg_ef%eflx_grnd_lake, & 
     !$acc veg_ef%eflx_sh_tot, & 
     !$acc veg_ef%eflx_soil_grnd )
     !$acc update self(& 
     !$acc lakestate_vars%lake_icefrac_col )
     !$acc update self(& 
     !$acc aerosol_vars%mss_dsttot_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_octot_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_dst_top_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst_col_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_oc_top_col, & 
     !$acc aerosol_vars%mss_bctot_col, & 
     !$acc aerosol_vars%mss_bc_col_col, & 
     !$acc aerosol_vars%mss_oc_col_col, & 
     !$acc aerosol_vars%mss_bc_top_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%flx_bc_dep_wet_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry1_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry2_col, & 
     !$acc aerosol_vars%flx_bc_dep_dry_col, & 
     !$acc aerosol_vars%flx_oc_dep_wet_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet1_col, & 
     !$acc aerosol_vars%flx_bc_dep_pho_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet4_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry4_col, & 
     !$acc aerosol_vars%flx_oc_dep_dry_col, & 
     !$acc aerosol_vars%flx_bc_dep_phi_col, & 
     !$acc aerosol_vars%flx_oc_dep_phi_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%flx_dst_dep_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%flx_bc_dep_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet2_col, & 
     !$acc aerosol_vars%flx_dst_dep_wet3_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%flx_oc_dep_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%flx_oc_dep_pho_col, & 
     !$acc aerosol_vars%flx_dst_dep_dry3_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst2_col, & 
     !$acc aerosol_vars%mss_dst1_col, & 
     !$acc aerosol_vars%mss_bcpho_col, & 
     !$acc aerosol_vars%mss_bcphi_col, & 
     !$acc aerosol_vars%mss_dst4_col, & 
     !$acc aerosol_vars%mss_ocphi_col, & 
     !$acc aerosol_vars%mss_ocpho_col, & 
     !$acc aerosol_vars%mss_dst3_col, & 
     !$acc aerosol_vars%mss_dst2_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_pp%dz' 
     write(10,*) col_pp%dz
     write(10,*) 'col_pp%snl' 
     write(10,*) col_pp%snl
     write(10,*) 'col_pp%zi' 
     write(10,*) col_pp%zi
     write(10,*) 'col_pp%z' 
     write(10,*) col_pp%z
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%snot_top' 
     write(10,*) col_es%snot_top
     write(10,*) 'col_es%t_lake' 
     write(10,*) col_es%t_lake
     write(10,*) 'col_es%dTdz_top' 
     write(10,*) col_es%dTdz_top
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osno_top' 
     write(10,*) col_ws%h2osno_top
     write(10,*) 'col_ws%frac_iceold' 
     write(10,*) col_ws%frac_iceold
     write(10,*) 'col_ws%h2osoi_vol' 
     write(10,*) col_ws%h2osoi_vol
     write(10,*) 'col_ws%snw_rds_top' 
     write(10,*) col_ws%snw_rds_top
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%snowliq' 
     write(10,*) col_ws%snowliq
     write(10,*) 'col_ws%begwb' 
     write(10,*) col_ws%begwb
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%snowice' 
     write(10,*) col_ws%snowice
     write(10,*) 'col_ws%sno_liq_top' 
     write(10,*) col_ws%sno_liq_top
     write(10,*) 'col_ws%endwb' 
     write(10,*) col_ws%endwb
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%int_snow' 
     write(10,*) col_ws%int_snow
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%snow_depth' 
     write(10,*) col_ws%snow_depth
     write(10,*) 'col_ws%frac_sno_eff' 
     write(10,*) col_ws%frac_sno_eff
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%int_snow' 
     write(10,*) col_ws%int_snow
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%frac_sno' 
     write(10,*) col_ws%frac_sno
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_wf%qflx_irrig' 
     write(10,*) col_wf%qflx_irrig
     write(10,*) 'col_wf%qflx_h2osfc_surf' 
     write(10,*) col_wf%qflx_h2osfc_surf
     write(10,*) 'col_wf%qflx_drain_perched' 
     write(10,*) col_wf%qflx_drain_perched
     write(10,*) 'col_wf%qflx_sub_snow' 
     write(10,*) col_wf%qflx_sub_snow
     write(10,*) 'col_wf%qflx_rsub_sat' 
     write(10,*) col_wf%qflx_rsub_sat
     write(10,*) 'col_wf%qflx_infl' 
     write(10,*) col_wf%qflx_infl
     write(10,*) 'col_wf%qflx_floodc' 
     write(10,*) col_wf%qflx_floodc
     write(10,*) 'col_wf%qflx_runoff' 
     write(10,*) col_wf%qflx_runoff
     write(10,*) 'col_wf%qflx_snomelt' 
     write(10,*) col_wf%qflx_snomelt
     write(10,*) 'col_wf%qflx_snow_melt' 
     write(10,*) col_wf%qflx_snow_melt
     write(10,*) 'col_wf%qflx_dew_grnd' 
     write(10,*) col_wf%qflx_dew_grnd
     write(10,*) 'col_wf%qflx_evap_tot' 
     write(10,*) col_wf%qflx_evap_tot
     write(10,*) 'col_wf%qflx_drain' 
     write(10,*) col_wf%qflx_drain
     write(10,*) 'col_wf%qflx_rain_grnd' 
     write(10,*) col_wf%qflx_rain_grnd
     write(10,*) 'col_wf%qflx_surf' 
     write(10,*) col_wf%qflx_surf
     write(10,*) 'col_wf%qflx_leafdrip' 
     write(10,*) col_wf%qflx_leafdrip
     write(10,*) 'col_wf%qflx_snow_grnd' 
     write(10,*) col_wf%qflx_snow_grnd
     write(10,*) 'col_wf%qflx_prec_grnd' 
     write(10,*) col_wf%qflx_prec_grnd
     write(10,*) 'col_wf%qflx_snwcp_ice' 
     write(10,*) col_wf%qflx_snwcp_ice
     write(10,*) 'col_wf%qflx_top_soil' 
     write(10,*) col_wf%qflx_top_soil
     write(10,*) 'col_wf%qflx_snwcp_liq' 
     write(10,*) col_wf%qflx_snwcp_liq
     write(10,*) 'col_wf%qflx_evap_grnd' 
     write(10,*) col_wf%qflx_evap_grnd
     write(10,*) 'col_wf%qflx_qrgwl' 
     write(10,*) col_wf%qflx_qrgwl
     write(10,*) 'col_wf%qflx_dirct_rain' 
     write(10,*) col_wf%qflx_dirct_rain
     write(10,*) 'col_wf%qflx_dew_snow' 
     write(10,*) col_wf%qflx_dew_snow
     write(10,*) 'col_wf%qflx_sl_top_soil' 
     write(10,*) col_wf%qflx_sl_top_soil
     write(10,*) 'col_wf%qflx_snow_melt' 
     write(10,*) col_wf%qflx_snow_melt
     write(10,*) 'col_wf%qflx_top_soil' 
     write(10,*) col_wf%qflx_top_soil
     write(10,*) 'col_wf%mflx_neg_snow_1d' 
     write(10,*) col_wf%mflx_neg_snow_1d
     write(10,*) 'col_wf%mflx_snowlyr' 
     write(10,*) col_wf%mflx_snowlyr
     write(10,*) 'col_wf%qflx_snow2topsoi' 
     write(10,*) col_wf%qflx_snow2topsoi
     write(10,*) 'col_wf%qflx_sl_top_soil' 
     write(10,*) col_wf%qflx_sl_top_soil
     write(10,*) 'veg_wf%qflx_rain_grnd' 
     write(10,*) veg_wf%qflx_rain_grnd
     write(10,*) 'veg_wf%qflx_irrig_patch' 
     write(10,*) veg_wf%qflx_irrig_patch
     write(10,*) 'veg_wf%qflx_snwcp_ice' 
     write(10,*) veg_wf%qflx_snwcp_ice
     write(10,*) 'veg_wf%qflx_sub_snow' 
     write(10,*) veg_wf%qflx_sub_snow
     write(10,*) 'veg_wf%qflx_snwcp_liq' 
     write(10,*) veg_wf%qflx_snwcp_liq
     write(10,*) 'veg_wf%qflx_evap_grnd' 
     write(10,*) veg_wf%qflx_evap_grnd
     write(10,*) 'veg_wf%qflx_dew_grnd' 
     write(10,*) veg_wf%qflx_dew_grnd
     write(10,*) 'veg_wf%qflx_dirct_rain' 
     write(10,*) veg_wf%qflx_dirct_rain
     write(10,*) 'veg_wf%qflx_dew_snow' 
     write(10,*) veg_wf%qflx_dew_snow
     write(10,*) 'veg_wf%qflx_leafdrip' 
     write(10,*) veg_wf%qflx_leafdrip
     write(10,*) 'veg_wf%qflx_snow_grnd' 
     write(10,*) veg_wf%qflx_snow_grnd
     write(10,*) 'veg_wf%qflx_prec_grnd' 
     write(10,*) veg_wf%qflx_prec_grnd
     write(10,*) 'col_ef%eflx_snomelt' 
     write(10,*) col_ef%eflx_snomelt
     write(10,*) 'veg_ef%eflx_sh_grnd' 
     write(10,*) veg_ef%eflx_sh_grnd
     write(10,*) 'veg_ef%eflx_gnet' 
     write(10,*) veg_ef%eflx_gnet
     write(10,*) 'veg_ef%eflx_grnd_lake' 
     write(10,*) veg_ef%eflx_grnd_lake
     write(10,*) 'veg_ef%eflx_sh_tot' 
     write(10,*) veg_ef%eflx_sh_tot
     write(10,*) 'veg_ef%eflx_soil_grnd' 
     write(10,*) veg_ef%eflx_soil_grnd
     write(10,*) 'lakestate_vars%lake_icefrac_col' 
     write(10,*) lakestate_vars%lake_icefrac_col
     write(10,*) 'aerosol_vars%mss_dsttot_col' 
     write(10,*) aerosol_vars%mss_dsttot_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_octot_col' 
     write(10,*) aerosol_vars%mss_octot_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_dst_top_col' 
     write(10,*) aerosol_vars%mss_dst_top_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst_col_col' 
     write(10,*) aerosol_vars%mss_dst_col_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_oc_top_col' 
     write(10,*) aerosol_vars%mss_oc_top_col
     write(10,*) 'aerosol_vars%mss_bctot_col' 
     write(10,*) aerosol_vars%mss_bctot_col
     write(10,*) 'aerosol_vars%mss_bc_col_col' 
     write(10,*) aerosol_vars%mss_bc_col_col
     write(10,*) 'aerosol_vars%mss_oc_col_col' 
     write(10,*) aerosol_vars%mss_oc_col_col
     write(10,*) 'aerosol_vars%mss_bc_top_col' 
     write(10,*) aerosol_vars%mss_bc_top_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%flx_bc_dep_wet_col' 
     write(10,*) aerosol_vars%flx_bc_dep_wet_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry1_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry1_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry2_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry2_col
     write(10,*) 'aerosol_vars%flx_bc_dep_dry_col' 
     write(10,*) aerosol_vars%flx_bc_dep_dry_col
     write(10,*) 'aerosol_vars%flx_oc_dep_wet_col' 
     write(10,*) aerosol_vars%flx_oc_dep_wet_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet1_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet1_col
     write(10,*) 'aerosol_vars%flx_bc_dep_pho_col' 
     write(10,*) aerosol_vars%flx_bc_dep_pho_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet4_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet4_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry4_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry4_col
     write(10,*) 'aerosol_vars%flx_oc_dep_dry_col' 
     write(10,*) aerosol_vars%flx_oc_dep_dry_col
     write(10,*) 'aerosol_vars%flx_bc_dep_phi_col' 
     write(10,*) aerosol_vars%flx_bc_dep_phi_col
     write(10,*) 'aerosol_vars%flx_oc_dep_phi_col' 
     write(10,*) aerosol_vars%flx_oc_dep_phi_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%flx_dst_dep_col' 
     write(10,*) aerosol_vars%flx_dst_dep_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%flx_bc_dep_col' 
     write(10,*) aerosol_vars%flx_bc_dep_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet2_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet2_col
     write(10,*) 'aerosol_vars%flx_dst_dep_wet3_col' 
     write(10,*) aerosol_vars%flx_dst_dep_wet3_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%flx_oc_dep_col' 
     write(10,*) aerosol_vars%flx_oc_dep_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%flx_oc_dep_pho_col' 
     write(10,*) aerosol_vars%flx_oc_dep_pho_col
     write(10,*) 'aerosol_vars%flx_dst_dep_dry3_col' 
     write(10,*) aerosol_vars%flx_dst_dep_dry3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     write(10,*) 'aerosol_vars%mss_dst1_col' 
     write(10,*) aerosol_vars%mss_dst1_col
     write(10,*) 'aerosol_vars%mss_bcpho_col' 
     write(10,*) aerosol_vars%mss_bcpho_col
     write(10,*) 'aerosol_vars%mss_bcphi_col' 
     write(10,*) aerosol_vars%mss_bcphi_col
     write(10,*) 'aerosol_vars%mss_dst4_col' 
     write(10,*) aerosol_vars%mss_dst4_col
     write(10,*) 'aerosol_vars%mss_ocphi_col' 
     write(10,*) aerosol_vars%mss_ocphi_col
     write(10,*) 'aerosol_vars%mss_ocpho_col' 
     write(10,*) aerosol_vars%mss_ocpho_col
     write(10,*) 'aerosol_vars%mss_dst3_col' 
     write(10,*) aerosol_vars%mss_dst3_col
     write(10,*) 'aerosol_vars%mss_dst2_col' 
     write(10,*) aerosol_vars%mss_dst2_col
     close(10)
end subroutine 
subroutine update_vars_SnowAge_grain(gpu)
     use ColumnDataType, only : col_ws 
     use ColumnDataType, only : col_es 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SnowAge_grain.txt"
     else
          file='cpu_SnowAge_grain.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%sno_liq_top, & 
     !$acc col_ws%snw_rds, & 
     !$acc col_ws%snw_rds_top )
     !$acc update self(& 
     !$acc col_es%snot_top, & 
     !$acc col_es%dTdz_top )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%sno_liq_top' 
     write(10,*) col_ws%sno_liq_top
     write(10,*) 'col_ws%snw_rds' 
     write(10,*) col_ws%snw_rds
     write(10,*) 'col_ws%snw_rds_top' 
     write(10,*) col_ws%snw_rds_top
     write(10,*) 'col_es%snot_top' 
     write(10,*) col_es%snot_top
     write(10,*) 'col_es%dTdz_top' 
     write(10,*) col_es%dTdz_top
     close(10)
end subroutine 
subroutine update_vars_EcosystemDynNoLeaching1(gpu)
     use ColumnDataType, only : col_nf 
     use VegetationDataType, only : veg_nf 
     use VegetationDataType, only : veg_cf 
     use ColumnDataType, only : col_pf 
     use VegetationDataType, only : veg_pf 
     use ColumnDataType, only : col_ps 
     use ColumnDataType, only : col_cf 
     use CNDecompCascadeConType, only : decomp_cascade_con 
     use clm_instMod, only : cnstate_vars 
     use VegetationDataType, only : veg_ns 
     use VegetationDataType, only : c14_veg_cf 
     use VegetationDataType, only : c13_veg_cf 
     use clm_varctl        , only : use_c13, use_c14
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_EcosystemDynNoLeaching1.txt"
     else
          file='cpu_EcosystemDynNoLeaching1.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_nf%ndep_to_sminn, & 
     !$acc col_nf%nfix_to_sminn, & 
     !$acc col_nf%nfix_to_sminn, & 
     !$acc col_nf%nfix_to_ecosysn, & 
     !$acc col_nf%fert_to_sminn, & 
     !$acc col_nf%soyfixn_to_sminn, & 
     !$acc col_nf%plant_ndemand, & 
     !$acc col_nf%plant_ndemand_vr )
     !$acc update self(& 
     !$acc veg_nf%nfix_to_plantn, & 
     !$acc veg_nf%fert, & 
     !$acc veg_nf%soyfixn, & 
     !$acc veg_nf%retransn_to_npool, & 
     !$acc veg_nf%livestemn_to_retransn, & 
     !$acc veg_nf%plant_ndemand, & 
     !$acc veg_nf%frootn_to_retransn, & 
     !$acc veg_nf%leafn_to_retransn, & 
     !$acc veg_nf%avail_retransn )
     !$acc update self(& 
     !$acc veg_cf%xr, & 
     !$acc veg_cf%livecroot_mr, & 
     !$acc veg_cf%froot_mr, & 
     !$acc veg_cf%livestem_mr, & 
     !$acc veg_cf%grain_mr, & 
     !$acc veg_cf%leaf_mr, & 
     !$acc veg_cf%froot_curmr, & 
     !$acc veg_cf%xsmrpool_recover, & 
     !$acc veg_cf%availc, & 
     !$acc veg_cf%grain_xsmr, & 
     !$acc veg_cf%froot_xsmr, & 
     !$acc veg_cf%livecroot_curmr, & 
     !$acc veg_cf%grain_curmr, & 
     !$acc veg_cf%psnshade_to_cpool, & 
     !$acc veg_cf%gpp_before_downreg, & 
     !$acc veg_cf%psnsun_to_cpool, & 
     !$acc veg_cf%cpool_to_xsmrpool, & 
     !$acc veg_cf%leaf_curmr, & 
     !$acc veg_cf%leaf_xsmr, & 
     !$acc veg_cf%livestem_curmr, & 
     !$acc veg_cf%livecroot_xsmr, & 
     !$acc veg_cf%livestem_xsmr )
     !$acc update self(& 
     !$acc col_pf%primp_to_labilep_vr, & 
     !$acc col_pf%biochem_pmin_vr, & 
     !$acc col_pf%biochem_pmin_ppools_vr, & 
     !$acc col_pf%biochem_pmin_to_ecosysp_vr, & 
     !$acc col_pf%biochem_pmin_ppools_vr, & 
     !$acc col_pf%biochem_pmin_vr, & 
     !$acc col_pf%pdep_to_sminp, & 
     !$acc col_pf%plant_pdemand_vr, & 
     !$acc col_pf%plant_pdemand )
     !$acc update self(& 
     !$acc veg_pf%biochem_pmin_to_plant, & 
     !$acc veg_pf%retransp_to_ppool, & 
     !$acc veg_pf%avail_retransp, & 
     !$acc veg_pf%plant_pdemand )
     !$acc update self(& 
     !$acc col_ps%decomp_ppools_vr )
     !$acc update self(& 
     !$acc col_cf%o_scalar, & 
     !$acc col_cf%w_scalar, & 
     !$acc col_cf%decomp_k, & 
     !$acc col_cf%t_scalar, & 
     !$acc col_cf%o_scalar, & 
     !$acc col_cf%w_scalar, & 
     !$acc col_cf%decomp_k, & 
     !$acc col_cf%t_scalar )
     !$acc update self(& 
     !$acc decomp_cascade_con%decomp_k_pools, & 
     !$acc decomp_cascade_con%decomp_k_pools )
     !$acc update self(& 
     !$acc cnstate_vars%scalaravg_col, & 
     !$acc cnstate_vars%scalaravg_col, & 
     !$acc cnstate_vars%stem_prof_patch, & 
     !$acc cnstate_vars%croot_prof_patch, & 
     !$acc cnstate_vars%pdep_prof_col, & 
     !$acc cnstate_vars%nfixation_prof_col, & 
     !$acc cnstate_vars%ndep_prof_col, & 
     !$acc cnstate_vars%leaf_prof_patch, & 
     !$acc cnstate_vars%froot_prof_patch, & 
     !$acc cnstate_vars%grain_flag_patch, & 
     !$acc cnstate_vars%astem_patch, & 
     !$acc cnstate_vars%p_allometry_patch, & 
     !$acc cnstate_vars%tempmax_retransp_patch, & 
     !$acc cnstate_vars%tempmax_retransn_patch, & 
     !$acc cnstate_vars%c_allometry_patch, & 
     !$acc cnstate_vars%n_allometry_patch, & 
     !$acc cnstate_vars%tempsum_potential_gpp_patch, & 
     !$acc cnstate_vars%aleaf_patch, & 
     !$acc cnstate_vars%aleafi_patch, & 
     !$acc cnstate_vars%astemi_patch )
     !$acc update self(& 
     !$acc veg_ns%benefit_pgpp_pleafc )
     !$acc update self(& 
     !$acc c14_veg_cf%psnshade_to_cpool, & 
     !$acc c14_veg_cf%psnsun_to_cpool )
     !$acc update self(& 
     !$acc c13_veg_cf%psnshade_to_cpool, & 
     !$acc c13_veg_cf%psnsun_to_cpool )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_nf%ndep_to_sminn' 
     write(10,*) col_nf%ndep_to_sminn
     write(10,*) 'col_nf%nfix_to_sminn' 
     write(10,*) col_nf%nfix_to_sminn
     write(10,*) 'col_nf%nfix_to_sminn' 
     write(10,*) col_nf%nfix_to_sminn
     write(10,*) 'col_nf%nfix_to_ecosysn' 
     write(10,*) col_nf%nfix_to_ecosysn
     write(10,*) 'col_nf%fert_to_sminn' 
     write(10,*) col_nf%fert_to_sminn
     write(10,*) 'col_nf%soyfixn_to_sminn' 
     write(10,*) col_nf%soyfixn_to_sminn
     write(10,*) 'col_nf%plant_ndemand' 
     write(10,*) col_nf%plant_ndemand
     write(10,*) 'col_nf%plant_ndemand_vr' 
     write(10,*) col_nf%plant_ndemand_vr
     write(10,*) 'veg_nf%nfix_to_plantn' 
     write(10,*) veg_nf%nfix_to_plantn
     write(10,*) 'veg_nf%fert' 
     write(10,*) veg_nf%fert
     write(10,*) 'veg_nf%soyfixn' 
     write(10,*) veg_nf%soyfixn
     write(10,*) 'veg_nf%retransn_to_npool' 
     write(10,*) veg_nf%retransn_to_npool
     write(10,*) 'veg_nf%livestemn_to_retransn' 
     write(10,*) veg_nf%livestemn_to_retransn
     write(10,*) 'veg_nf%plant_ndemand' 
     write(10,*) veg_nf%plant_ndemand
     write(10,*) 'veg_nf%frootn_to_retransn' 
     write(10,*) veg_nf%frootn_to_retransn
     write(10,*) 'veg_nf%leafn_to_retransn' 
     write(10,*) veg_nf%leafn_to_retransn
     write(10,*) 'veg_nf%avail_retransn' 
     write(10,*) veg_nf%avail_retransn
     write(10,*) 'veg_cf%xr' 
     write(10,*) veg_cf%xr
     write(10,*) 'veg_cf%livecroot_mr' 
     write(10,*) veg_cf%livecroot_mr
     write(10,*) 'veg_cf%froot_mr' 
     write(10,*) veg_cf%froot_mr
     write(10,*) 'veg_cf%livestem_mr' 
     write(10,*) veg_cf%livestem_mr
     write(10,*) 'veg_cf%grain_mr' 
     write(10,*) veg_cf%grain_mr
     write(10,*) 'veg_cf%leaf_mr' 
     write(10,*) veg_cf%leaf_mr
     write(10,*) 'veg_cf%froot_curmr' 
     write(10,*) veg_cf%froot_curmr
     write(10,*) 'veg_cf%xsmrpool_recover' 
     write(10,*) veg_cf%xsmrpool_recover
     write(10,*) 'veg_cf%availc' 
     write(10,*) veg_cf%availc
     write(10,*) 'veg_cf%grain_xsmr' 
     write(10,*) veg_cf%grain_xsmr
     write(10,*) 'veg_cf%froot_xsmr' 
     write(10,*) veg_cf%froot_xsmr
     write(10,*) 'veg_cf%livecroot_curmr' 
     write(10,*) veg_cf%livecroot_curmr
     write(10,*) 'veg_cf%grain_curmr' 
     write(10,*) veg_cf%grain_curmr
     write(10,*) 'veg_cf%psnshade_to_cpool' 
     write(10,*) veg_cf%psnshade_to_cpool
     write(10,*) 'veg_cf%gpp_before_downreg' 
     write(10,*) veg_cf%gpp_before_downreg
     write(10,*) 'veg_cf%psnsun_to_cpool' 
     write(10,*) veg_cf%psnsun_to_cpool
     write(10,*) 'veg_cf%cpool_to_xsmrpool' 
     write(10,*) veg_cf%cpool_to_xsmrpool
     write(10,*) 'veg_cf%leaf_curmr' 
     write(10,*) veg_cf%leaf_curmr
     write(10,*) 'veg_cf%leaf_xsmr' 
     write(10,*) veg_cf%leaf_xsmr
     write(10,*) 'veg_cf%livestem_curmr' 
     write(10,*) veg_cf%livestem_curmr
     write(10,*) 'veg_cf%livecroot_xsmr' 
     write(10,*) veg_cf%livecroot_xsmr
     write(10,*) 'veg_cf%livestem_xsmr' 
     write(10,*) veg_cf%livestem_xsmr
     write(10,*) 'col_pf%primp_to_labilep_vr' 
     write(10,*) col_pf%primp_to_labilep_vr
     write(10,*) 'col_pf%biochem_pmin_vr' 
     write(10,*) col_pf%biochem_pmin_vr
     write(10,*) 'col_pf%biochem_pmin_ppools_vr' 
     write(10,*) col_pf%biochem_pmin_ppools_vr
     write(10,*) 'col_pf%biochem_pmin_to_ecosysp_vr' 
     write(10,*) col_pf%biochem_pmin_to_ecosysp_vr
     write(10,*) 'col_pf%biochem_pmin_ppools_vr' 
     write(10,*) col_pf%biochem_pmin_ppools_vr
     write(10,*) 'col_pf%biochem_pmin_vr' 
     write(10,*) col_pf%biochem_pmin_vr
     write(10,*) 'col_pf%pdep_to_sminp' 
     write(10,*) col_pf%pdep_to_sminp
     write(10,*) 'col_pf%plant_pdemand_vr' 
     write(10,*) col_pf%plant_pdemand_vr
     write(10,*) 'col_pf%plant_pdemand' 
     write(10,*) col_pf%plant_pdemand
     write(10,*) 'veg_pf%biochem_pmin_to_plant' 
     write(10,*) veg_pf%biochem_pmin_to_plant
     write(10,*) 'veg_pf%retransp_to_ppool' 
     write(10,*) veg_pf%retransp_to_ppool
     write(10,*) 'veg_pf%avail_retransp' 
     write(10,*) veg_pf%avail_retransp
     write(10,*) 'veg_pf%plant_pdemand' 
     write(10,*) veg_pf%plant_pdemand
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_cf%o_scalar' 
     write(10,*) col_cf%o_scalar
     write(10,*) 'col_cf%w_scalar' 
     write(10,*) col_cf%w_scalar
     write(10,*) 'col_cf%decomp_k' 
     write(10,*) col_cf%decomp_k
     write(10,*) 'col_cf%t_scalar' 
     write(10,*) col_cf%t_scalar
     write(10,*) 'col_cf%o_scalar' 
     write(10,*) col_cf%o_scalar
     write(10,*) 'col_cf%w_scalar' 
     write(10,*) col_cf%w_scalar
     write(10,*) 'col_cf%decomp_k' 
     write(10,*) col_cf%decomp_k
     write(10,*) 'col_cf%t_scalar' 
     write(10,*) col_cf%t_scalar
     write(10,*) 'decomp_cascade_con%decomp_k_pools' 
     write(10,*) decomp_cascade_con%decomp_k_pools
     write(10,*) 'decomp_cascade_con%decomp_k_pools' 
     write(10,*) decomp_cascade_con%decomp_k_pools
     write(10,*) 'cnstate_vars%scalaravg_col' 
     write(10,*) cnstate_vars%scalaravg_col
     write(10,*) 'cnstate_vars%scalaravg_col' 
     write(10,*) cnstate_vars%scalaravg_col
     write(10,*) 'cnstate_vars%stem_prof_patch' 
     write(10,*) cnstate_vars%stem_prof_patch
     write(10,*) 'cnstate_vars%croot_prof_patch' 
     write(10,*) cnstate_vars%croot_prof_patch
     write(10,*) 'cnstate_vars%pdep_prof_col' 
     write(10,*) cnstate_vars%pdep_prof_col
     write(10,*) 'cnstate_vars%nfixation_prof_col' 
     write(10,*) cnstate_vars%nfixation_prof_col
     write(10,*) 'cnstate_vars%ndep_prof_col' 
     write(10,*) cnstate_vars%ndep_prof_col
     write(10,*) 'cnstate_vars%leaf_prof_patch' 
     write(10,*) cnstate_vars%leaf_prof_patch
     write(10,*) 'cnstate_vars%froot_prof_patch' 
     write(10,*) cnstate_vars%froot_prof_patch
     write(10,*) 'cnstate_vars%grain_flag_patch' 
     write(10,*) cnstate_vars%grain_flag_patch
     write(10,*) 'cnstate_vars%astem_patch' 
     write(10,*) cnstate_vars%astem_patch
     write(10,*) 'cnstate_vars%p_allometry_patch' 
     write(10,*) cnstate_vars%p_allometry_patch
     write(10,*) 'cnstate_vars%tempmax_retransp_patch' 
     write(10,*) cnstate_vars%tempmax_retransp_patch
     write(10,*) 'cnstate_vars%tempmax_retransn_patch' 
     write(10,*) cnstate_vars%tempmax_retransn_patch
     write(10,*) 'cnstate_vars%c_allometry_patch' 
     write(10,*) cnstate_vars%c_allometry_patch
     write(10,*) 'cnstate_vars%n_allometry_patch' 
     write(10,*) cnstate_vars%n_allometry_patch
     write(10,*) 'cnstate_vars%tempsum_potential_gpp_patch' 
     write(10,*) cnstate_vars%tempsum_potential_gpp_patch
     write(10,*) 'cnstate_vars%aleaf_patch' 
     write(10,*) cnstate_vars%aleaf_patch
     write(10,*) 'cnstate_vars%aleafi_patch' 
     write(10,*) cnstate_vars%aleafi_patch
     write(10,*) 'cnstate_vars%astemi_patch' 
     write(10,*) cnstate_vars%astemi_patch
     write(10,*) 'veg_ns%benefit_pgpp_pleafc' 
     write(10,*) veg_ns%benefit_pgpp_pleafc
     if ( use_c14 ) then
     write(10,*) 'c14_veg_cf%psnshade_to_cpool' 
     write(10,*) c14_veg_cf%psnshade_to_cpool
     write(10,*) 'c14_veg_cf%psnsun_to_cpool' 
     write(10,*) c14_veg_cf%psnsun_to_cpool
     endif
     if ( use_c13 ) then
     write(10,*) 'c13_veg_cf%psnshade_to_cpool' 
     write(10,*) c13_veg_cf%psnshade_to_cpool
     write(10,*) 'c13_veg_cf%psnsun_to_cpool' 
     write(10,*) c13_veg_cf%psnsun_to_cpool
     endif
     close(10)
end subroutine 
subroutine update_vars_EcosystemDynNoLeaching2(gpu)
     use ColumnDataType, only : col_nf 
     use ColumnDataType, only : col_pf 
     use ColumnDataType, only : col_cf 
     use clm_instMod, only : cnstate_vars 
     use VegetationDataType, only : veg_nf 
     use VegetationDataType, only : veg_pf 
     use VegetationDataType, only : veg_ns 
     use ColumnDataType, only : col_ns 
     use VegetationType, only : veg_pp 
     use VegetationPropertiesType, only : veg_vp 
     use VegetationDataType, only : veg_cf 
     use clm_instMod, only : canopystate_vars 
     use clm_varctl        , only : use_c13, use_c14
     use VegetationDataType, only : c13_veg_cf 
     use VegetationDataType, only : c14_veg_cf 
     use VegetationDataType, only : veg_es 
     use VegetationDataType, only : veg_cs 
     use VegetationDataType, only : veg_ps 
     use TopounitDataType, only : top_as 
     use VegetationDataType, only : veg_ef 
     use clm_instMod, only : crop_vars 
     use clm_instMod, only : soilstate_vars 
     use ColumnType, only : col_pp 
     use CNDecompCascadeConType, only : decomp_cascade_con 
     use ColumnDataType, only : col_cs 
     use ColumnDataType, only : col_ps 
     use ColumnDataType, only : c13_col_cs 
     use ColumnDataType, only : c14_col_cs 
     use ColumnDataType, only : c13_col_cf 
     use ColumnDataType, only : c14_col_cf 
     use VegetationDataType, only : c14_veg_cs 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_EcosystemDynNoLeaching2.txt"
     else
          file='cpu_EcosystemDynNoLeaching2.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_nf%pmnf_decomp_cascade, & 
     !$acc col_nf%soil_n_immob_flux_vr, & 
     !$acc col_nf%potential_immob_vr, & 
     !$acc col_nf%net_nmin_vr, & 
     !$acc col_nf%soil_n_immob_flux, & 
     !$acc col_nf%gross_nmin_vr, & 
     !$acc col_nf%soil_n_grossmin_flux, & 
     !$acc col_nf%decomp_cascade_ntransfer_vr, & 
     !$acc col_nf%decomp_cascade_sminn_flux_vr, & 
     !$acc col_nf%sminn_to_denit_decomp_cascade_vr, & 
     !$acc col_nf%fmax_denit_nitrate_vr, & 
     !$acc col_nf%r_psi, & 
     !$acc col_nf%k_nitr_t_vr, & 
     !$acc col_nf%diffus, & 
     !$acc col_nf%fr_WFPS, & 
     !$acc col_nf%fmax_denit_carbonsubstrate_vr, & 
     !$acc col_nf%ratio_no3_co2, & 
     !$acc col_nf%pot_f_denit_vr, & 
     !$acc col_nf%pot_f_nit_vr, & 
     !$acc col_nf%ratio_k1, & 
     !$acc col_nf%n2_n2o_ratio_denit_vr, & 
     !$acc col_nf%f_denit_base_vr, & 
     !$acc col_nf%soil_bulkdensity, & 
     !$acc col_nf%k_nitr_h2o_vr, & 
     !$acc col_nf%k_nitr_ph_vr, & 
     !$acc col_nf%anaerobic_frac, & 
     !$acc col_nf%soil_co2_prod, & 
     !$acc col_nf%k_nitr_vr, & 
     !$acc col_nf%wfps_vr, & 
     !$acc col_nf%smin_no3_massdens_vr, & 
     !$acc col_nf%actual_immob, & 
     !$acc col_nf%actual_immob_nh4_vr, & 
     !$acc col_nf%actual_immob_no3_vr, & 
     !$acc col_nf%f_n2o_nit_vr, & 
     !$acc col_nf%supplement_to_sminn_vr, & 
     !$acc col_nf%actual_immob_vr, & 
     !$acc col_nf%col_plant_no3demand_vr, & 
     !$acc col_nf%actual_immob_no3, & 
     !$acc col_nf%actual_immob_nh4, & 
     !$acc col_nf%sminn_to_plant, & 
     !$acc col_nf%sminn_to_denit_excess_vr, & 
     !$acc col_nf%smin_no3_to_plant_vr, & 
     !$acc col_nf%f_denit_vr, & 
     !$acc col_nf%f_nit_vr, & 
     !$acc col_nf%col_plant_nh4demand_vr, & 
     !$acc col_nf%f_n2o_denit_vr, & 
     !$acc col_nf%sminn_to_plant_vr, & 
     !$acc col_nf%potential_immob, & 
     !$acc col_nf%col_plant_ndemand_vr, & 
     !$acc col_nf%smin_nh4_to_plant_vr, & 
     !$acc col_nf%actual_immob, & 
     !$acc col_nf%gross_nmin, & 
     !$acc col_nf%sminn_to_plant, & 
     !$acc col_nf%net_nmin, & 
     !$acc col_nf%potential_immob, & 
     !$acc col_nf%sminn_to_plant, & 
     !$acc col_nf%smin_no3_to_plant_vr, & 
     !$acc col_nf%plant_n_uptake_flux, & 
     !$acc col_nf%sminn_to_plant_vr, & 
     !$acc col_nf%smin_nh4_to_plant_vr, & 
     !$acc col_nf%hrv_cropn_to_prod1n, & 
     !$acc col_nf%phenology_n_to_litr_cel_n, & 
     !$acc col_nf%phenology_n_to_litr_lig_n, & 
     !$acc col_nf%phenology_n_to_litr_met_n, & 
     !$acc col_nf%decomp_npools_sourcesink, & 
     !$acc col_nf%decomp_npools_transport_tendency, & 
     !$acc col_nf%gap_mortality_n_to_litr_met_n, & 
     !$acc col_nf%gap_mortality_n_to_litr_lig_n, & 
     !$acc col_nf%gap_mortality_n_to_litr_cel_n, & 
     !$acc col_nf%gap_mortality_n_to_cwdn, & 
     !$acc col_nf%hrv_deadstemn_to_prod100n, & 
     !$acc col_nf%hrv_deadstemn_to_prod10n, & 
     !$acc col_nf%harvest_n_to_litr_cel_n, & 
     !$acc col_nf%harvest_n_to_litr_met_n, & 
     !$acc col_nf%harvest_n_to_litr_lig_n, & 
     !$acc col_nf%harvest_n_to_cwdn, & 
     !$acc col_nf%prod10n_loss, & 
     !$acc col_nf%prod100n_loss, & 
     !$acc col_nf%prod1n_loss, & 
     !$acc col_nf%m_n_to_litr_cel_fire, & 
     !$acc col_nf%fire_mortality_n_to_cwdn, & 
     !$acc col_nf%m_decomp_npools_to_fire_vr, & 
     !$acc col_nf%m_n_to_litr_lig_fire, & 
     !$acc col_nf%m_n_to_litr_met_fire, & 
     !$acc col_nf%decomp_npools_deposit, & 
     !$acc col_nf%decomp_npools_yield_vr, & 
     !$acc col_nf%decomp_npools_erode )
     !$acc update self(& 
     !$acc col_pf%pmpf_decomp_cascade, & 
     !$acc col_pf%soil_p_immob_flux_vr, & 
     !$acc col_pf%potential_immob_p_vr, & 
     !$acc col_pf%net_pmin_vr, & 
     !$acc col_pf%decomp_cascade_sminp_flux_vr, & 
     !$acc col_pf%gross_pmin_vr, & 
     !$acc col_pf%soil_p_grossmin_flux, & 
     !$acc col_pf%decomp_cascade_ptransfer_vr, & 
     !$acc col_pf%soil_p_immob_flux, & 
     !$acc col_pf%sminp_to_plant_vr, & 
     !$acc col_pf%sminp_to_plant, & 
     !$acc col_pf%adsorb_to_labilep_vr, & 
     !$acc col_pf%actual_immob_p, & 
     !$acc col_pf%desorb_to_solutionp_vr, & 
     !$acc col_pf%col_plant_pdemand_vr, & 
     !$acc col_pf%actual_immob_p_vr, & 
     !$acc col_pf%potential_immob_p, & 
     !$acc col_pf%supplement_to_sminp_vr, & 
     !$acc col_pf%gross_pmin, & 
     !$acc col_pf%net_pmin, & 
     !$acc col_pf%sminp_to_plant, & 
     !$acc col_pf%plant_p_uptake_flux, & 
     !$acc col_pf%sminp_to_plant_vr, & 
     !$acc col_pf%hrv_cropp_to_prod1p, & 
     !$acc col_pf%phenology_p_to_litr_cel_p, & 
     !$acc col_pf%phenology_p_to_litr_lig_p, & 
     !$acc col_pf%phenology_p_to_litr_met_p, & 
     !$acc col_pf%decomp_ppools_sourcesink, & 
     !$acc col_pf%decomp_ppools_transport_tendency, & 
     !$acc col_pf%gap_mortality_p_to_litr_met_p, & 
     !$acc col_pf%gap_mortality_p_to_litr_cel_p, & 
     !$acc col_pf%gap_mortality_p_to_litr_lig_p, & 
     !$acc col_pf%gap_mortality_p_to_cwdp, & 
     !$acc col_pf%hrv_deadstemp_to_prod10p, & 
     !$acc col_pf%harvest_p_to_cwdp, & 
     !$acc col_pf%hrv_deadstemp_to_prod100p, & 
     !$acc col_pf%harvest_p_to_litr_cel_p, & 
     !$acc col_pf%harvest_p_to_litr_lig_p, & 
     !$acc col_pf%harvest_p_to_litr_met_p, & 
     !$acc col_pf%prod10p_loss, & 
     !$acc col_pf%prod100p_loss, & 
     !$acc col_pf%prod1p_loss, & 
     !$acc col_pf%m_p_to_litr_lig_fire, & 
     !$acc col_pf%fire_mortality_p_to_cwdp, & 
     !$acc col_pf%m_p_to_litr_met_fire, & 
     !$acc col_pf%m_decomp_ppools_to_fire_vr, & 
     !$acc col_pf%m_p_to_litr_cel_fire, & 
     !$acc col_pf%secondp_deposit, & 
     !$acc col_pf%primp_yield_vr, & 
     !$acc col_pf%occlp_erode, & 
     !$acc col_pf%decomp_ppools_deposit, & 
     !$acc col_pf%occlp_yield_vr, & 
     !$acc col_pf%decomp_ppools_yield_vr, & 
     !$acc col_pf%labilep_yield_vr, & 
     !$acc col_pf%primp_deposit, & 
     !$acc col_pf%labilep_erode, & 
     !$acc col_pf%occlp_deposit, & 
     !$acc col_pf%labilep_deposit, & 
     !$acc col_pf%secondp_erode, & 
     !$acc col_pf%secondp_yield_vr, & 
     !$acc col_pf%decomp_ppools_erode, & 
     !$acc col_pf%primp_erode )
     !$acc update self(& 
     !$acc col_cf%fphr, & 
     !$acc col_cf%decomp_cascade_hr_vr, & 
     !$acc col_cf%phr_vr, & 
     !$acc col_cf%decomp_cascade_ctransfer_vr, & 
     !$acc col_cf%fphr, & 
     !$acc col_cf%phr_vr, & 
     !$acc col_cf%hrv_cropc_to_prod1c, & 
     !$acc col_cf%phenology_c_to_litr_lig_c, & 
     !$acc col_cf%phenology_c_to_litr_cel_c, & 
     !$acc col_cf%phenology_c_to_litr_met_c, & 
     !$acc col_cf%phenology_c_to_litr_lig_c, & 
     !$acc col_cf%phenology_c_to_litr_cel_c, & 
     !$acc col_cf%phenology_c_to_litr_met_c, & 
     !$acc col_cf%decomp_cpools_sourcesink, & 
     !$acc col_cf%decomp_cpools_transport_tendency, & 
     !$acc col_cf%gap_mortality_c_to_litr_lig_c, & 
     !$acc col_cf%gap_mortality_c_to_litr_cel_c, & 
     !$acc col_cf%gap_mortality_c_to_cwdc, & 
     !$acc col_cf%gap_mortality_c_to_litr_met_c, & 
     !$acc col_cf%gap_mortality_c_to_litr_lig_c, & 
     !$acc col_cf%gap_mortality_c_to_litr_cel_c, & 
     !$acc col_cf%gap_mortality_c_to_cwdc, & 
     !$acc col_cf%gap_mortality_c_to_litr_met_c, & 
     !$acc col_cf%harvest_c_to_litr_lig_c, & 
     !$acc col_cf%hrv_deadstemc_to_prod10c, & 
     !$acc col_cf%hrv_deadstemc_to_prod100c, & 
     !$acc col_cf%harvest_c_to_litr_cel_c, & 
     !$acc col_cf%harvest_c_to_cwdc, & 
     !$acc col_cf%harvest_c_to_litr_met_c, & 
     !$acc col_cf%harvest_c_to_litr_lig_c, & 
     !$acc col_cf%hrv_deadstemc_to_prod10c, & 
     !$acc col_cf%hrv_deadstemc_to_prod100c, & 
     !$acc col_cf%harvest_c_to_litr_cel_c, & 
     !$acc col_cf%harvest_c_to_cwdc, & 
     !$acc col_cf%harvest_c_to_litr_met_c, & 
     !$acc col_cf%prod10c_loss, & 
     !$acc col_cf%prod100c_loss, & 
     !$acc col_cf%prod1c_loss, & 
     !$acc col_cf%fire_mortality_c_to_cwdc, & 
     !$acc col_cf%m_c_to_litr_met_fire, & 
     !$acc col_cf%m_c_to_litr_cel_fire, & 
     !$acc col_cf%m_c_to_litr_lig_fire, & 
     !$acc col_cf%m_decomp_cpools_to_fire_vr, & 
     !$acc col_cf%somc_fire, & 
     !$acc col_cf%decomp_cpools_yield_vr, & 
     !$acc col_cf%decomp_cpools_deposit, & 
     !$acc col_cf%decomp_cpools_erode )
     !$acc update self(& 
     !$acc cnstate_vars%fpg_p_vr_col, & 
     !$acc cnstate_vars%cn_scalar, & 
     !$acc cnstate_vars%fpi_p_vr_col, & 
     !$acc cnstate_vars%fpg_p_col, & 
     !$acc cnstate_vars%fpi_vr_col, & 
     !$acc cnstate_vars%fpg_nh4_vr_col, & 
     !$acc cnstate_vars%fpg_col, & 
     !$acc cnstate_vars%cp_scalar, & 
     !$acc cnstate_vars%fpi_p_col, & 
     !$acc cnstate_vars%fpi_col, & 
     !$acc cnstate_vars%fpg_no3_vr_col, & 
     !$acc cnstate_vars%fpg_vr_col, & 
     !$acc cnstate_vars%fpi_vr_col, & 
     !$acc cnstate_vars%fpg_col, & 
     !$acc cnstate_vars%fpi_col, & 
     !$acc cnstate_vars%downreg_patch, & 
     !$acc cnstate_vars%n_allometry_patch, & 
     !$acc cnstate_vars%p_allometry_patch, & 
     !$acc cnstate_vars%c_allometry_patch, & 
     !$acc cnstate_vars%tempavg_t2m_patch, & 
     !$acc cnstate_vars%bglfr_froot_patch, & 
     !$acc cnstate_vars%lgsf_patch, & 
     !$acc cnstate_vars%bglfr_leaf_patch, & 
     !$acc cnstate_vars%bgtr_patch, & 
     !$acc cnstate_vars%lgsf_patch, & 
     !$acc cnstate_vars%onset_flag_patch, & 
     !$acc cnstate_vars%days_active_patch, & 
     !$acc cnstate_vars%bgtr_patch, & 
     !$acc cnstate_vars%onset_counter_patch, & 
     !$acc cnstate_vars%onset_gdd_patch, & 
     !$acc cnstate_vars%offset_counter_patch, & 
     !$acc cnstate_vars%onset_gddflag_patch, & 
     !$acc cnstate_vars%bglfr_froot_patch, & 
     !$acc cnstate_vars%dormant_flag_patch, & 
     !$acc cnstate_vars%bglfr_leaf_patch, & 
     !$acc cnstate_vars%offset_flag_patch, & 
     !$acc cnstate_vars%offset_fdd_patch, & 
     !$acc cnstate_vars%lgsf_patch, & 
     !$acc cnstate_vars%onset_flag_patch, & 
     !$acc cnstate_vars%days_active_patch, & 
     !$acc cnstate_vars%bgtr_patch, & 
     !$acc cnstate_vars%onset_counter_patch, & 
     !$acc cnstate_vars%onset_swi_patch, & 
     !$acc cnstate_vars%onset_gdd_patch, & 
     !$acc cnstate_vars%offset_counter_patch, & 
     !$acc cnstate_vars%onset_gddflag_patch, & 
     !$acc cnstate_vars%offset_swi_patch, & 
     !$acc cnstate_vars%bglfr_froot_patch, & 
     !$acc cnstate_vars%dormant_flag_patch, & 
     !$acc cnstate_vars%bglfr_leaf_patch, & 
     !$acc cnstate_vars%onset_fdd_patch, & 
     !$acc cnstate_vars%offset_flag_patch, & 
     !$acc cnstate_vars%lgsf_patch, & 
     !$acc cnstate_vars%huigrain_patch, & 
     !$acc cnstate_vars%bgtr_patch, & 
     !$acc cnstate_vars%huileaf_patch, & 
     !$acc cnstate_vars%onset_flag_patch, & 
     !$acc cnstate_vars%onset_counter_patch, & 
     !$acc cnstate_vars%idop_patch, & 
     !$acc cnstate_vars%bglfr_froot_patch, & 
     !$acc cnstate_vars%gddmaturity_patch, & 
     !$acc cnstate_vars%offset_counter_patch, & 
     !$acc cnstate_vars%hdidx_patch, & 
     !$acc cnstate_vars%bglfr_leaf_patch, & 
     !$acc cnstate_vars%offset_flag_patch, & 
     !$acc cnstate_vars%cumvd_patch, & 
     !$acc cnstate_vars%hdidx_patch, & 
     !$acc cnstate_vars%huigrain_patch, & 
     !$acc cnstate_vars%cumvd_patch, & 
     !$acc cnstate_vars%gddmaturity_patch, & 
     !$acc cnstate_vars%r_mort_cal_patch, & 
     !$acc cnstate_vars%nfire_col, & 
     !$acc cnstate_vars%lgdp1_col, & 
     !$acc cnstate_vars%lfc_col, & 
     !$acc cnstate_vars%dtrotr_col, & 
     !$acc cnstate_vars%lfwt_col, & 
     !$acc cnstate_vars%burndate_patch, & 
     !$acc cnstate_vars%trotr1_col, & 
     !$acc cnstate_vars%fsr_col, & 
     !$acc cnstate_vars%baf_crop_col, & 
     !$acc cnstate_vars%baf_peatf_col, & 
     !$acc cnstate_vars%fbac1_col, & 
     !$acc cnstate_vars%cropf_col, & 
     !$acc cnstate_vars%fd_col, & 
     !$acc cnstate_vars%lpop_col, & 
     !$acc cnstate_vars%farea_burned_col, & 
     !$acc cnstate_vars%trotr2_col, & 
     !$acc cnstate_vars%fbac_col, & 
     !$acc cnstate_vars%wtlf_col, & 
     !$acc cnstate_vars%lgdp_col, & 
     !$acc cnstate_vars%lfc_col, & 
     !$acc cnstate_vars%lfc2_col, & 
     !$acc cnstate_vars%rc14_atm_patch )
     !$acc update self(& 
     !$acc veg_nf%smin_no3_to_plant, & 
     !$acc veg_nf%plant_ndemand_vr, & 
     !$acc veg_nf%sminn_to_plant, & 
     !$acc veg_nf%plant_no3demand_vr, & 
     !$acc veg_nf%smin_nh4_to_plant, & 
     !$acc veg_nf%plant_nh4demand_vr, & 
     !$acc veg_nf%npool_to_frootn, & 
     !$acc veg_nf%retransn_to_npool, & 
     !$acc veg_nf%npool_to_livecrootn, & 
     !$acc veg_nf%npool_to_livecrootn_storage, & 
     !$acc veg_nf%npool_to_deadcrootn_storage, & 
     !$acc veg_nf%npool_to_livestemn_storage, & 
     !$acc veg_nf%npool_to_leafn_storage, & 
     !$acc veg_nf%npool_to_frootn_storage, & 
     !$acc veg_nf%plant_nalloc, & 
     !$acc veg_nf%npool_to_grainn, & 
     !$acc veg_nf%npool_to_leafn, & 
     !$acc veg_nf%npool_to_deadstemn, & 
     !$acc veg_nf%npool_to_deadstemn_storage, & 
     !$acc veg_nf%npool_to_deadcrootn, & 
     !$acc veg_nf%supplement_to_plantn, & 
     !$acc veg_nf%sminn_to_npool, & 
     !$acc veg_nf%npool_to_livestemn, & 
     !$acc veg_nf%avail_retransn, & 
     !$acc veg_nf%npool_to_grainn_storage, & 
     !$acc veg_nf%leafn_storage_to_xfer, & 
     !$acc veg_nf%livestemn_storage_to_xfer, & 
     !$acc veg_nf%deadcrootn_xfer_to_deadcrootn, & 
     !$acc veg_nf%frootn_xfer_to_frootn, & 
     !$acc veg_nf%frootn_storage_to_xfer, & 
     !$acc veg_nf%deadstemn_xfer_to_deadstemn, & 
     !$acc veg_nf%deadcrootn_storage_to_xfer, & 
     !$acc veg_nf%leafn_xfer_to_leafn, & 
     !$acc veg_nf%deadstemn_storage_to_xfer, & 
     !$acc veg_nf%livecrootn_xfer_to_livecrootn, & 
     !$acc veg_nf%livestemn_xfer_to_livestemn, & 
     !$acc veg_nf%livecrootn_storage_to_xfer, & 
     !$acc veg_nf%leafn_storage_to_xfer, & 
     !$acc veg_nf%livestemn_storage_to_xfer, & 
     !$acc veg_nf%deadcrootn_xfer_to_deadcrootn, & 
     !$acc veg_nf%frootn_xfer_to_frootn, & 
     !$acc veg_nf%frootn_storage_to_xfer, & 
     !$acc veg_nf%deadstemn_xfer_to_deadstemn, & 
     !$acc veg_nf%deadcrootn_storage_to_xfer, & 
     !$acc veg_nf%leafn_xfer_to_leafn, & 
     !$acc veg_nf%deadstemn_storage_to_xfer, & 
     !$acc veg_nf%livecrootn_xfer_to_livecrootn, & 
     !$acc veg_nf%livestemn_xfer_to_livestemn, & 
     !$acc veg_nf%livecrootn_storage_to_xfer, & 
     !$acc veg_nf%fert_counter, & 
     !$acc veg_nf%crop_seedn_to_leaf, & 
     !$acc veg_nf%fert, & 
     !$acc veg_nf%deadcrootn_xfer_to_deadcrootn, & 
     !$acc veg_nf%frootn_xfer_to_frootn, & 
     !$acc veg_nf%deadstemn_xfer_to_deadstemn, & 
     !$acc veg_nf%leafn_xfer_to_leafn, & 
     !$acc veg_nf%livecrootn_xfer_to_livecrootn, & 
     !$acc veg_nf%livestemn_xfer_to_livestemn, & 
     !$acc veg_nf%hrv_livestemn_to_prod1n, & 
     !$acc veg_nf%hrv_grainn_to_prod1n, & 
     !$acc veg_nf%hrv_leafn_to_prod1n, & 
     !$acc veg_nf%hrv_cropn_to_prod1n, & 
     !$acc veg_nf%leafn_to_retransn, & 
     !$acc veg_nf%frootn_to_litter, & 
     !$acc veg_nf%leafn_to_litter, & 
     !$acc veg_nf%livestemn_to_litter, & 
     !$acc veg_nf%frootn_to_litter, & 
     !$acc veg_nf%leafn_to_litter, & 
     !$acc veg_nf%leafn_to_retransn, & 
     !$acc veg_nf%livecrootn_to_retransn, & 
     !$acc veg_nf%livestemn_to_deadstemn, & 
     !$acc veg_nf%livecrootn_to_deadcrootn, & 
     !$acc veg_nf%livestemn_to_retransn, & 
     !$acc veg_nf%npool_to_livecrootn, & 
     !$acc veg_nf%deadcrootn_xfer_to_deadcrootn, & 
     !$acc veg_nf%livestemn_to_retransn, & 
     !$acc veg_nf%npool_to_frootn_storage, & 
     !$acc veg_nf%livecrootn_xfer_to_livecrootn, & 
     !$acc veg_nf%grainn_to_food, & 
     !$acc veg_nf%npool_to_grainn_storage, & 
     !$acc veg_nf%npool_to_deadcrootn_storage, & 
     !$acc veg_nf%npool_to_deadstemn_storage, & 
     !$acc veg_nf%npool_to_deadcrootn, & 
     !$acc veg_nf%livestemn_to_litter, & 
     !$acc veg_nf%deadstemn_storage_to_xfer, & 
     !$acc veg_nf%livecrootn_to_retransn, & 
     !$acc veg_nf%livestemn_storage_to_xfer, & 
     !$acc veg_nf%npool_to_livestemn_storage, & 
     !$acc veg_nf%npool_to_leafn_storage, & 
     !$acc veg_nf%npool_to_livestemn, & 
     !$acc veg_nf%grainn_xfer_to_grainn, & 
     !$acc veg_nf%livecrootn_to_deadcrootn, & 
     !$acc veg_nf%deadstemn_xfer_to_deadstemn, & 
     !$acc veg_nf%deadcrootn_storage_to_xfer, & 
     !$acc veg_nf%livestemn_xfer_to_livestemn, & 
     !$acc veg_nf%livestemn_to_deadstemn, & 
     !$acc veg_nf%npool_to_frootn, & 
     !$acc veg_nf%npool_to_livecrootn_storage, & 
     !$acc veg_nf%npool_to_grainn, & 
     !$acc veg_nf%frootn_to_retransn, & 
     !$acc veg_nf%npool_to_deadstemn, & 
     !$acc veg_nf%livecrootn_storage_to_xfer, & 
     !$acc veg_nf%npool_to_leafn, & 
     !$acc veg_nf%m_deadstemn_storage_to_litter, & 
     !$acc veg_nf%m_deadstemn_xfer_to_litter, & 
     !$acc veg_nf%m_livecrootn_to_litter, & 
     !$acc veg_nf%m_leafn_to_litter, & 
     !$acc veg_nf%m_deadcrootn_storage_to_litter, & 
     !$acc veg_nf%m_deadstemn_to_litter, & 
     !$acc veg_nf%m_npool_to_litter, & 
     !$acc veg_nf%m_frootn_storage_to_litter, & 
     !$acc veg_nf%m_deadcrootn_xfer_to_litter, & 
     !$acc veg_nf%m_livestemn_xfer_to_litter, & 
     !$acc veg_nf%m_livecrootn_storage_to_litter, & 
     !$acc veg_nf%m_livestemn_to_litter, & 
     !$acc veg_nf%m_frootn_to_litter, & 
     !$acc veg_nf%m_frootn_xfer_to_litter, & 
     !$acc veg_nf%m_retransn_to_litter, & 
     !$acc veg_nf%m_livecrootn_xfer_to_litter, & 
     !$acc veg_nf%m_deadcrootn_to_litter, & 
     !$acc veg_nf%m_leafn_xfer_to_litter, & 
     !$acc veg_nf%m_livestemn_storage_to_litter, & 
     !$acc veg_nf%m_leafn_storage_to_litter, & 
     !$acc veg_nf%m_deadstemn_storage_to_litter, & 
     !$acc veg_nf%m_leafn_to_litter, & 
     !$acc veg_nf%m_deadcrootn_storage_to_litter, & 
     !$acc veg_nf%m_frootn_storage_to_litter, & 
     !$acc veg_nf%m_frootn_to_litter, & 
     !$acc veg_nf%m_livecrootn_storage_to_litter, & 
     !$acc veg_nf%m_leafn_xfer_to_litter, & 
     !$acc veg_nf%m_retransn_to_litter, & 
     !$acc veg_nf%m_livestemn_storage_to_litter, & 
     !$acc veg_nf%m_frootn_xfer_to_litter, & 
     !$acc veg_nf%m_leafn_storage_to_litter, & 
     !$acc veg_nf%hrv_leafn_to_litter, & 
     !$acc veg_nf%hrv_frootn_xfer_to_litter, & 
     !$acc veg_nf%hrv_livestemn_xfer_to_litter, & 
     !$acc veg_nf%hrv_deadstemn_storage_to_litter, & 
     !$acc veg_nf%hrv_deadcrootn_xfer_to_litter, & 
     !$acc veg_nf%hrv_deadstemn_to_prod100n, & 
     !$acc veg_nf%hrv_livecrootn_storage_to_litter, & 
     !$acc veg_nf%hrv_deadcrootn_to_litter, & 
     !$acc veg_nf%hrv_livestemn_to_litter, & 
     !$acc veg_nf%hrv_retransn_to_litter, & 
     !$acc veg_nf%hrv_livecrootn_to_litter, & 
     !$acc veg_nf%hrv_npool_to_litter, & 
     !$acc veg_nf%hrv_deadcrootn_storage_to_litter, & 
     !$acc veg_nf%hrv_deadstemn_xfer_to_litter, & 
     !$acc veg_nf%hrv_livestemn_storage_to_litter, & 
     !$acc veg_nf%hrv_deadstemn_to_prod10n, & 
     !$acc veg_nf%hrv_leafn_xfer_to_litter, & 
     !$acc veg_nf%hrv_livecrootn_xfer_to_litter, & 
     !$acc veg_nf%hrv_frootn_to_litter, & 
     !$acc veg_nf%hrv_leafn_storage_to_litter, & 
     !$acc veg_nf%hrv_frootn_storage_to_litter, & 
     !$acc veg_nf%m_livestemn_to_deadstemn_fire, & 
     !$acc veg_nf%m_npool_to_litter_fire, & 
     !$acc veg_nf%m_leafn_to_fire, & 
     !$acc veg_nf%m_frootn_xfer_to_litter_fire, & 
     !$acc veg_nf%m_livecrootn_to_fire, & 
     !$acc veg_nf%m_leafn_storage_to_litter_fire, & 
     !$acc veg_nf%m_retransn_to_litter_fire, & 
     !$acc veg_nf%m_livestemn_xfer_to_fire, & 
     !$acc veg_nf%m_livestemn_storage_to_fire, & 
     !$acc veg_nf%m_deadcrootn_xfer_to_fire, & 
     !$acc veg_nf%m_frootn_xfer_to_fire, & 
     !$acc veg_nf%m_retransn_to_fire, & 
     !$acc veg_nf%m_livestemn_to_fire, & 
     !$acc veg_nf%m_deadstemn_xfer_to_fire, & 
     !$acc veg_nf%m_livestemn_storage_to_litter_fire, & 
     !$acc veg_nf%m_leafn_storage_to_fire, & 
     !$acc veg_nf%m_leafn_to_litter_fire, & 
     !$acc veg_nf%m_livestemn_to_litter_fire, & 
     !$acc veg_nf%m_deadstemn_storage_to_fire, & 
     !$acc veg_nf%m_livecrootn_xfer_to_litter_fire, & 
     !$acc veg_nf%m_leafn_xfer_to_litter_fire, & 
     !$acc veg_nf%m_npool_to_fire, & 
     !$acc veg_nf%m_frootn_to_fire, & 
     !$acc veg_nf%m_livecrootn_xfer_to_fire, & 
     !$acc veg_nf%m_leafn_xfer_to_fire, & 
     !$acc veg_nf%m_deadstemn_to_fire, & 
     !$acc veg_nf%m_deadstemn_xfer_to_litter_fire, & 
     !$acc veg_nf%m_deadcrootn_storage_to_litter_fire, & 
     !$acc veg_nf%m_deadstemn_to_litter_fire, & 
     !$acc veg_nf%m_deadcrootn_storage_to_fire, & 
     !$acc veg_nf%m_frootn_storage_to_fire, & 
     !$acc veg_nf%m_livecrootn_to_litter_fire, & 
     !$acc veg_nf%m_livecrootn_storage_to_litter_fire, & 
     !$acc veg_nf%m_deadcrootn_to_litter_fire, & 
     !$acc veg_nf%m_frootn_to_litter_fire, & 
     !$acc veg_nf%m_livestemn_xfer_to_litter_fire, & 
     !$acc veg_nf%m_livecrootn_to_deadcrootn_fire, & 
     !$acc veg_nf%m_frootn_storage_to_litter_fire, & 
     !$acc veg_nf%m_deadcrootn_to_fire, & 
     !$acc veg_nf%m_deadstemn_storage_to_litter_fire, & 
     !$acc veg_nf%m_livecrootn_storage_to_fire, & 
     !$acc veg_nf%m_deadcrootn_xfer_to_litter_fire )
     !$acc update self(& 
     !$acc veg_pf%sminp_to_plant, & 
     !$acc veg_pf%plant_pdemand_vr, & 
     !$acc veg_pf%ppool_to_grainp, & 
     !$acc veg_pf%ppool_to_livestemp_storage, & 
     !$acc veg_pf%ppool_to_leafp, & 
     !$acc veg_pf%plant_palloc, & 
     !$acc veg_pf%ppool_to_livecrootp, & 
     !$acc veg_pf%ppool_to_deadcrootp_storage, & 
     !$acc veg_pf%avail_retransp, & 
     !$acc veg_pf%ppool_to_frootp_storage, & 
     !$acc veg_pf%ppool_to_deadstemp_storage, & 
     !$acc veg_pf%retransp_to_ppool, & 
     !$acc veg_pf%ppool_to_deadstemp, & 
     !$acc veg_pf%ppool_to_grainp_storage, & 
     !$acc veg_pf%ppool_to_leafp_storage, & 
     !$acc veg_pf%supplement_to_plantp, & 
     !$acc veg_pf%ppool_to_deadcrootp, & 
     !$acc veg_pf%sminp_to_ppool, & 
     !$acc veg_pf%ppool_to_livecrootp_storage, & 
     !$acc veg_pf%ppool_to_frootp, & 
     !$acc veg_pf%ppool_to_livestemp, & 
     !$acc veg_pf%frootp_xfer_to_frootp, & 
     !$acc veg_pf%deadstemp_xfer_to_deadstemp, & 
     !$acc veg_pf%leafp_xfer_to_leafp, & 
     !$acc veg_pf%deadstemp_storage_to_xfer, & 
     !$acc veg_pf%livestemp_xfer_to_livestemp, & 
     !$acc veg_pf%deadcrootp_xfer_to_deadcrootp, & 
     !$acc veg_pf%livestemp_storage_to_xfer, & 
     !$acc veg_pf%leafp_storage_to_xfer, & 
     !$acc veg_pf%frootp_storage_to_xfer, & 
     !$acc veg_pf%livecrootp_xfer_to_livecrootp, & 
     !$acc veg_pf%deadcrootp_storage_to_xfer, & 
     !$acc veg_pf%livecrootp_storage_to_xfer, & 
     !$acc veg_pf%frootp_xfer_to_frootp, & 
     !$acc veg_pf%deadstemp_xfer_to_deadstemp, & 
     !$acc veg_pf%leafp_xfer_to_leafp, & 
     !$acc veg_pf%deadstemp_storage_to_xfer, & 
     !$acc veg_pf%livestemp_xfer_to_livestemp, & 
     !$acc veg_pf%deadcrootp_xfer_to_deadcrootp, & 
     !$acc veg_pf%livestemp_storage_to_xfer, & 
     !$acc veg_pf%leafp_storage_to_xfer, & 
     !$acc veg_pf%frootp_storage_to_xfer, & 
     !$acc veg_pf%livecrootp_xfer_to_livecrootp, & 
     !$acc veg_pf%deadcrootp_storage_to_xfer, & 
     !$acc veg_pf%livecrootp_storage_to_xfer, & 
     !$acc veg_pf%crop_seedp_to_leaf, & 
     !$acc veg_pf%frootp_xfer_to_frootp, & 
     !$acc veg_pf%deadstemp_xfer_to_deadstemp, & 
     !$acc veg_pf%leafp_xfer_to_leafp, & 
     !$acc veg_pf%livestemp_xfer_to_livestemp, & 
     !$acc veg_pf%deadcrootp_xfer_to_deadcrootp, & 
     !$acc veg_pf%livecrootp_xfer_to_livecrootp, & 
     !$acc veg_pf%hrv_livestemp_to_prod1p, & 
     !$acc veg_pf%hrv_grainp_to_prod1p, & 
     !$acc veg_pf%hrv_leafp_to_prod1p, & 
     !$acc veg_pf%hrv_cropp_to_prod1p, & 
     !$acc veg_pf%leafp_to_litter, & 
     !$acc veg_pf%frootp_to_litter, & 
     !$acc veg_pf%livestemp_to_litter, & 
     !$acc veg_pf%leafp_to_retransp, & 
     !$acc veg_pf%leafp_to_litter, & 
     !$acc veg_pf%frootp_to_litter, & 
     !$acc veg_pf%leafp_to_retransp, & 
     !$acc veg_pf%livestemp_to_deadstemp, & 
     !$acc veg_pf%livecrootp_to_deadcrootp, & 
     !$acc veg_pf%livestemp_to_retransp, & 
     !$acc veg_pf%livecrootp_to_retransp, & 
     !$acc veg_pf%ppool_to_livestemp_storage, & 
     !$acc veg_pf%deadstemp_storage_to_xfer, & 
     !$acc veg_pf%livestemp_to_deadstemp, & 
     !$acc veg_pf%retransp_to_ppool, & 
     !$acc veg_pf%livecrootp_to_retransp, & 
     !$acc veg_pf%livecrootp_xfer_to_livecrootp, & 
     !$acc veg_pf%ppool_to_grainp, & 
     !$acc veg_pf%livestemp_xfer_to_livestemp, & 
     !$acc veg_pf%ppool_to_deadcrootp_storage, & 
     !$acc veg_pf%ppool_to_grainp_storage, & 
     !$acc veg_pf%ppool_to_leafp_storage, & 
     !$acc veg_pf%leafp_to_retransp, & 
     !$acc veg_pf%ppool_to_deadcrootp, & 
     !$acc veg_pf%livestemp_to_litter, & 
     !$acc veg_pf%livecrootp_to_deadcrootp, & 
     !$acc veg_pf%deadstemp_xfer_to_deadstemp, & 
     !$acc veg_pf%frootp_to_litter, & 
     !$acc veg_pf%ppool_to_frootp_storage, & 
     !$acc veg_pf%grainp_xfer_to_grainp, & 
     !$acc veg_pf%frootp_to_retransp, & 
     !$acc veg_pf%livestemp_to_retransp, & 
     !$acc veg_pf%supplement_to_plantp, & 
     !$acc veg_pf%livecrootp_storage_to_xfer, & 
     !$acc veg_pf%ppool_to_livestemp, & 
     !$acc veg_pf%frootp_xfer_to_frootp, & 
     !$acc veg_pf%leafp_xfer_to_leafp, & 
     !$acc veg_pf%ppool_to_leafp, & 
     !$acc veg_pf%ppool_to_livecrootp, & 
     !$acc veg_pf%grainp_to_food, & 
     !$acc veg_pf%deadcrootp_xfer_to_deadcrootp, & 
     !$acc veg_pf%livestemp_storage_to_xfer, & 
     !$acc veg_pf%ppool_to_deadstemp_storage, & 
     !$acc veg_pf%leafp_to_litter, & 
     !$acc veg_pf%ppool_to_deadstemp, & 
     !$acc veg_pf%frootp_storage_to_xfer, & 
     !$acc veg_pf%ppool_to_livecrootp_storage, & 
     !$acc veg_pf%ppool_to_frootp, & 
     !$acc veg_pf%deadcrootp_storage_to_xfer, & 
     !$acc veg_pf%m_frootp_storage_to_litter, & 
     !$acc veg_pf%m_frootp_to_litter, & 
     !$acc veg_pf%m_deadstemp_to_litter, & 
     !$acc veg_pf%m_livecrootp_storage_to_litter, & 
     !$acc veg_pf%m_livestemp_xfer_to_litter, & 
     !$acc veg_pf%m_leafp_storage_to_litter, & 
     !$acc veg_pf%m_deadstemp_xfer_to_litter, & 
     !$acc veg_pf%m_livestemp_storage_to_litter, & 
     !$acc veg_pf%m_livecrootp_to_litter, & 
     !$acc veg_pf%m_deadstemp_storage_to_litter, & 
     !$acc veg_pf%m_livestemp_to_litter, & 
     !$acc veg_pf%m_livecrootp_xfer_to_litter, & 
     !$acc veg_pf%m_deadcrootp_storage_to_litter, & 
     !$acc veg_pf%m_leafp_to_litter, & 
     !$acc veg_pf%m_retransp_to_litter, & 
     !$acc veg_pf%m_deadcrootp_xfer_to_litter, & 
     !$acc veg_pf%m_deadcrootp_to_litter, & 
     !$acc veg_pf%m_leafp_xfer_to_litter, & 
     !$acc veg_pf%m_ppool_to_litter, & 
     !$acc veg_pf%m_frootp_xfer_to_litter, & 
     !$acc veg_pf%hrv_leafp_to_litter, & 
     !$acc veg_pf%hrv_retransp_to_litter, & 
     !$acc veg_pf%hrv_deadcrootp_to_litter, & 
     !$acc veg_pf%hrv_leafp_xfer_to_litter, & 
     !$acc veg_pf%hrv_ppool_to_litter, & 
     !$acc veg_pf%hrv_deadstemp_storage_to_litter, & 
     !$acc veg_pf%hrv_livestemp_storage_to_litter, & 
     !$acc veg_pf%hrv_deadcrootp_xfer_to_litter, & 
     !$acc veg_pf%hrv_frootp_storage_to_litter, & 
     !$acc veg_pf%hrv_deadcrootp_storage_to_litter, & 
     !$acc veg_pf%hrv_deadstemp_to_prod10p, & 
     !$acc veg_pf%hrv_deadstemp_xfer_to_litter, & 
     !$acc veg_pf%hrv_livestemp_xfer_to_litter, & 
     !$acc veg_pf%hrv_livecrootp_to_litter, & 
     !$acc veg_pf%hrv_frootp_to_litter, & 
     !$acc veg_pf%hrv_livecrootp_storage_to_litter, & 
     !$acc veg_pf%hrv_deadstemp_to_prod100p, & 
     !$acc veg_pf%hrv_leafp_storage_to_litter, & 
     !$acc veg_pf%hrv_livestemp_to_litter, & 
     !$acc veg_pf%hrv_frootp_xfer_to_litter, & 
     !$acc veg_pf%hrv_livecrootp_xfer_to_litter, & 
     !$acc veg_pf%m_leafp_to_fire, & 
     !$acc veg_pf%m_deadstemp_xfer_to_litter_fire, & 
     !$acc veg_pf%m_livestemp_storage_to_litter_fire, & 
     !$acc veg_pf%m_livestemp_xfer_to_litter_fire, & 
     !$acc veg_pf%m_frootp_xfer_to_litter_fire, & 
     !$acc veg_pf%m_retransp_to_litter_fire, & 
     !$acc veg_pf%m_deadcrootp_xfer_to_litter_fire, & 
     !$acc veg_pf%m_deadcrootp_xfer_to_fire, & 
     !$acc veg_pf%m_deadcrootp_to_fire, & 
     !$acc veg_pf%m_livecrootp_storage_to_litter_fire, & 
     !$acc veg_pf%m_livecrootp_storage_to_fire, & 
     !$acc veg_pf%m_livestemp_to_deadstemp_fire, & 
     !$acc veg_pf%m_frootp_to_fire, & 
     !$acc veg_pf%m_leafp_to_litter_fire, & 
     !$acc veg_pf%m_livecrootp_to_deadcrootp_fire, & 
     !$acc veg_pf%m_deadstemp_to_fire, & 
     !$acc veg_pf%m_deadstemp_storage_to_fire, & 
     !$acc veg_pf%m_ppool_to_fire, & 
     !$acc veg_pf%m_livestemp_to_fire, & 
     !$acc veg_pf%m_frootp_storage_to_fire, & 
     !$acc veg_pf%m_livecrootp_to_litter_fire, & 
     !$acc veg_pf%m_frootp_storage_to_litter_fire, & 
     !$acc veg_pf%m_livecrootp_xfer_to_fire, & 
     !$acc veg_pf%m_ppool_to_litter_fire, & 
     !$acc veg_pf%m_deadstemp_storage_to_litter_fire, & 
     !$acc veg_pf%m_frootp_xfer_to_fire, & 
     !$acc veg_pf%m_livestemp_storage_to_fire, & 
     !$acc veg_pf%m_livestemp_to_litter_fire, & 
     !$acc veg_pf%m_deadstemp_xfer_to_fire, & 
     !$acc veg_pf%m_retransp_to_fire, & 
     !$acc veg_pf%m_livecrootp_to_fire, & 
     !$acc veg_pf%m_leafp_storage_to_litter_fire, & 
     !$acc veg_pf%m_leafp_xfer_to_fire, & 
     !$acc veg_pf%m_frootp_to_litter_fire, & 
     !$acc veg_pf%m_leafp_storage_to_fire, & 
     !$acc veg_pf%m_deadcrootp_storage_to_fire, & 
     !$acc veg_pf%m_deadcrootp_storage_to_litter_fire, & 
     !$acc veg_pf%m_livecrootp_xfer_to_litter_fire, & 
     !$acc veg_pf%m_deadstemp_to_litter_fire, & 
     !$acc veg_pf%m_deadcrootp_to_litter_fire, & 
     !$acc veg_pf%m_leafp_xfer_to_litter_fire, & 
     !$acc veg_pf%m_livestemp_xfer_to_fire )
     !$acc update self(& 
     !$acc veg_ns%pnup_pfrootc, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn, & 
     !$acc veg_ns%deadstemn, & 
     !$acc veg_ns%grainn_storage, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%retransn, & 
     !$acc veg_ns%npool, & 
     !$acc veg_ns%cropseedn_deficit, & 
     !$acc veg_ns%grainn_xfer, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn_storage, & 
     !$acc veg_ns%frootn, & 
     !$acc veg_ns%livecrootn, & 
     !$acc veg_ns%leafn, & 
     !$acc veg_ns%livecrootn_storage, & 
     !$acc veg_ns%frootn_storage, & 
     !$acc veg_ns%leafn_storage, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%deadstemn_storage, & 
     !$acc veg_ns%deadcrootn_storage, & 
     !$acc veg_ns%deadcrootn, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%grainn, & 
     !$acc veg_ns%livestemn, & 
     !$acc veg_ns%deadcrootn_storage, & 
     !$acc veg_ns%deadstemn, & 
     !$acc veg_ns%npool, & 
     !$acc veg_ns%frootn, & 
     !$acc veg_ns%deadstemn_storage, & 
     !$acc veg_ns%livecrootn, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%livestemn_storage, & 
     !$acc veg_ns%deadcrootn, & 
     !$acc veg_ns%livecrootn_storage, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%frootn_storage, & 
     !$acc veg_ns%leafn_storage, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%retransn, & 
     !$acc veg_ns%leafn, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn, & 
     !$acc veg_ns%deadstemn, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%retransn, & 
     !$acc veg_ns%npool, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn_storage, & 
     !$acc veg_ns%frootn, & 
     !$acc veg_ns%leafn, & 
     !$acc veg_ns%livecrootn_storage, & 
     !$acc veg_ns%frootn_storage, & 
     !$acc veg_ns%leafn_storage, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%deadstemn_storage, & 
     !$acc veg_ns%grainn, & 
     !$acc veg_ns%deadcrootn_storage, & 
     !$acc veg_ns%deadcrootn, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%livecrootn )
     !$acc update self(& 
     !$acc col_ns%smin_no3_vr, & 
     !$acc col_ns%smin_nh4_vr, & 
     !$acc col_ns%smin_nh4_vr, & 
     !$acc col_ns%smin_no3_vr, & 
     !$acc col_ns%sminn_vr, & 
     !$acc col_ns%decomp_npools_vr, & 
     !$acc col_ns%decomp_npools_vr, & 
     !$acc col_ns%decomp_npools_vr, & 
     !$acc col_ns%prod100n, & 
     !$acc col_ns%prod10n, & 
     !$acc col_ns%prod1n )
     !$acc update self(& 
     !$acc veg_pp%itype, & 
     !$acc veg_pp%itype )
     !$acc update self(& 
     !$acc veg_vp%woody )
     !$acc update self(& 
     !$acc veg_cf%cpool_to_leafc_storage, & 
     !$acc veg_cf%xsmrpool_turnover, & 
     !$acc veg_cf%grain_curmr, & 
     !$acc veg_cf%froot_mr, & 
     !$acc veg_cf%cpool_to_deadcrootc_storage, & 
     !$acc veg_cf%livestem_curmr, & 
     !$acc veg_cf%grain_mr, & 
     !$acc veg_cf%livestem_xsmr, & 
     !$acc veg_cf%allocation_stem, & 
     !$acc veg_cf%cpool_to_grainc_storage, & 
     !$acc veg_cf%froot_curmr, & 
     !$acc veg_cf%availc, & 
     !$acc veg_cf%livecroot_curmr, & 
     !$acc veg_cf%livecroot_mr, & 
     !$acc veg_cf%plant_calloc, & 
     !$acc veg_cf%cpool_to_grainc, & 
     !$acc veg_cf%livestem_mr, & 
     !$acc veg_cf%leaf_curmr, & 
     !$acc veg_cf%cpool_to_gresp_storage, & 
     !$acc veg_cf%cpool_to_livestemc_storage, & 
     !$acc veg_cf%grain_xsmr, & 
     !$acc veg_cf%cpool_to_deadstemc_storage, & 
     !$acc veg_cf%froot_xsmr, & 
     !$acc veg_cf%cpool_to_frootc_storage, & 
     !$acc veg_cf%excess_cflux, & 
     !$acc veg_cf%psnshade_to_cpool, & 
     !$acc veg_cf%psnsun_to_cpool, & 
     !$acc veg_cf%cpool_to_leafc, & 
     !$acc veg_cf%cpool_to_livecrootc, & 
     !$acc veg_cf%livecroot_xsmr, & 
     !$acc veg_cf%allocation_froot, & 
     !$acc veg_cf%leaf_mr, & 
     !$acc veg_cf%cpool_to_frootc, & 
     !$acc veg_cf%xsmrpool_recover, & 
     !$acc veg_cf%allocation_leaf, & 
     !$acc veg_cf%cpool_to_deadstemc, & 
     !$acc veg_cf%cpool_to_xsmrpool, & 
     !$acc veg_cf%cpool_to_deadcrootc, & 
     !$acc veg_cf%leaf_xsmr, & 
     !$acc veg_cf%cpool_to_livecrootc_storage, & 
     !$acc veg_cf%cpool_to_livestemc, & 
     !$acc veg_cf%prev_frootc_to_litter, & 
     !$acc veg_cf%frootc_storage_to_xfer, & 
     !$acc veg_cf%leafc_xfer_to_leafc, & 
     !$acc veg_cf%livecrootc_storage_to_xfer, & 
     !$acc veg_cf%livestemc_xfer_to_livestemc, & 
     !$acc veg_cf%livestemc_storage_to_xfer, & 
     !$acc veg_cf%deadstemc_xfer_to_deadstemc, & 
     !$acc veg_cf%deadstemc_storage_to_xfer, & 
     !$acc veg_cf%leafc_storage_to_xfer, & 
     !$acc veg_cf%gresp_storage_to_xfer, & 
     !$acc veg_cf%prev_leafc_to_litter, & 
     !$acc veg_cf%deadcrootc_xfer_to_deadcrootc, & 
     !$acc veg_cf%livecrootc_xfer_to_livecrootc, & 
     !$acc veg_cf%frootc_xfer_to_frootc, & 
     !$acc veg_cf%deadcrootc_storage_to_xfer, & 
     !$acc veg_cf%prev_frootc_to_litter, & 
     !$acc veg_cf%frootc_storage_to_xfer, & 
     !$acc veg_cf%leafc_xfer_to_leafc, & 
     !$acc veg_cf%livecrootc_storage_to_xfer, & 
     !$acc veg_cf%livestemc_xfer_to_livestemc, & 
     !$acc veg_cf%livestemc_storage_to_xfer, & 
     !$acc veg_cf%deadstemc_xfer_to_deadstemc, & 
     !$acc veg_cf%deadstemc_storage_to_xfer, & 
     !$acc veg_cf%leafc_storage_to_xfer, & 
     !$acc veg_cf%gresp_storage_to_xfer, & 
     !$acc veg_cf%prev_leafc_to_litter, & 
     !$acc veg_cf%deadcrootc_xfer_to_deadcrootc, & 
     !$acc veg_cf%livecrootc_xfer_to_livecrootc, & 
     !$acc veg_cf%frootc_xfer_to_frootc, & 
     !$acc veg_cf%deadcrootc_storage_to_xfer, & 
     !$acc veg_cf%crop_seedc_to_leaf, & 
     !$acc veg_cf%leafc_xfer_to_leafc, & 
     !$acc veg_cf%livestemc_xfer_to_livestemc, & 
     !$acc veg_cf%deadstemc_xfer_to_deadstemc, & 
     !$acc veg_cf%deadcrootc_xfer_to_deadcrootc, & 
     !$acc veg_cf%livecrootc_xfer_to_livecrootc, & 
     !$acc veg_cf%frootc_xfer_to_frootc, & 
     !$acc veg_cf%hrv_livestemc_to_prod1c, & 
     !$acc veg_cf%hrv_leafc_to_prod1c, & 
     !$acc veg_cf%hrv_grainc_to_prod1c, & 
     !$acc veg_cf%hrv_cropc_to_prod1c, & 
     !$acc veg_cf%prev_frootc_to_litter, & 
     !$acc veg_cf%frootc_to_litter, & 
     !$acc veg_cf%leafc_to_litter, & 
     !$acc veg_cf%livestemc_to_litter, & 
     !$acc veg_cf%prev_leafc_to_litter, & 
     !$acc veg_cf%frootc_to_litter, & 
     !$acc veg_cf%leafc_to_litter, & 
     !$acc veg_cf%livestemc_to_deadstemc, & 
     !$acc veg_cf%livecrootc_to_deadcrootc, & 
     !$acc veg_cf%cpool_livecroot_storage_gr, & 
     !$acc veg_cf%transfer_deadstem_gr, & 
     !$acc veg_cf%transfer_livecroot_gr, & 
     !$acc veg_cf%cpool_grain_storage_gr, & 
     !$acc veg_cf%cpool_leaf_storage_gr, & 
     !$acc veg_cf%transfer_livestem_gr, & 
     !$acc veg_cf%transfer_grain_gr, & 
     !$acc veg_cf%cpool_froot_storage_gr, & 
     !$acc veg_cf%cpool_deadstem_storage_gr, & 
     !$acc veg_cf%cpool_livestem_gr, & 
     !$acc veg_cf%cpool_livestem_storage_gr, & 
     !$acc veg_cf%cpool_livecroot_gr, & 
     !$acc veg_cf%cpool_deadcroot_gr, & 
     !$acc veg_cf%cpool_froot_gr, & 
     !$acc veg_cf%cpool_grain_gr, & 
     !$acc veg_cf%transfer_froot_gr, & 
     !$acc veg_cf%transfer_leaf_gr, & 
     !$acc veg_cf%cpool_leaf_gr, & 
     !$acc veg_cf%cpool_deadcroot_storage_gr, & 
     !$acc veg_cf%transfer_deadcroot_gr, & 
     !$acc veg_cf%cpool_deadstem_gr, & 
     !$acc veg_cf%cpool_livecroot_storage_gr, & 
     !$acc veg_cf%grainc_xfer_to_grainc, & 
     !$acc veg_cf%cpool_to_leafc_storage, & 
     !$acc veg_cf%livestemc_xfer_to_livestemc, & 
     !$acc veg_cf%grain_curmr, & 
     !$acc veg_cf%deadstemc_xfer_to_deadstemc, & 
     !$acc veg_cf%cpool_to_deadcrootc_storage, & 
     !$acc veg_cf%cpool_grain_storage_gr, & 
     !$acc veg_cf%cpool_leaf_storage_gr, & 
     !$acc veg_cf%livestem_curmr, & 
     !$acc veg_cf%livestemc_to_deadstemc, & 
     !$acc veg_cf%deadcrootc_storage_to_xfer, & 
     !$acc veg_cf%cpool_to_grainc_storage, & 
     !$acc veg_cf%froot_curmr, & 
     !$acc veg_cf%cpool_froot_storage_gr, & 
     !$acc veg_cf%cpool_deadstem_storage_gr, & 
     !$acc veg_cf%cpool_livestem_gr, & 
     !$acc veg_cf%livecroot_curmr, & 
     !$acc veg_cf%cpool_livestem_storage_gr, & 
     !$acc veg_cf%livestemc_storage_to_xfer, & 
     !$acc veg_cf%cpool_deadcroot_gr, & 
     !$acc veg_cf%cpool_livecroot_gr, & 
     !$acc veg_cf%cpool_to_grainc, & 
     !$acc veg_cf%gresp_storage_to_xfer, & 
     !$acc veg_cf%leaf_curmr, & 
     !$acc veg_cf%grainc_storage_to_xfer, & 
     !$acc veg_cf%livecrootc_to_deadcrootc, & 
     !$acc veg_cf%deadcrootc_xfer_to_deadcrootc, & 
     !$acc veg_cf%cpool_to_gresp_storage, & 
     !$acc veg_cf%cpool_to_livestemc_storage, & 
     !$acc veg_cf%frootc_xfer_to_frootc, & 
     !$acc veg_cf%frootc_storage_to_xfer, & 
     !$acc veg_cf%leafc_xfer_to_leafc, & 
     !$acc veg_cf%cpool_to_deadstemc_storage, & 
     !$acc veg_cf%cpool_to_frootc_storage, & 
     !$acc veg_cf%cpool_froot_gr, & 
     !$acc veg_cf%grainc_to_food, & 
     !$acc veg_cf%frootc_to_litter, & 
     !$acc veg_cf%cpool_to_leafc, & 
     !$acc veg_cf%cpool_to_livecrootc, & 
     !$acc veg_cf%cpool_grain_gr, & 
     !$acc veg_cf%livecrootc_xfer_to_livecrootc, & 
     !$acc veg_cf%cpool_to_frootc, & 
     !$acc veg_cf%livecrootc_storage_to_xfer, & 
     !$acc veg_cf%cpool_leaf_gr, & 
     !$acc veg_cf%cpool_deadcroot_storage_gr, & 
     !$acc veg_cf%deadstemc_storage_to_xfer, & 
     !$acc veg_cf%leafc_storage_to_xfer, & 
     !$acc veg_cf%cpool_to_xsmrpool, & 
     !$acc veg_cf%leafc_to_litter, & 
     !$acc veg_cf%cpool_to_deadcrootc, & 
     !$acc veg_cf%cpool_to_deadstemc, & 
     !$acc veg_cf%livestemc_to_litter, & 
     !$acc veg_cf%cpool_to_livecrootc_storage, & 
     !$acc veg_cf%cpool_deadstem_gr, & 
     !$acc veg_cf%cpool_to_livestemc, & 
     !$acc veg_cf%xsmrpool_to_atm, & 
     !$acc veg_cf%m_livecrootc_storage_to_litter, & 
     !$acc veg_cf%m_gresp_xfer_to_litter, & 
     !$acc veg_cf%m_deadcrootc_to_litter, & 
     !$acc veg_cf%m_deadcrootc_xfer_to_litter, & 
     !$acc veg_cf%m_livecrootc_to_litter, & 
     !$acc veg_cf%m_deadstemc_storage_to_litter, & 
     !$acc veg_cf%m_livestemc_to_litter, & 
     !$acc veg_cf%m_gresp_storage_to_litter, & 
     !$acc veg_cf%m_leafc_to_litter, & 
     !$acc veg_cf%m_frootc_storage_to_litter, & 
     !$acc veg_cf%m_cpool_to_litter, & 
     !$acc veg_cf%m_livestemc_storage_to_litter, & 
     !$acc veg_cf%m_frootc_xfer_to_litter, & 
     !$acc veg_cf%m_deadstemc_to_litter, & 
     !$acc veg_cf%m_leafc_xfer_to_litter, & 
     !$acc veg_cf%m_frootc_to_litter, & 
     !$acc veg_cf%m_deadstemc_xfer_to_litter, & 
     !$acc veg_cf%m_leafc_storage_to_litter, & 
     !$acc veg_cf%m_livestemc_xfer_to_litter, & 
     !$acc veg_cf%m_deadcrootc_storage_to_litter, & 
     !$acc veg_cf%m_livecrootc_xfer_to_litter, & 
     !$acc veg_cf%hrv_deadstemc_storage_to_litter, & 
     !$acc veg_cf%hrv_cpool_to_litter, & 
     !$acc veg_cf%hrv_livestemc_storage_to_litter, & 
     !$acc veg_cf%hrv_leafc_storage_to_litter, & 
     !$acc veg_cf%hrv_deadcrootc_to_litter, & 
     !$acc veg_cf%hrv_deadstemc_to_prod10c, & 
     !$acc veg_cf%hrv_deadstemc_to_prod100c, & 
     !$acc veg_cf%hrv_leafc_xfer_to_litter, & 
     !$acc veg_cf%hrv_leafc_to_litter, & 
     !$acc veg_cf%hrv_livecrootc_xfer_to_litter, & 
     !$acc veg_cf%hrv_frootc_to_litter, & 
     !$acc veg_cf%hrv_deadcrootc_xfer_to_litter, & 
     !$acc veg_cf%hrv_xsmrpool_to_atm, & 
     !$acc veg_cf%hrv_livestemc_to_litter, & 
     !$acc veg_cf%hrv_deadcrootc_storage_to_litter, & 
     !$acc veg_cf%hrv_deadstemc_xfer_to_litter, & 
     !$acc veg_cf%hrv_gresp_storage_to_litter, & 
     !$acc veg_cf%hrv_livecrootc_storage_to_litter, & 
     !$acc veg_cf%hrv_livestemc_xfer_to_litter, & 
     !$acc veg_cf%hrv_frootc_xfer_to_litter, & 
     !$acc veg_cf%hrv_livecrootc_to_litter, & 
     !$acc veg_cf%hrv_frootc_storage_to_litter, & 
     !$acc veg_cf%hrv_gresp_xfer_to_litter, & 
     !$acc veg_cf%m_livecrootc_xfer_to_fire, & 
     !$acc veg_cf%m_leafc_xfer_to_fire, & 
     !$acc veg_cf%m_livestemc_to_fire, & 
     !$acc veg_cf%m_livestemc_xfer_to_litter_fire, & 
     !$acc veg_cf%m_frootc_storage_to_fire, & 
     !$acc veg_cf%m_frootc_to_litter_fire, & 
     !$acc veg_cf%m_deadcrootc_to_litter_fire, & 
     !$acc veg_cf%m_deadstemc_storage_to_litter_fire, & 
     !$acc veg_cf%m_livestemc_storage_to_litter_fire, & 
     !$acc veg_cf%m_livestemc_to_litter_fire, & 
     !$acc veg_cf%m_deadstemc_xfer_to_fire, & 
     !$acc veg_cf%m_gresp_xfer_to_litter_fire, & 
     !$acc veg_cf%m_deadstemc_to_litter_fire, & 
     !$acc veg_cf%m_leafc_to_fire, & 
     !$acc veg_cf%m_deadstemc_storage_to_fire, & 
     !$acc veg_cf%m_deadcrootc_storage_to_fire, & 
     !$acc veg_cf%m_deadcrootc_xfer_to_litter_fire, & 
     !$acc veg_cf%m_livestemc_to_deadstemc_fire, & 
     !$acc veg_cf%m_frootc_storage_to_litter_fire, & 
     !$acc veg_cf%m_frootc_to_fire, & 
     !$acc veg_cf%m_deadcrootc_xfer_to_fire, & 
     !$acc veg_cf%m_gresp_storage_to_litter_fire, & 
     !$acc veg_cf%m_livestemc_xfer_to_fire, & 
     !$acc veg_cf%m_frootc_xfer_to_fire, & 
     !$acc veg_cf%m_cpool_to_litter_fire, & 
     !$acc veg_cf%m_livecrootc_storage_to_litter_fire, & 
     !$acc veg_cf%m_deadcrootc_storage_to_litter_fire, & 
     !$acc veg_cf%m_deadstemc_xfer_to_litter_fire, & 
     !$acc veg_cf%m_deadcrootc_to_fire, & 
     !$acc veg_cf%m_livecrootc_storage_to_fire, & 
     !$acc veg_cf%m_gresp_xfer_to_fire, & 
     !$acc veg_cf%m_leafc_to_litter_fire, & 
     !$acc veg_cf%m_frootc_xfer_to_litter_fire, & 
     !$acc veg_cf%m_livestemc_storage_to_fire, & 
     !$acc veg_cf%m_livecrootc_xfer_to_litter_fire, & 
     !$acc veg_cf%m_livecrootc_to_litter_fire, & 
     !$acc veg_cf%m_leafc_xfer_to_litter_fire, & 
     !$acc veg_cf%m_deadstemc_to_fire, & 
     !$acc veg_cf%m_livecrootc_to_fire, & 
     !$acc veg_cf%m_gresp_storage_to_fire, & 
     !$acc veg_cf%m_leafc_storage_to_litter_fire, & 
     !$acc veg_cf%m_leafc_storage_to_fire, & 
     !$acc veg_cf%m_cpool_to_fire, & 
     !$acc veg_cf%m_livecrootc_to_deadcrootc_fire )
     !$acc update self(& 
     !$acc canopystate_vars%laisun_patch, & 
     !$acc canopystate_vars%laisha_patch )
     !$acc update self(& 
     !$acc c13_veg_cf%psnshade_to_cpool, & 
     !$acc c13_veg_cf%psnsun_to_cpool )
     !$acc update self(& 
     !$acc c14_veg_cf%psnshade_to_cpool, & 
     !$acc c14_veg_cf%psnsun_to_cpool )
     !$acc update self(& 
     !$acc veg_es%gdd820, & 
     !$acc veg_es%gdd1020, & 
     !$acc veg_es%gdd020, & 
     !$acc veg_es%t_ref2m )
     !$acc update self(& 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer, & 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%cpool, & 
     !$acc veg_cs%cpool, & 
     !$acc veg_cs%gresp_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%grainc, & 
     !$acc veg_cs%leafc_storage, & 
     !$acc veg_cs%deadstemc, & 
     !$acc veg_cs%deadcrootc, & 
     !$acc veg_cs%livecrootc_storage, & 
     !$acc veg_cs%livestemc_storage, & 
     !$acc veg_cs%frootc_storage, & 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%cropseedc_deficit, & 
     !$acc veg_cs%xsmrpool, & 
     !$acc veg_cs%grainc_storage, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%gresp_storage, & 
     !$acc veg_cs%grainc_xfer, & 
     !$acc veg_cs%frootc, & 
     !$acc veg_cs%livestemc, & 
     !$acc veg_cs%leafc, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%deadstemc_storage, & 
     !$acc veg_cs%livecrootc, & 
     !$acc veg_cs%deadcrootc_storage, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer, & 
     !$acc veg_cs%cpool, & 
     !$acc veg_cs%gresp_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%leafc_storage, & 
     !$acc veg_cs%deadstemc, & 
     !$acc veg_cs%deadcrootc, & 
     !$acc veg_cs%livecrootc_storage, & 
     !$acc veg_cs%livestemc_storage, & 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%frootc_storage, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%gresp_storage, & 
     !$acc veg_cs%frootc, & 
     !$acc veg_cs%livestemc, & 
     !$acc veg_cs%leafc, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%deadstemc_storage, & 
     !$acc veg_cs%livecrootc, & 
     !$acc veg_cs%deadcrootc_storage, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer, & 
     !$acc veg_cs%cpool, & 
     !$acc veg_cs%gresp_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%grainc, & 
     !$acc veg_cs%leafc_storage, & 
     !$acc veg_cs%deadstemc, & 
     !$acc veg_cs%deadcrootc, & 
     !$acc veg_cs%livecrootc_storage, & 
     !$acc veg_cs%livestemc_storage, & 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%frootc_storage, & 
     !$acc veg_cs%xsmrpool, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%gresp_storage, & 
     !$acc veg_cs%frootc, & 
     !$acc veg_cs%livestemc, & 
     !$acc veg_cs%leafc, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%deadstemc_storage, & 
     !$acc veg_cs%livecrootc, & 
     !$acc veg_cs%deadcrootc_storage, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer, & 
     !$acc veg_cs%totvegc, & 
     !$acc veg_cs%deadstemc, & 
     !$acc veg_cs%leafc, & 
     !$acc veg_cs%cpool, & 
     !$acc veg_cs%gresp_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%leafc_storage, & 
     !$acc veg_cs%deadstemc, & 
     !$acc veg_cs%deadcrootc, & 
     !$acc veg_cs%livecrootc_storage, & 
     !$acc veg_cs%livestemc_storage, & 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%frootc_storage, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%gresp_storage, & 
     !$acc veg_cs%frootc, & 
     !$acc veg_cs%livestemc, & 
     !$acc veg_cs%leafc, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%deadstemc_storage, & 
     !$acc veg_cs%livecrootc, & 
     !$acc veg_cs%deadcrootc_storage, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer )
     !$acc update self(& 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%deadcrootp_storage, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%livecrootp, & 
     !$acc veg_ps%cropseedp_deficit, & 
     !$acc veg_ps%retransp, & 
     !$acc veg_ps%leafp_storage, & 
     !$acc veg_ps%deadstemp_storage, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%ppool, & 
     !$acc veg_ps%livecrootp_storage, & 
     !$acc veg_ps%leafp, & 
     !$acc veg_ps%grainp, & 
     !$acc veg_ps%livestemp_storage, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%grainp_storage, & 
     !$acc veg_ps%livestemp, & 
     !$acc veg_ps%grainp_xfer, & 
     !$acc veg_ps%deadcrootp, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%frootp_storage, & 
     !$acc veg_ps%frootp, & 
     !$acc veg_ps%deadstemp, & 
     !$acc veg_ps%leafp, & 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%deadcrootp_storage, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%livecrootp, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%retransp, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%deadstemp_storage, & 
     !$acc veg_ps%livestemp, & 
     !$acc veg_ps%leafp_storage, & 
     !$acc veg_ps%deadstemp, & 
     !$acc veg_ps%deadcrootp, & 
     !$acc veg_ps%frootp_storage, & 
     !$acc veg_ps%frootp, & 
     !$acc veg_ps%livecrootp_storage, & 
     !$acc veg_ps%livestemp_storage, & 
     !$acc veg_ps%ppool, & 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%deadcrootp_storage, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%livecrootp, & 
     !$acc veg_ps%retransp, & 
     !$acc veg_ps%leafp_storage, & 
     !$acc veg_ps%deadstemp_storage, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%livecrootp_storage, & 
     !$acc veg_ps%ppool, & 
     !$acc veg_ps%leafp, & 
     !$acc veg_ps%grainp, & 
     !$acc veg_ps%livestemp_storage, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%livestemp, & 
     !$acc veg_ps%deadcrootp, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%frootp_storage, & 
     !$acc veg_ps%frootp, & 
     !$acc veg_ps%deadstemp )
     !$acc update self(& 
     !$acc top_as%windbot, & 
     !$acc top_as%pbot, & 
     !$acc top_as%rhbot )
     !$acc update self(& 
     !$acc veg_ef%eflx_soil_grnd, & 
     !$acc veg_ef%netrad )
     !$acc update self(& 
     !$acc crop_vars%prev_xt_bar_patch, & 
     !$acc crop_vars%xt_patch, & 
     !$acc crop_vars%p2ETo_patch, & 
     !$acc crop_vars%cvp_patch, & 
     !$acc crop_vars%ETo_patch, & 
     !$acc crop_vars%xp_patch, & 
     !$acc crop_vars%cvt_patch, & 
     !$acc crop_vars%plantmonth_patch, & 
     !$acc crop_vars%xt_bar_patch, & 
     !$acc crop_vars%p2ETo_bar_patch, & 
     !$acc crop_vars%prev_p2ETo_bar_patch, & 
     !$acc crop_vars%prev_xp_bar_patch, & 
     !$acc crop_vars%P2E_rm_patch, & 
     !$acc crop_vars%xp_bar_patch, & 
     !$acc crop_vars%plantday_patch, & 
     !$acc crop_vars%vf_patch, & 
     !$acc crop_vars%crpyld_patch, & 
     !$acc crop_vars%cropplant_patch, & 
     !$acc crop_vars%croplive_patch, & 
     !$acc crop_vars%harvdate_patch, & 
     !$acc crop_vars%dmyield_patch, & 
     !$acc crop_vars%harvday_patch, & 
     !$acc crop_vars%vf_patch, & 
     !$acc crop_vars%crpyld_patch, & 
     !$acc crop_vars%dmyield_patch )
     !$acc update self(& 
     !$acc soilstate_vars%rootfr_patch, & 
     !$acc soilstate_vars%root_depth_patch )
     !$acc update self(& 
     !$acc col_pp%nlevbed )
     !$acc update self(& 
     !$acc decomp_cascade_con%cascade_donor_pool, & 
     !$acc decomp_cascade_con%cascade_receiver_pool, & 
     !$acc decomp_cascade_con%cascade_donor_pool, & 
     !$acc decomp_cascade_con%cascade_receiver_pool, & 
     !$acc decomp_cascade_con%cascade_donor_pool, & 
     !$acc decomp_cascade_con%cascade_receiver_pool )
     !$acc update self(& 
     !$acc col_cs%decomp_som2c_vr, & 
     !$acc col_cs%decomp_cpools_vr, & 
     !$acc col_cs%decomp_cpools_vr, & 
     !$acc col_cs%decomp_cpools_vr, & 
     !$acc col_cs%prod100c, & 
     !$acc col_cs%prod10c, & 
     !$acc col_cs%prod1c, & 
     !$acc col_cs%leafc, & 
     !$acc col_cs%fuelc, & 
     !$acc col_cs%deadstemc, & 
     !$acc col_cs%totvegc, & 
     !$acc col_cs%rootc, & 
     !$acc col_cs%fuelc_crop, & 
     !$acc col_cs%decomp_cpools_vr )
     !$acc update self(& 
     !$acc col_ps%solutionp_vr, & 
     !$acc col_ps%decomp_ppools_vr, & 
     !$acc col_ps%decomp_ppools_vr, & 
     !$acc col_ps%decomp_ppools_vr, & 
     !$acc col_ps%prod100p, & 
     !$acc col_ps%prod10p, & 
     !$acc col_ps%prod1p )
     !$acc update self(& 
     !$acc c13_col_cs%decomp_cpools_vr, & 
     !$acc c13_col_cs%prod100c, & 
     !$acc c13_col_cs%prod10c, & 
     !$acc c13_col_cs%prod1c )
     !$acc update self(& 
     !$acc c14_col_cs%decomp_cpools_vr, & 
     !$acc c14_col_cs%prod100c, & 
     !$acc c14_col_cs%prod10c, & 
     !$acc c14_col_cs%prod1c, & 
     !$acc c14_col_cs%seedc, & 
     !$acc c14_col_cs%decomp_cpools_vr )
     !$acc update self(& 
     !$acc c13_col_cf%decomp_cpools_transport_tendency, & 
     !$acc c13_col_cf%prod10c_loss, & 
     !$acc c13_col_cf%prod100c_loss, & 
     !$acc c13_col_cf%prod1c_loss )
     !$acc update self(& 
     !$acc c14_col_cf%decomp_cpools_transport_tendency, & 
     !$acc c14_col_cf%prod10c_loss, & 
     !$acc c14_col_cf%prod100c_loss, & 
     !$acc c14_col_cf%prod1c_loss )
     !$acc update self(& 
     !$acc c14_veg_cs%cpool, & 
     !$acc c14_veg_cs%gresp_xfer, & 
     !$acc c14_veg_cs%deadstemc_xfer, & 
     !$acc c14_veg_cs%leafc_storage, & 
     !$acc c14_veg_cs%deadstemc, & 
     !$acc c14_veg_cs%deadcrootc, & 
     !$acc c14_veg_cs%livecrootc_storage, & 
     !$acc c14_veg_cs%livestemc_storage, & 
     !$acc c14_veg_cs%frootc_storage, & 
     !$acc c14_veg_cs%livestemc_xfer, & 
     !$acc c14_veg_cs%xsmrpool, & 
     !$acc c14_veg_cs%leafc_xfer, & 
     !$acc c14_veg_cs%gresp_storage, & 
     !$acc c14_veg_cs%ctrunc, & 
     !$acc c14_veg_cs%frootc, & 
     !$acc c14_veg_cs%livestemc, & 
     !$acc c14_veg_cs%leafc, & 
     !$acc c14_veg_cs%livecrootc_xfer, & 
     !$acc c14_veg_cs%deadstemc_storage, & 
     !$acc c14_veg_cs%livecrootc, & 
     !$acc c14_veg_cs%deadcrootc_storage, & 
     !$acc c14_veg_cs%frootc_xfer, & 
     !$acc c14_veg_cs%deadcrootc_xfer )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_nf%pmnf_decomp_cascade' 
     write(10,*) col_nf%pmnf_decomp_cascade
     write(10,*) 'col_nf%soil_n_immob_flux_vr' 
     write(10,*) col_nf%soil_n_immob_flux_vr
     write(10,*) 'col_nf%potential_immob_vr' 
     write(10,*) col_nf%potential_immob_vr
     write(10,*) 'col_nf%net_nmin_vr' 
     write(10,*) col_nf%net_nmin_vr
     write(10,*) 'col_nf%soil_n_immob_flux' 
     write(10,*) col_nf%soil_n_immob_flux
     write(10,*) 'col_nf%gross_nmin_vr' 
     write(10,*) col_nf%gross_nmin_vr
     write(10,*) 'col_nf%soil_n_grossmin_flux' 
     write(10,*) col_nf%soil_n_grossmin_flux
     write(10,*) 'col_nf%decomp_cascade_ntransfer_vr' 
     write(10,*) col_nf%decomp_cascade_ntransfer_vr
     write(10,*) 'col_nf%decomp_cascade_sminn_flux_vr' 
     write(10,*) col_nf%decomp_cascade_sminn_flux_vr
     write(10,*) 'col_nf%sminn_to_denit_decomp_cascade_vr' 
     write(10,*) col_nf%sminn_to_denit_decomp_cascade_vr
     write(10,*) 'col_nf%fmax_denit_nitrate_vr' 
     write(10,*) col_nf%fmax_denit_nitrate_vr
     write(10,*) 'col_nf%r_psi' 
     write(10,*) col_nf%r_psi
     write(10,*) 'col_nf%k_nitr_t_vr' 
     write(10,*) col_nf%k_nitr_t_vr
     write(10,*) 'col_nf%diffus' 
     write(10,*) col_nf%diffus
     write(10,*) 'col_nf%fr_WFPS' 
     write(10,*) col_nf%fr_WFPS
     write(10,*) 'col_nf%fmax_denit_carbonsubstrate_vr' 
     write(10,*) col_nf%fmax_denit_carbonsubstrate_vr
     write(10,*) 'col_nf%ratio_no3_co2' 
     write(10,*) col_nf%ratio_no3_co2
     write(10,*) 'col_nf%pot_f_denit_vr' 
     write(10,*) col_nf%pot_f_denit_vr
     write(10,*) 'col_nf%pot_f_nit_vr' 
     write(10,*) col_nf%pot_f_nit_vr
     write(10,*) 'col_nf%ratio_k1' 
     write(10,*) col_nf%ratio_k1
     write(10,*) 'col_nf%n2_n2o_ratio_denit_vr' 
     write(10,*) col_nf%n2_n2o_ratio_denit_vr
     write(10,*) 'col_nf%f_denit_base_vr' 
     write(10,*) col_nf%f_denit_base_vr
     write(10,*) 'col_nf%soil_bulkdensity' 
     write(10,*) col_nf%soil_bulkdensity
     write(10,*) 'col_nf%k_nitr_h2o_vr' 
     write(10,*) col_nf%k_nitr_h2o_vr
     write(10,*) 'col_nf%k_nitr_ph_vr' 
     write(10,*) col_nf%k_nitr_ph_vr
     write(10,*) 'col_nf%anaerobic_frac' 
     write(10,*) col_nf%anaerobic_frac
     write(10,*) 'col_nf%soil_co2_prod' 
     write(10,*) col_nf%soil_co2_prod
     write(10,*) 'col_nf%k_nitr_vr' 
     write(10,*) col_nf%k_nitr_vr
     write(10,*) 'col_nf%wfps_vr' 
     write(10,*) col_nf%wfps_vr
     write(10,*) 'col_nf%smin_no3_massdens_vr' 
     write(10,*) col_nf%smin_no3_massdens_vr
     write(10,*) 'col_nf%actual_immob' 
     write(10,*) col_nf%actual_immob
     write(10,*) 'col_nf%actual_immob_nh4_vr' 
     write(10,*) col_nf%actual_immob_nh4_vr
     write(10,*) 'col_nf%actual_immob_no3_vr' 
     write(10,*) col_nf%actual_immob_no3_vr
     write(10,*) 'col_nf%f_n2o_nit_vr' 
     write(10,*) col_nf%f_n2o_nit_vr
     write(10,*) 'col_nf%supplement_to_sminn_vr' 
     write(10,*) col_nf%supplement_to_sminn_vr
     write(10,*) 'col_nf%actual_immob_vr' 
     write(10,*) col_nf%actual_immob_vr
     write(10,*) 'col_nf%col_plant_no3demand_vr' 
     write(10,*) col_nf%col_plant_no3demand_vr
     write(10,*) 'col_nf%actual_immob_no3' 
     write(10,*) col_nf%actual_immob_no3
     write(10,*) 'col_nf%actual_immob_nh4' 
     write(10,*) col_nf%actual_immob_nh4
     write(10,*) 'col_nf%sminn_to_plant' 
     write(10,*) col_nf%sminn_to_plant
     write(10,*) 'col_nf%sminn_to_denit_excess_vr' 
     write(10,*) col_nf%sminn_to_denit_excess_vr
     write(10,*) 'col_nf%smin_no3_to_plant_vr' 
     write(10,*) col_nf%smin_no3_to_plant_vr
     write(10,*) 'col_nf%f_denit_vr' 
     write(10,*) col_nf%f_denit_vr
     write(10,*) 'col_nf%f_nit_vr' 
     write(10,*) col_nf%f_nit_vr
     write(10,*) 'col_nf%col_plant_nh4demand_vr' 
     write(10,*) col_nf%col_plant_nh4demand_vr
     write(10,*) 'col_nf%f_n2o_denit_vr' 
     write(10,*) col_nf%f_n2o_denit_vr
     write(10,*) 'col_nf%sminn_to_plant_vr' 
     write(10,*) col_nf%sminn_to_plant_vr
     write(10,*) 'col_nf%potential_immob' 
     write(10,*) col_nf%potential_immob
     write(10,*) 'col_nf%col_plant_ndemand_vr' 
     write(10,*) col_nf%col_plant_ndemand_vr
     write(10,*) 'col_nf%smin_nh4_to_plant_vr' 
     write(10,*) col_nf%smin_nh4_to_plant_vr
     write(10,*) 'col_nf%actual_immob' 
     write(10,*) col_nf%actual_immob
     write(10,*) 'col_nf%gross_nmin' 
     write(10,*) col_nf%gross_nmin
     write(10,*) 'col_nf%sminn_to_plant' 
     write(10,*) col_nf%sminn_to_plant
     write(10,*) 'col_nf%net_nmin' 
     write(10,*) col_nf%net_nmin
     write(10,*) 'col_nf%potential_immob' 
     write(10,*) col_nf%potential_immob
     write(10,*) 'col_nf%sminn_to_plant' 
     write(10,*) col_nf%sminn_to_plant
     write(10,*) 'col_nf%smin_no3_to_plant_vr' 
     write(10,*) col_nf%smin_no3_to_plant_vr
     write(10,*) 'col_nf%plant_n_uptake_flux' 
     write(10,*) col_nf%plant_n_uptake_flux
     write(10,*) 'col_nf%sminn_to_plant_vr' 
     write(10,*) col_nf%sminn_to_plant_vr
     write(10,*) 'col_nf%smin_nh4_to_plant_vr' 
     write(10,*) col_nf%smin_nh4_to_plant_vr
     write(10,*) 'col_nf%hrv_cropn_to_prod1n' 
     write(10,*) col_nf%hrv_cropn_to_prod1n
     write(10,*) 'col_nf%phenology_n_to_litr_cel_n' 
     write(10,*) col_nf%phenology_n_to_litr_cel_n
     write(10,*) 'col_nf%phenology_n_to_litr_lig_n' 
     write(10,*) col_nf%phenology_n_to_litr_lig_n
     write(10,*) 'col_nf%phenology_n_to_litr_met_n' 
     write(10,*) col_nf%phenology_n_to_litr_met_n
     write(10,*) 'col_nf%decomp_npools_sourcesink' 
     write(10,*) col_nf%decomp_npools_sourcesink
     write(10,*) 'col_nf%decomp_npools_transport_tendency' 
     write(10,*) col_nf%decomp_npools_transport_tendency
     write(10,*) 'col_nf%gap_mortality_n_to_litr_met_n' 
     write(10,*) col_nf%gap_mortality_n_to_litr_met_n
     write(10,*) 'col_nf%gap_mortality_n_to_litr_lig_n' 
     write(10,*) col_nf%gap_mortality_n_to_litr_lig_n
     write(10,*) 'col_nf%gap_mortality_n_to_litr_cel_n' 
     write(10,*) col_nf%gap_mortality_n_to_litr_cel_n
     write(10,*) 'col_nf%gap_mortality_n_to_cwdn' 
     write(10,*) col_nf%gap_mortality_n_to_cwdn
     write(10,*) 'col_nf%hrv_deadstemn_to_prod100n' 
     write(10,*) col_nf%hrv_deadstemn_to_prod100n
     write(10,*) 'col_nf%hrv_deadstemn_to_prod10n' 
     write(10,*) col_nf%hrv_deadstemn_to_prod10n
     write(10,*) 'col_nf%harvest_n_to_litr_cel_n' 
     write(10,*) col_nf%harvest_n_to_litr_cel_n
     write(10,*) 'col_nf%harvest_n_to_litr_met_n' 
     write(10,*) col_nf%harvest_n_to_litr_met_n
     write(10,*) 'col_nf%harvest_n_to_litr_lig_n' 
     write(10,*) col_nf%harvest_n_to_litr_lig_n
     write(10,*) 'col_nf%harvest_n_to_cwdn' 
     write(10,*) col_nf%harvest_n_to_cwdn
     write(10,*) 'col_nf%prod10n_loss' 
     write(10,*) col_nf%prod10n_loss
     write(10,*) 'col_nf%prod100n_loss' 
     write(10,*) col_nf%prod100n_loss
     write(10,*) 'col_nf%prod1n_loss' 
     write(10,*) col_nf%prod1n_loss
     write(10,*) 'col_nf%m_n_to_litr_cel_fire' 
     write(10,*) col_nf%m_n_to_litr_cel_fire
     write(10,*) 'col_nf%fire_mortality_n_to_cwdn' 
     write(10,*) col_nf%fire_mortality_n_to_cwdn
     write(10,*) 'col_nf%m_decomp_npools_to_fire_vr' 
     write(10,*) col_nf%m_decomp_npools_to_fire_vr
     write(10,*) 'col_nf%m_n_to_litr_lig_fire' 
     write(10,*) col_nf%m_n_to_litr_lig_fire
     write(10,*) 'col_nf%m_n_to_litr_met_fire' 
     write(10,*) col_nf%m_n_to_litr_met_fire
     write(10,*) 'col_nf%decomp_npools_deposit' 
     write(10,*) col_nf%decomp_npools_deposit
     write(10,*) 'col_nf%decomp_npools_yield_vr' 
     write(10,*) col_nf%decomp_npools_yield_vr
     write(10,*) 'col_nf%decomp_npools_erode' 
     write(10,*) col_nf%decomp_npools_erode
     write(10,*) 'col_pf%pmpf_decomp_cascade' 
     write(10,*) col_pf%pmpf_decomp_cascade
     write(10,*) 'col_pf%soil_p_immob_flux_vr' 
     write(10,*) col_pf%soil_p_immob_flux_vr
     write(10,*) 'col_pf%potential_immob_p_vr' 
     write(10,*) col_pf%potential_immob_p_vr
     write(10,*) 'col_pf%net_pmin_vr' 
     write(10,*) col_pf%net_pmin_vr
     write(10,*) 'col_pf%decomp_cascade_sminp_flux_vr' 
     write(10,*) col_pf%decomp_cascade_sminp_flux_vr
     write(10,*) 'col_pf%gross_pmin_vr' 
     write(10,*) col_pf%gross_pmin_vr
     write(10,*) 'col_pf%soil_p_grossmin_flux' 
     write(10,*) col_pf%soil_p_grossmin_flux
     write(10,*) 'col_pf%decomp_cascade_ptransfer_vr' 
     write(10,*) col_pf%decomp_cascade_ptransfer_vr
     write(10,*) 'col_pf%soil_p_immob_flux' 
     write(10,*) col_pf%soil_p_immob_flux
     write(10,*) 'col_pf%sminp_to_plant_vr' 
     write(10,*) col_pf%sminp_to_plant_vr
     write(10,*) 'col_pf%sminp_to_plant' 
     write(10,*) col_pf%sminp_to_plant
     write(10,*) 'col_pf%adsorb_to_labilep_vr' 
     write(10,*) col_pf%adsorb_to_labilep_vr
     write(10,*) 'col_pf%actual_immob_p' 
     write(10,*) col_pf%actual_immob_p
     write(10,*) 'col_pf%desorb_to_solutionp_vr' 
     write(10,*) col_pf%desorb_to_solutionp_vr
     write(10,*) 'col_pf%col_plant_pdemand_vr' 
     write(10,*) col_pf%col_plant_pdemand_vr
     write(10,*) 'col_pf%actual_immob_p_vr' 
     write(10,*) col_pf%actual_immob_p_vr
     write(10,*) 'col_pf%potential_immob_p' 
     write(10,*) col_pf%potential_immob_p
     write(10,*) 'col_pf%supplement_to_sminp_vr' 
     write(10,*) col_pf%supplement_to_sminp_vr
     write(10,*) 'col_pf%gross_pmin' 
     write(10,*) col_pf%gross_pmin
     write(10,*) 'col_pf%net_pmin' 
     write(10,*) col_pf%net_pmin
     write(10,*) 'col_pf%sminp_to_plant' 
     write(10,*) col_pf%sminp_to_plant
     write(10,*) 'col_pf%plant_p_uptake_flux' 
     write(10,*) col_pf%plant_p_uptake_flux
     write(10,*) 'col_pf%sminp_to_plant_vr' 
     write(10,*) col_pf%sminp_to_plant_vr
     write(10,*) 'col_pf%hrv_cropp_to_prod1p' 
     write(10,*) col_pf%hrv_cropp_to_prod1p
     write(10,*) 'col_pf%phenology_p_to_litr_cel_p' 
     write(10,*) col_pf%phenology_p_to_litr_cel_p
     write(10,*) 'col_pf%phenology_p_to_litr_lig_p' 
     write(10,*) col_pf%phenology_p_to_litr_lig_p
     write(10,*) 'col_pf%phenology_p_to_litr_met_p' 
     write(10,*) col_pf%phenology_p_to_litr_met_p
     write(10,*) 'col_pf%decomp_ppools_sourcesink' 
     write(10,*) col_pf%decomp_ppools_sourcesink
     write(10,*) 'col_pf%decomp_ppools_transport_tendency' 
     write(10,*) col_pf%decomp_ppools_transport_tendency
     write(10,*) 'col_pf%gap_mortality_p_to_litr_met_p' 
     write(10,*) col_pf%gap_mortality_p_to_litr_met_p
     write(10,*) 'col_pf%gap_mortality_p_to_litr_cel_p' 
     write(10,*) col_pf%gap_mortality_p_to_litr_cel_p
     write(10,*) 'col_pf%gap_mortality_p_to_litr_lig_p' 
     write(10,*) col_pf%gap_mortality_p_to_litr_lig_p
     write(10,*) 'col_pf%gap_mortality_p_to_cwdp' 
     write(10,*) col_pf%gap_mortality_p_to_cwdp
     write(10,*) 'col_pf%hrv_deadstemp_to_prod10p' 
     write(10,*) col_pf%hrv_deadstemp_to_prod10p
     write(10,*) 'col_pf%harvest_p_to_cwdp' 
     write(10,*) col_pf%harvest_p_to_cwdp
     write(10,*) 'col_pf%hrv_deadstemp_to_prod100p' 
     write(10,*) col_pf%hrv_deadstemp_to_prod100p
     write(10,*) 'col_pf%harvest_p_to_litr_cel_p' 
     write(10,*) col_pf%harvest_p_to_litr_cel_p
     write(10,*) 'col_pf%harvest_p_to_litr_lig_p' 
     write(10,*) col_pf%harvest_p_to_litr_lig_p
     write(10,*) 'col_pf%harvest_p_to_litr_met_p' 
     write(10,*) col_pf%harvest_p_to_litr_met_p
     write(10,*) 'col_pf%prod10p_loss' 
     write(10,*) col_pf%prod10p_loss
     write(10,*) 'col_pf%prod100p_loss' 
     write(10,*) col_pf%prod100p_loss
     write(10,*) 'col_pf%prod1p_loss' 
     write(10,*) col_pf%prod1p_loss
     write(10,*) 'col_pf%m_p_to_litr_lig_fire' 
     write(10,*) col_pf%m_p_to_litr_lig_fire
     write(10,*) 'col_pf%fire_mortality_p_to_cwdp' 
     write(10,*) col_pf%fire_mortality_p_to_cwdp
     write(10,*) 'col_pf%m_p_to_litr_met_fire' 
     write(10,*) col_pf%m_p_to_litr_met_fire
     write(10,*) 'col_pf%m_decomp_ppools_to_fire_vr' 
     write(10,*) col_pf%m_decomp_ppools_to_fire_vr
     write(10,*) 'col_pf%m_p_to_litr_cel_fire' 
     write(10,*) col_pf%m_p_to_litr_cel_fire
     write(10,*) 'col_pf%secondp_deposit' 
     write(10,*) col_pf%secondp_deposit
     write(10,*) 'col_pf%primp_yield_vr' 
     write(10,*) col_pf%primp_yield_vr
     write(10,*) 'col_pf%occlp_erode' 
     write(10,*) col_pf%occlp_erode
     write(10,*) 'col_pf%decomp_ppools_deposit' 
     write(10,*) col_pf%decomp_ppools_deposit
     write(10,*) 'col_pf%occlp_yield_vr' 
     write(10,*) col_pf%occlp_yield_vr
     write(10,*) 'col_pf%decomp_ppools_yield_vr' 
     write(10,*) col_pf%decomp_ppools_yield_vr
     write(10,*) 'col_pf%labilep_yield_vr' 
     write(10,*) col_pf%labilep_yield_vr
     write(10,*) 'col_pf%primp_deposit' 
     write(10,*) col_pf%primp_deposit
     write(10,*) 'col_pf%labilep_erode' 
     write(10,*) col_pf%labilep_erode
     write(10,*) 'col_pf%occlp_deposit' 
     write(10,*) col_pf%occlp_deposit
     write(10,*) 'col_pf%labilep_deposit' 
     write(10,*) col_pf%labilep_deposit
     write(10,*) 'col_pf%secondp_erode' 
     write(10,*) col_pf%secondp_erode
     write(10,*) 'col_pf%secondp_yield_vr' 
     write(10,*) col_pf%secondp_yield_vr
     write(10,*) 'col_pf%decomp_ppools_erode' 
     write(10,*) col_pf%decomp_ppools_erode
     write(10,*) 'col_pf%primp_erode' 
     write(10,*) col_pf%primp_erode
     write(10,*) 'col_cf%fphr' 
     write(10,*) col_cf%fphr
     write(10,*) 'col_cf%decomp_cascade_hr_vr' 
     write(10,*) col_cf%decomp_cascade_hr_vr
     write(10,*) 'col_cf%phr_vr' 
     write(10,*) col_cf%phr_vr
     write(10,*) 'col_cf%decomp_cascade_ctransfer_vr' 
     write(10,*) col_cf%decomp_cascade_ctransfer_vr
     write(10,*) 'col_cf%fphr' 
     write(10,*) col_cf%fphr
     write(10,*) 'col_cf%phr_vr' 
     write(10,*) col_cf%phr_vr
     write(10,*) 'col_cf%hrv_cropc_to_prod1c' 
     write(10,*) col_cf%hrv_cropc_to_prod1c
     write(10,*) 'col_cf%phenology_c_to_litr_lig_c' 
     write(10,*) col_cf%phenology_c_to_litr_lig_c
     write(10,*) 'col_cf%phenology_c_to_litr_cel_c' 
     write(10,*) col_cf%phenology_c_to_litr_cel_c
     write(10,*) 'col_cf%phenology_c_to_litr_met_c' 
     write(10,*) col_cf%phenology_c_to_litr_met_c
     write(10,*) 'col_cf%phenology_c_to_litr_lig_c' 
     write(10,*) col_cf%phenology_c_to_litr_lig_c
     write(10,*) 'col_cf%phenology_c_to_litr_cel_c' 
     write(10,*) col_cf%phenology_c_to_litr_cel_c
     write(10,*) 'col_cf%phenology_c_to_litr_met_c' 
     write(10,*) col_cf%phenology_c_to_litr_met_c
     write(10,*) 'col_cf%decomp_cpools_sourcesink' 
     write(10,*) col_cf%decomp_cpools_sourcesink
     write(10,*) 'col_cf%decomp_cpools_transport_tendency' 
     write(10,*) col_cf%decomp_cpools_transport_tendency
     write(10,*) 'col_cf%gap_mortality_c_to_litr_lig_c' 
     write(10,*) col_cf%gap_mortality_c_to_litr_lig_c
     write(10,*) 'col_cf%gap_mortality_c_to_litr_cel_c' 
     write(10,*) col_cf%gap_mortality_c_to_litr_cel_c
     write(10,*) 'col_cf%gap_mortality_c_to_cwdc' 
     write(10,*) col_cf%gap_mortality_c_to_cwdc
     write(10,*) 'col_cf%gap_mortality_c_to_litr_met_c' 
     write(10,*) col_cf%gap_mortality_c_to_litr_met_c
     write(10,*) 'col_cf%gap_mortality_c_to_litr_lig_c' 
     write(10,*) col_cf%gap_mortality_c_to_litr_lig_c
     write(10,*) 'col_cf%gap_mortality_c_to_litr_cel_c' 
     write(10,*) col_cf%gap_mortality_c_to_litr_cel_c
     write(10,*) 'col_cf%gap_mortality_c_to_cwdc' 
     write(10,*) col_cf%gap_mortality_c_to_cwdc
     write(10,*) 'col_cf%gap_mortality_c_to_litr_met_c' 
     write(10,*) col_cf%gap_mortality_c_to_litr_met_c
     write(10,*) 'col_cf%harvest_c_to_litr_lig_c' 
     write(10,*) col_cf%harvest_c_to_litr_lig_c
     write(10,*) 'col_cf%hrv_deadstemc_to_prod10c' 
     write(10,*) col_cf%hrv_deadstemc_to_prod10c
     write(10,*) 'col_cf%hrv_deadstemc_to_prod100c' 
     write(10,*) col_cf%hrv_deadstemc_to_prod100c
     write(10,*) 'col_cf%harvest_c_to_litr_cel_c' 
     write(10,*) col_cf%harvest_c_to_litr_cel_c
     write(10,*) 'col_cf%harvest_c_to_cwdc' 
     write(10,*) col_cf%harvest_c_to_cwdc
     write(10,*) 'col_cf%harvest_c_to_litr_met_c' 
     write(10,*) col_cf%harvest_c_to_litr_met_c
     write(10,*) 'col_cf%harvest_c_to_litr_lig_c' 
     write(10,*) col_cf%harvest_c_to_litr_lig_c
     write(10,*) 'col_cf%hrv_deadstemc_to_prod10c' 
     write(10,*) col_cf%hrv_deadstemc_to_prod10c
     write(10,*) 'col_cf%hrv_deadstemc_to_prod100c' 
     write(10,*) col_cf%hrv_deadstemc_to_prod100c
     write(10,*) 'col_cf%harvest_c_to_litr_cel_c' 
     write(10,*) col_cf%harvest_c_to_litr_cel_c
     write(10,*) 'col_cf%harvest_c_to_cwdc' 
     write(10,*) col_cf%harvest_c_to_cwdc
     write(10,*) 'col_cf%harvest_c_to_litr_met_c' 
     write(10,*) col_cf%harvest_c_to_litr_met_c
     write(10,*) 'col_cf%prod10c_loss' 
     write(10,*) col_cf%prod10c_loss
     write(10,*) 'col_cf%prod100c_loss' 
     write(10,*) col_cf%prod100c_loss
     write(10,*) 'col_cf%prod1c_loss' 
     write(10,*) col_cf%prod1c_loss
     write(10,*) 'col_cf%fire_mortality_c_to_cwdc' 
     write(10,*) col_cf%fire_mortality_c_to_cwdc
     write(10,*) 'col_cf%m_c_to_litr_met_fire' 
     write(10,*) col_cf%m_c_to_litr_met_fire
     write(10,*) 'col_cf%m_c_to_litr_cel_fire' 
     write(10,*) col_cf%m_c_to_litr_cel_fire
     write(10,*) 'col_cf%m_c_to_litr_lig_fire' 
     write(10,*) col_cf%m_c_to_litr_lig_fire
     write(10,*) 'col_cf%m_decomp_cpools_to_fire_vr' 
     write(10,*) col_cf%m_decomp_cpools_to_fire_vr
     write(10,*) 'col_cf%somc_fire' 
     write(10,*) col_cf%somc_fire
     write(10,*) 'col_cf%decomp_cpools_yield_vr' 
     write(10,*) col_cf%decomp_cpools_yield_vr
     write(10,*) 'col_cf%decomp_cpools_deposit' 
     write(10,*) col_cf%decomp_cpools_deposit
     write(10,*) 'col_cf%decomp_cpools_erode' 
     write(10,*) col_cf%decomp_cpools_erode
     write(10,*) 'cnstate_vars%fpg_p_vr_col' 
     write(10,*) cnstate_vars%fpg_p_vr_col
     write(10,*) 'cnstate_vars%cn_scalar' 
     write(10,*) cnstate_vars%cn_scalar
     write(10,*) 'cnstate_vars%fpi_p_vr_col' 
     write(10,*) cnstate_vars%fpi_p_vr_col
     write(10,*) 'cnstate_vars%fpg_p_col' 
     write(10,*) cnstate_vars%fpg_p_col
     write(10,*) 'cnstate_vars%fpi_vr_col' 
     write(10,*) cnstate_vars%fpi_vr_col
     write(10,*) 'cnstate_vars%fpg_nh4_vr_col' 
     write(10,*) cnstate_vars%fpg_nh4_vr_col
     write(10,*) 'cnstate_vars%fpg_col' 
     write(10,*) cnstate_vars%fpg_col
     write(10,*) 'cnstate_vars%cp_scalar' 
     write(10,*) cnstate_vars%cp_scalar
     write(10,*) 'cnstate_vars%fpi_p_col' 
     write(10,*) cnstate_vars%fpi_p_col
     write(10,*) 'cnstate_vars%fpi_col' 
     write(10,*) cnstate_vars%fpi_col
     write(10,*) 'cnstate_vars%fpg_no3_vr_col' 
     write(10,*) cnstate_vars%fpg_no3_vr_col
     write(10,*) 'cnstate_vars%fpg_vr_col' 
     write(10,*) cnstate_vars%fpg_vr_col
     write(10,*) 'cnstate_vars%fpi_vr_col' 
     write(10,*) cnstate_vars%fpi_vr_col
     write(10,*) 'cnstate_vars%fpg_col' 
     write(10,*) cnstate_vars%fpg_col
     write(10,*) 'cnstate_vars%fpi_col' 
     write(10,*) cnstate_vars%fpi_col
     write(10,*) 'cnstate_vars%downreg_patch' 
     write(10,*) cnstate_vars%downreg_patch
     write(10,*) 'cnstate_vars%n_allometry_patch' 
     write(10,*) cnstate_vars%n_allometry_patch
     write(10,*) 'cnstate_vars%p_allometry_patch' 
     write(10,*) cnstate_vars%p_allometry_patch
     write(10,*) 'cnstate_vars%c_allometry_patch' 
     write(10,*) cnstate_vars%c_allometry_patch
     write(10,*) 'cnstate_vars%tempavg_t2m_patch' 
     write(10,*) cnstate_vars%tempavg_t2m_patch
     write(10,*) 'cnstate_vars%bglfr_froot_patch' 
     write(10,*) cnstate_vars%bglfr_froot_patch
     write(10,*) 'cnstate_vars%lgsf_patch' 
     write(10,*) cnstate_vars%lgsf_patch
     write(10,*) 'cnstate_vars%bglfr_leaf_patch' 
     write(10,*) cnstate_vars%bglfr_leaf_patch
     write(10,*) 'cnstate_vars%bgtr_patch' 
     write(10,*) cnstate_vars%bgtr_patch
     write(10,*) 'cnstate_vars%lgsf_patch' 
     write(10,*) cnstate_vars%lgsf_patch
     write(10,*) 'cnstate_vars%onset_flag_patch' 
     write(10,*) cnstate_vars%onset_flag_patch
     write(10,*) 'cnstate_vars%days_active_patch' 
     write(10,*) cnstate_vars%days_active_patch
     write(10,*) 'cnstate_vars%bgtr_patch' 
     write(10,*) cnstate_vars%bgtr_patch
     write(10,*) 'cnstate_vars%onset_counter_patch' 
     write(10,*) cnstate_vars%onset_counter_patch
     write(10,*) 'cnstate_vars%onset_gdd_patch' 
     write(10,*) cnstate_vars%onset_gdd_patch
     write(10,*) 'cnstate_vars%offset_counter_patch' 
     write(10,*) cnstate_vars%offset_counter_patch
     write(10,*) 'cnstate_vars%onset_gddflag_patch' 
     write(10,*) cnstate_vars%onset_gddflag_patch
     write(10,*) 'cnstate_vars%bglfr_froot_patch' 
     write(10,*) cnstate_vars%bglfr_froot_patch
     write(10,*) 'cnstate_vars%dormant_flag_patch' 
     write(10,*) cnstate_vars%dormant_flag_patch
     write(10,*) 'cnstate_vars%bglfr_leaf_patch' 
     write(10,*) cnstate_vars%bglfr_leaf_patch
     write(10,*) 'cnstate_vars%offset_flag_patch' 
     write(10,*) cnstate_vars%offset_flag_patch
     write(10,*) 'cnstate_vars%offset_fdd_patch' 
     write(10,*) cnstate_vars%offset_fdd_patch
     write(10,*) 'cnstate_vars%lgsf_patch' 
     write(10,*) cnstate_vars%lgsf_patch
     write(10,*) 'cnstate_vars%onset_flag_patch' 
     write(10,*) cnstate_vars%onset_flag_patch
     write(10,*) 'cnstate_vars%days_active_patch' 
     write(10,*) cnstate_vars%days_active_patch
     write(10,*) 'cnstate_vars%bgtr_patch' 
     write(10,*) cnstate_vars%bgtr_patch
     write(10,*) 'cnstate_vars%onset_counter_patch' 
     write(10,*) cnstate_vars%onset_counter_patch
     write(10,*) 'cnstate_vars%onset_swi_patch' 
     write(10,*) cnstate_vars%onset_swi_patch
     write(10,*) 'cnstate_vars%onset_gdd_patch' 
     write(10,*) cnstate_vars%onset_gdd_patch
     write(10,*) 'cnstate_vars%offset_counter_patch' 
     write(10,*) cnstate_vars%offset_counter_patch
     write(10,*) 'cnstate_vars%onset_gddflag_patch' 
     write(10,*) cnstate_vars%onset_gddflag_patch
     write(10,*) 'cnstate_vars%offset_swi_patch' 
     write(10,*) cnstate_vars%offset_swi_patch
     write(10,*) 'cnstate_vars%bglfr_froot_patch' 
     write(10,*) cnstate_vars%bglfr_froot_patch
     write(10,*) 'cnstate_vars%dormant_flag_patch' 
     write(10,*) cnstate_vars%dormant_flag_patch
     write(10,*) 'cnstate_vars%bglfr_leaf_patch' 
     write(10,*) cnstate_vars%bglfr_leaf_patch
     write(10,*) 'cnstate_vars%onset_fdd_patch' 
     write(10,*) cnstate_vars%onset_fdd_patch
     write(10,*) 'cnstate_vars%offset_flag_patch' 
     write(10,*) cnstate_vars%offset_flag_patch
     write(10,*) 'cnstate_vars%lgsf_patch' 
     write(10,*) cnstate_vars%lgsf_patch
     write(10,*) 'cnstate_vars%huigrain_patch' 
     write(10,*) cnstate_vars%huigrain_patch
     write(10,*) 'cnstate_vars%bgtr_patch' 
     write(10,*) cnstate_vars%bgtr_patch
     write(10,*) 'cnstate_vars%huileaf_patch' 
     write(10,*) cnstate_vars%huileaf_patch
     write(10,*) 'cnstate_vars%onset_flag_patch' 
     write(10,*) cnstate_vars%onset_flag_patch
     write(10,*) 'cnstate_vars%onset_counter_patch' 
     write(10,*) cnstate_vars%onset_counter_patch
     write(10,*) 'cnstate_vars%idop_patch' 
     write(10,*) cnstate_vars%idop_patch
     write(10,*) 'cnstate_vars%bglfr_froot_patch' 
     write(10,*) cnstate_vars%bglfr_froot_patch
     write(10,*) 'cnstate_vars%gddmaturity_patch' 
     write(10,*) cnstate_vars%gddmaturity_patch
     write(10,*) 'cnstate_vars%offset_counter_patch' 
     write(10,*) cnstate_vars%offset_counter_patch
     write(10,*) 'cnstate_vars%hdidx_patch' 
     write(10,*) cnstate_vars%hdidx_patch
     write(10,*) 'cnstate_vars%bglfr_leaf_patch' 
     write(10,*) cnstate_vars%bglfr_leaf_patch
     write(10,*) 'cnstate_vars%offset_flag_patch' 
     write(10,*) cnstate_vars%offset_flag_patch
     write(10,*) 'cnstate_vars%cumvd_patch' 
     write(10,*) cnstate_vars%cumvd_patch
     write(10,*) 'cnstate_vars%hdidx_patch' 
     write(10,*) cnstate_vars%hdidx_patch
     write(10,*) 'cnstate_vars%huigrain_patch' 
     write(10,*) cnstate_vars%huigrain_patch
     write(10,*) 'cnstate_vars%cumvd_patch' 
     write(10,*) cnstate_vars%cumvd_patch
     write(10,*) 'cnstate_vars%gddmaturity_patch' 
     write(10,*) cnstate_vars%gddmaturity_patch
     write(10,*) 'cnstate_vars%r_mort_cal_patch' 
     write(10,*) cnstate_vars%r_mort_cal_patch
     write(10,*) 'cnstate_vars%nfire_col' 
     write(10,*) cnstate_vars%nfire_col
     write(10,*) 'cnstate_vars%lgdp1_col' 
     write(10,*) cnstate_vars%lgdp1_col
     write(10,*) 'cnstate_vars%lfc_col' 
     write(10,*) cnstate_vars%lfc_col
     write(10,*) 'cnstate_vars%dtrotr_col' 
     write(10,*) cnstate_vars%dtrotr_col
     write(10,*) 'cnstate_vars%lfwt_col' 
     write(10,*) cnstate_vars%lfwt_col
     write(10,*) 'cnstate_vars%burndate_patch' 
     write(10,*) cnstate_vars%burndate_patch
     write(10,*) 'cnstate_vars%trotr1_col' 
     write(10,*) cnstate_vars%trotr1_col
     write(10,*) 'cnstate_vars%fsr_col' 
     write(10,*) cnstate_vars%fsr_col
     write(10,*) 'cnstate_vars%baf_crop_col' 
     write(10,*) cnstate_vars%baf_crop_col
     write(10,*) 'cnstate_vars%baf_peatf_col' 
     write(10,*) cnstate_vars%baf_peatf_col
     write(10,*) 'cnstate_vars%fbac1_col' 
     write(10,*) cnstate_vars%fbac1_col
     write(10,*) 'cnstate_vars%cropf_col' 
     write(10,*) cnstate_vars%cropf_col
     write(10,*) 'cnstate_vars%fd_col' 
     write(10,*) cnstate_vars%fd_col
     write(10,*) 'cnstate_vars%lpop_col' 
     write(10,*) cnstate_vars%lpop_col
     write(10,*) 'cnstate_vars%farea_burned_col' 
     write(10,*) cnstate_vars%farea_burned_col
     write(10,*) 'cnstate_vars%trotr2_col' 
     write(10,*) cnstate_vars%trotr2_col
     write(10,*) 'cnstate_vars%fbac_col' 
     write(10,*) cnstate_vars%fbac_col
     write(10,*) 'cnstate_vars%wtlf_col' 
     write(10,*) cnstate_vars%wtlf_col
     write(10,*) 'cnstate_vars%lgdp_col' 
     write(10,*) cnstate_vars%lgdp_col
     write(10,*) 'cnstate_vars%lfc_col' 
     write(10,*) cnstate_vars%lfc_col
     write(10,*) 'cnstate_vars%lfc2_col' 
     write(10,*) cnstate_vars%lfc2_col
     write(10,*) 'cnstate_vars%rc14_atm_patch' 
     write(10,*) cnstate_vars%rc14_atm_patch
     write(10,*) 'veg_nf%smin_no3_to_plant' 
     write(10,*) veg_nf%smin_no3_to_plant
     write(10,*) 'veg_nf%plant_ndemand_vr' 
     write(10,*) veg_nf%plant_ndemand_vr
     write(10,*) 'veg_nf%sminn_to_plant' 
     write(10,*) veg_nf%sminn_to_plant
     write(10,*) 'veg_nf%plant_no3demand_vr' 
     write(10,*) veg_nf%plant_no3demand_vr
     write(10,*) 'veg_nf%smin_nh4_to_plant' 
     write(10,*) veg_nf%smin_nh4_to_plant
     write(10,*) 'veg_nf%plant_nh4demand_vr' 
     write(10,*) veg_nf%plant_nh4demand_vr
     write(10,*) 'veg_nf%npool_to_frootn' 
     write(10,*) veg_nf%npool_to_frootn
     write(10,*) 'veg_nf%retransn_to_npool' 
     write(10,*) veg_nf%retransn_to_npool
     write(10,*) 'veg_nf%npool_to_livecrootn' 
     write(10,*) veg_nf%npool_to_livecrootn
     write(10,*) 'veg_nf%npool_to_livecrootn_storage' 
     write(10,*) veg_nf%npool_to_livecrootn_storage
     write(10,*) 'veg_nf%npool_to_deadcrootn_storage' 
     write(10,*) veg_nf%npool_to_deadcrootn_storage
     write(10,*) 'veg_nf%npool_to_livestemn_storage' 
     write(10,*) veg_nf%npool_to_livestemn_storage
     write(10,*) 'veg_nf%npool_to_leafn_storage' 
     write(10,*) veg_nf%npool_to_leafn_storage
     write(10,*) 'veg_nf%npool_to_frootn_storage' 
     write(10,*) veg_nf%npool_to_frootn_storage
     write(10,*) 'veg_nf%plant_nalloc' 
     write(10,*) veg_nf%plant_nalloc
     write(10,*) 'veg_nf%npool_to_grainn' 
     write(10,*) veg_nf%npool_to_grainn
     write(10,*) 'veg_nf%npool_to_leafn' 
     write(10,*) veg_nf%npool_to_leafn
     write(10,*) 'veg_nf%npool_to_deadstemn' 
     write(10,*) veg_nf%npool_to_deadstemn
     write(10,*) 'veg_nf%npool_to_deadstemn_storage' 
     write(10,*) veg_nf%npool_to_deadstemn_storage
     write(10,*) 'veg_nf%npool_to_deadcrootn' 
     write(10,*) veg_nf%npool_to_deadcrootn
     write(10,*) 'veg_nf%supplement_to_plantn' 
     write(10,*) veg_nf%supplement_to_plantn
     write(10,*) 'veg_nf%sminn_to_npool' 
     write(10,*) veg_nf%sminn_to_npool
     write(10,*) 'veg_nf%npool_to_livestemn' 
     write(10,*) veg_nf%npool_to_livestemn
     write(10,*) 'veg_nf%avail_retransn' 
     write(10,*) veg_nf%avail_retransn
     write(10,*) 'veg_nf%npool_to_grainn_storage' 
     write(10,*) veg_nf%npool_to_grainn_storage
     write(10,*) 'veg_nf%leafn_storage_to_xfer' 
     write(10,*) veg_nf%leafn_storage_to_xfer
     write(10,*) 'veg_nf%livestemn_storage_to_xfer' 
     write(10,*) veg_nf%livestemn_storage_to_xfer
     write(10,*) 'veg_nf%deadcrootn_xfer_to_deadcrootn' 
     write(10,*) veg_nf%deadcrootn_xfer_to_deadcrootn
     write(10,*) 'veg_nf%frootn_xfer_to_frootn' 
     write(10,*) veg_nf%frootn_xfer_to_frootn
     write(10,*) 'veg_nf%frootn_storage_to_xfer' 
     write(10,*) veg_nf%frootn_storage_to_xfer
     write(10,*) 'veg_nf%deadstemn_xfer_to_deadstemn' 
     write(10,*) veg_nf%deadstemn_xfer_to_deadstemn
     write(10,*) 'veg_nf%deadcrootn_storage_to_xfer' 
     write(10,*) veg_nf%deadcrootn_storage_to_xfer
     write(10,*) 'veg_nf%leafn_xfer_to_leafn' 
     write(10,*) veg_nf%leafn_xfer_to_leafn
     write(10,*) 'veg_nf%deadstemn_storage_to_xfer' 
     write(10,*) veg_nf%deadstemn_storage_to_xfer
     write(10,*) 'veg_nf%livecrootn_xfer_to_livecrootn' 
     write(10,*) veg_nf%livecrootn_xfer_to_livecrootn
     write(10,*) 'veg_nf%livestemn_xfer_to_livestemn' 
     write(10,*) veg_nf%livestemn_xfer_to_livestemn
     write(10,*) 'veg_nf%livecrootn_storage_to_xfer' 
     write(10,*) veg_nf%livecrootn_storage_to_xfer
     write(10,*) 'veg_nf%leafn_storage_to_xfer' 
     write(10,*) veg_nf%leafn_storage_to_xfer
     write(10,*) 'veg_nf%livestemn_storage_to_xfer' 
     write(10,*) veg_nf%livestemn_storage_to_xfer
     write(10,*) 'veg_nf%deadcrootn_xfer_to_deadcrootn' 
     write(10,*) veg_nf%deadcrootn_xfer_to_deadcrootn
     write(10,*) 'veg_nf%frootn_xfer_to_frootn' 
     write(10,*) veg_nf%frootn_xfer_to_frootn
     write(10,*) 'veg_nf%frootn_storage_to_xfer' 
     write(10,*) veg_nf%frootn_storage_to_xfer
     write(10,*) 'veg_nf%deadstemn_xfer_to_deadstemn' 
     write(10,*) veg_nf%deadstemn_xfer_to_deadstemn
     write(10,*) 'veg_nf%deadcrootn_storage_to_xfer' 
     write(10,*) veg_nf%deadcrootn_storage_to_xfer
     write(10,*) 'veg_nf%leafn_xfer_to_leafn' 
     write(10,*) veg_nf%leafn_xfer_to_leafn
     write(10,*) 'veg_nf%deadstemn_storage_to_xfer' 
     write(10,*) veg_nf%deadstemn_storage_to_xfer
     write(10,*) 'veg_nf%livecrootn_xfer_to_livecrootn' 
     write(10,*) veg_nf%livecrootn_xfer_to_livecrootn
     write(10,*) 'veg_nf%livestemn_xfer_to_livestemn' 
     write(10,*) veg_nf%livestemn_xfer_to_livestemn
     write(10,*) 'veg_nf%livecrootn_storage_to_xfer' 
     write(10,*) veg_nf%livecrootn_storage_to_xfer
     write(10,*) 'veg_nf%fert_counter' 
     write(10,*) veg_nf%fert_counter
     write(10,*) 'veg_nf%crop_seedn_to_leaf' 
     write(10,*) veg_nf%crop_seedn_to_leaf
     write(10,*) 'veg_nf%fert' 
     write(10,*) veg_nf%fert
     write(10,*) 'veg_nf%deadcrootn_xfer_to_deadcrootn' 
     write(10,*) veg_nf%deadcrootn_xfer_to_deadcrootn
     write(10,*) 'veg_nf%frootn_xfer_to_frootn' 
     write(10,*) veg_nf%frootn_xfer_to_frootn
     write(10,*) 'veg_nf%deadstemn_xfer_to_deadstemn' 
     write(10,*) veg_nf%deadstemn_xfer_to_deadstemn
     write(10,*) 'veg_nf%leafn_xfer_to_leafn' 
     write(10,*) veg_nf%leafn_xfer_to_leafn
     write(10,*) 'veg_nf%livecrootn_xfer_to_livecrootn' 
     write(10,*) veg_nf%livecrootn_xfer_to_livecrootn
     write(10,*) 'veg_nf%livestemn_xfer_to_livestemn' 
     write(10,*) veg_nf%livestemn_xfer_to_livestemn
     write(10,*) 'veg_nf%hrv_livestemn_to_prod1n' 
     write(10,*) veg_nf%hrv_livestemn_to_prod1n
     write(10,*) 'veg_nf%hrv_grainn_to_prod1n' 
     write(10,*) veg_nf%hrv_grainn_to_prod1n
     write(10,*) 'veg_nf%hrv_leafn_to_prod1n' 
     write(10,*) veg_nf%hrv_leafn_to_prod1n
     write(10,*) 'veg_nf%hrv_cropn_to_prod1n' 
     write(10,*) veg_nf%hrv_cropn_to_prod1n
     write(10,*) 'veg_nf%leafn_to_retransn' 
     write(10,*) veg_nf%leafn_to_retransn
     write(10,*) 'veg_nf%frootn_to_litter' 
     write(10,*) veg_nf%frootn_to_litter
     write(10,*) 'veg_nf%leafn_to_litter' 
     write(10,*) veg_nf%leafn_to_litter
     write(10,*) 'veg_nf%livestemn_to_litter' 
     write(10,*) veg_nf%livestemn_to_litter
     write(10,*) 'veg_nf%frootn_to_litter' 
     write(10,*) veg_nf%frootn_to_litter
     write(10,*) 'veg_nf%leafn_to_litter' 
     write(10,*) veg_nf%leafn_to_litter
     write(10,*) 'veg_nf%leafn_to_retransn' 
     write(10,*) veg_nf%leafn_to_retransn
     write(10,*) 'veg_nf%livecrootn_to_retransn' 
     write(10,*) veg_nf%livecrootn_to_retransn
     write(10,*) 'veg_nf%livestemn_to_deadstemn' 
     write(10,*) veg_nf%livestemn_to_deadstemn
     write(10,*) 'veg_nf%livecrootn_to_deadcrootn' 
     write(10,*) veg_nf%livecrootn_to_deadcrootn
     write(10,*) 'veg_nf%livestemn_to_retransn' 
     write(10,*) veg_nf%livestemn_to_retransn
     write(10,*) 'veg_nf%npool_to_livecrootn' 
     write(10,*) veg_nf%npool_to_livecrootn
     write(10,*) 'veg_nf%deadcrootn_xfer_to_deadcrootn' 
     write(10,*) veg_nf%deadcrootn_xfer_to_deadcrootn
     write(10,*) 'veg_nf%livestemn_to_retransn' 
     write(10,*) veg_nf%livestemn_to_retransn
     write(10,*) 'veg_nf%npool_to_frootn_storage' 
     write(10,*) veg_nf%npool_to_frootn_storage
     write(10,*) 'veg_nf%livecrootn_xfer_to_livecrootn' 
     write(10,*) veg_nf%livecrootn_xfer_to_livecrootn
     write(10,*) 'veg_nf%grainn_to_food' 
     write(10,*) veg_nf%grainn_to_food
     write(10,*) 'veg_nf%npool_to_grainn_storage' 
     write(10,*) veg_nf%npool_to_grainn_storage
     write(10,*) 'veg_nf%npool_to_deadcrootn_storage' 
     write(10,*) veg_nf%npool_to_deadcrootn_storage
     write(10,*) 'veg_nf%npool_to_deadstemn_storage' 
     write(10,*) veg_nf%npool_to_deadstemn_storage
     write(10,*) 'veg_nf%npool_to_deadcrootn' 
     write(10,*) veg_nf%npool_to_deadcrootn
     write(10,*) 'veg_nf%livestemn_to_litter' 
     write(10,*) veg_nf%livestemn_to_litter
     write(10,*) 'veg_nf%deadstemn_storage_to_xfer' 
     write(10,*) veg_nf%deadstemn_storage_to_xfer
     write(10,*) 'veg_nf%livecrootn_to_retransn' 
     write(10,*) veg_nf%livecrootn_to_retransn
     write(10,*) 'veg_nf%livestemn_storage_to_xfer' 
     write(10,*) veg_nf%livestemn_storage_to_xfer
     write(10,*) 'veg_nf%npool_to_livestemn_storage' 
     write(10,*) veg_nf%npool_to_livestemn_storage
     write(10,*) 'veg_nf%npool_to_leafn_storage' 
     write(10,*) veg_nf%npool_to_leafn_storage
     write(10,*) 'veg_nf%npool_to_livestemn' 
     write(10,*) veg_nf%npool_to_livestemn
     write(10,*) 'veg_nf%grainn_xfer_to_grainn' 
     write(10,*) veg_nf%grainn_xfer_to_grainn
     write(10,*) 'veg_nf%livecrootn_to_deadcrootn' 
     write(10,*) veg_nf%livecrootn_to_deadcrootn
     write(10,*) 'veg_nf%deadstemn_xfer_to_deadstemn' 
     write(10,*) veg_nf%deadstemn_xfer_to_deadstemn
     write(10,*) 'veg_nf%deadcrootn_storage_to_xfer' 
     write(10,*) veg_nf%deadcrootn_storage_to_xfer
     write(10,*) 'veg_nf%livestemn_xfer_to_livestemn' 
     write(10,*) veg_nf%livestemn_xfer_to_livestemn
     write(10,*) 'veg_nf%livestemn_to_deadstemn' 
     write(10,*) veg_nf%livestemn_to_deadstemn
     write(10,*) 'veg_nf%npool_to_frootn' 
     write(10,*) veg_nf%npool_to_frootn
     write(10,*) 'veg_nf%npool_to_livecrootn_storage' 
     write(10,*) veg_nf%npool_to_livecrootn_storage
     write(10,*) 'veg_nf%npool_to_grainn' 
     write(10,*) veg_nf%npool_to_grainn
     write(10,*) 'veg_nf%frootn_to_retransn' 
     write(10,*) veg_nf%frootn_to_retransn
     write(10,*) 'veg_nf%npool_to_deadstemn' 
     write(10,*) veg_nf%npool_to_deadstemn
     write(10,*) 'veg_nf%livecrootn_storage_to_xfer' 
     write(10,*) veg_nf%livecrootn_storage_to_xfer
     write(10,*) 'veg_nf%npool_to_leafn' 
     write(10,*) veg_nf%npool_to_leafn
     write(10,*) 'veg_nf%m_deadstemn_storage_to_litter' 
     write(10,*) veg_nf%m_deadstemn_storage_to_litter
     write(10,*) 'veg_nf%m_deadstemn_xfer_to_litter' 
     write(10,*) veg_nf%m_deadstemn_xfer_to_litter
     write(10,*) 'veg_nf%m_livecrootn_to_litter' 
     write(10,*) veg_nf%m_livecrootn_to_litter
     write(10,*) 'veg_nf%m_leafn_to_litter' 
     write(10,*) veg_nf%m_leafn_to_litter
     write(10,*) 'veg_nf%m_deadcrootn_storage_to_litter' 
     write(10,*) veg_nf%m_deadcrootn_storage_to_litter
     write(10,*) 'veg_nf%m_deadstemn_to_litter' 
     write(10,*) veg_nf%m_deadstemn_to_litter
     write(10,*) 'veg_nf%m_npool_to_litter' 
     write(10,*) veg_nf%m_npool_to_litter
     write(10,*) 'veg_nf%m_frootn_storage_to_litter' 
     write(10,*) veg_nf%m_frootn_storage_to_litter
     write(10,*) 'veg_nf%m_deadcrootn_xfer_to_litter' 
     write(10,*) veg_nf%m_deadcrootn_xfer_to_litter
     write(10,*) 'veg_nf%m_livestemn_xfer_to_litter' 
     write(10,*) veg_nf%m_livestemn_xfer_to_litter
     write(10,*) 'veg_nf%m_livecrootn_storage_to_litter' 
     write(10,*) veg_nf%m_livecrootn_storage_to_litter
     write(10,*) 'veg_nf%m_livestemn_to_litter' 
     write(10,*) veg_nf%m_livestemn_to_litter
     write(10,*) 'veg_nf%m_frootn_to_litter' 
     write(10,*) veg_nf%m_frootn_to_litter
     write(10,*) 'veg_nf%m_frootn_xfer_to_litter' 
     write(10,*) veg_nf%m_frootn_xfer_to_litter
     write(10,*) 'veg_nf%m_retransn_to_litter' 
     write(10,*) veg_nf%m_retransn_to_litter
     write(10,*) 'veg_nf%m_livecrootn_xfer_to_litter' 
     write(10,*) veg_nf%m_livecrootn_xfer_to_litter
     write(10,*) 'veg_nf%m_deadcrootn_to_litter' 
     write(10,*) veg_nf%m_deadcrootn_to_litter
     write(10,*) 'veg_nf%m_leafn_xfer_to_litter' 
     write(10,*) veg_nf%m_leafn_xfer_to_litter
     write(10,*) 'veg_nf%m_livestemn_storage_to_litter' 
     write(10,*) veg_nf%m_livestemn_storage_to_litter
     write(10,*) 'veg_nf%m_leafn_storage_to_litter' 
     write(10,*) veg_nf%m_leafn_storage_to_litter
     write(10,*) 'veg_nf%m_deadstemn_storage_to_litter' 
     write(10,*) veg_nf%m_deadstemn_storage_to_litter
     write(10,*) 'veg_nf%m_leafn_to_litter' 
     write(10,*) veg_nf%m_leafn_to_litter
     write(10,*) 'veg_nf%m_deadcrootn_storage_to_litter' 
     write(10,*) veg_nf%m_deadcrootn_storage_to_litter
     write(10,*) 'veg_nf%m_frootn_storage_to_litter' 
     write(10,*) veg_nf%m_frootn_storage_to_litter
     write(10,*) 'veg_nf%m_frootn_to_litter' 
     write(10,*) veg_nf%m_frootn_to_litter
     write(10,*) 'veg_nf%m_livecrootn_storage_to_litter' 
     write(10,*) veg_nf%m_livecrootn_storage_to_litter
     write(10,*) 'veg_nf%m_leafn_xfer_to_litter' 
     write(10,*) veg_nf%m_leafn_xfer_to_litter
     write(10,*) 'veg_nf%m_retransn_to_litter' 
     write(10,*) veg_nf%m_retransn_to_litter
     write(10,*) 'veg_nf%m_livestemn_storage_to_litter' 
     write(10,*) veg_nf%m_livestemn_storage_to_litter
     write(10,*) 'veg_nf%m_frootn_xfer_to_litter' 
     write(10,*) veg_nf%m_frootn_xfer_to_litter
     write(10,*) 'veg_nf%m_leafn_storage_to_litter' 
     write(10,*) veg_nf%m_leafn_storage_to_litter
     write(10,*) 'veg_nf%hrv_leafn_to_litter' 
     write(10,*) veg_nf%hrv_leafn_to_litter
     write(10,*) 'veg_nf%hrv_frootn_xfer_to_litter' 
     write(10,*) veg_nf%hrv_frootn_xfer_to_litter
     write(10,*) 'veg_nf%hrv_livestemn_xfer_to_litter' 
     write(10,*) veg_nf%hrv_livestemn_xfer_to_litter
     write(10,*) 'veg_nf%hrv_deadstemn_storage_to_litter' 
     write(10,*) veg_nf%hrv_deadstemn_storage_to_litter
     write(10,*) 'veg_nf%hrv_deadcrootn_xfer_to_litter' 
     write(10,*) veg_nf%hrv_deadcrootn_xfer_to_litter
     write(10,*) 'veg_nf%hrv_deadstemn_to_prod100n' 
     write(10,*) veg_nf%hrv_deadstemn_to_prod100n
     write(10,*) 'veg_nf%hrv_livecrootn_storage_to_litter' 
     write(10,*) veg_nf%hrv_livecrootn_storage_to_litter
     write(10,*) 'veg_nf%hrv_deadcrootn_to_litter' 
     write(10,*) veg_nf%hrv_deadcrootn_to_litter
     write(10,*) 'veg_nf%hrv_livestemn_to_litter' 
     write(10,*) veg_nf%hrv_livestemn_to_litter
     write(10,*) 'veg_nf%hrv_retransn_to_litter' 
     write(10,*) veg_nf%hrv_retransn_to_litter
     write(10,*) 'veg_nf%hrv_livecrootn_to_litter' 
     write(10,*) veg_nf%hrv_livecrootn_to_litter
     write(10,*) 'veg_nf%hrv_npool_to_litter' 
     write(10,*) veg_nf%hrv_npool_to_litter
     write(10,*) 'veg_nf%hrv_deadcrootn_storage_to_litter' 
     write(10,*) veg_nf%hrv_deadcrootn_storage_to_litter
     write(10,*) 'veg_nf%hrv_deadstemn_xfer_to_litter' 
     write(10,*) veg_nf%hrv_deadstemn_xfer_to_litter
     write(10,*) 'veg_nf%hrv_livestemn_storage_to_litter' 
     write(10,*) veg_nf%hrv_livestemn_storage_to_litter
     write(10,*) 'veg_nf%hrv_deadstemn_to_prod10n' 
     write(10,*) veg_nf%hrv_deadstemn_to_prod10n
     write(10,*) 'veg_nf%hrv_leafn_xfer_to_litter' 
     write(10,*) veg_nf%hrv_leafn_xfer_to_litter
     write(10,*) 'veg_nf%hrv_livecrootn_xfer_to_litter' 
     write(10,*) veg_nf%hrv_livecrootn_xfer_to_litter
     write(10,*) 'veg_nf%hrv_frootn_to_litter' 
     write(10,*) veg_nf%hrv_frootn_to_litter
     write(10,*) 'veg_nf%hrv_leafn_storage_to_litter' 
     write(10,*) veg_nf%hrv_leafn_storage_to_litter
     write(10,*) 'veg_nf%hrv_frootn_storage_to_litter' 
     write(10,*) veg_nf%hrv_frootn_storage_to_litter
     write(10,*) 'veg_nf%m_livestemn_to_deadstemn_fire' 
     write(10,*) veg_nf%m_livestemn_to_deadstemn_fire
     write(10,*) 'veg_nf%m_npool_to_litter_fire' 
     write(10,*) veg_nf%m_npool_to_litter_fire
     write(10,*) 'veg_nf%m_leafn_to_fire' 
     write(10,*) veg_nf%m_leafn_to_fire
     write(10,*) 'veg_nf%m_frootn_xfer_to_litter_fire' 
     write(10,*) veg_nf%m_frootn_xfer_to_litter_fire
     write(10,*) 'veg_nf%m_livecrootn_to_fire' 
     write(10,*) veg_nf%m_livecrootn_to_fire
     write(10,*) 'veg_nf%m_leafn_storage_to_litter_fire' 
     write(10,*) veg_nf%m_leafn_storage_to_litter_fire
     write(10,*) 'veg_nf%m_retransn_to_litter_fire' 
     write(10,*) veg_nf%m_retransn_to_litter_fire
     write(10,*) 'veg_nf%m_livestemn_xfer_to_fire' 
     write(10,*) veg_nf%m_livestemn_xfer_to_fire
     write(10,*) 'veg_nf%m_livestemn_storage_to_fire' 
     write(10,*) veg_nf%m_livestemn_storage_to_fire
     write(10,*) 'veg_nf%m_deadcrootn_xfer_to_fire' 
     write(10,*) veg_nf%m_deadcrootn_xfer_to_fire
     write(10,*) 'veg_nf%m_frootn_xfer_to_fire' 
     write(10,*) veg_nf%m_frootn_xfer_to_fire
     write(10,*) 'veg_nf%m_retransn_to_fire' 
     write(10,*) veg_nf%m_retransn_to_fire
     write(10,*) 'veg_nf%m_livestemn_to_fire' 
     write(10,*) veg_nf%m_livestemn_to_fire
     write(10,*) 'veg_nf%m_deadstemn_xfer_to_fire' 
     write(10,*) veg_nf%m_deadstemn_xfer_to_fire
     write(10,*) 'veg_nf%m_livestemn_storage_to_litter_fire' 
     write(10,*) veg_nf%m_livestemn_storage_to_litter_fire
     write(10,*) 'veg_nf%m_leafn_storage_to_fire' 
     write(10,*) veg_nf%m_leafn_storage_to_fire
     write(10,*) 'veg_nf%m_leafn_to_litter_fire' 
     write(10,*) veg_nf%m_leafn_to_litter_fire
     write(10,*) 'veg_nf%m_livestemn_to_litter_fire' 
     write(10,*) veg_nf%m_livestemn_to_litter_fire
     write(10,*) 'veg_nf%m_deadstemn_storage_to_fire' 
     write(10,*) veg_nf%m_deadstemn_storage_to_fire
     write(10,*) 'veg_nf%m_livecrootn_xfer_to_litter_fire' 
     write(10,*) veg_nf%m_livecrootn_xfer_to_litter_fire
     write(10,*) 'veg_nf%m_leafn_xfer_to_litter_fire' 
     write(10,*) veg_nf%m_leafn_xfer_to_litter_fire
     write(10,*) 'veg_nf%m_npool_to_fire' 
     write(10,*) veg_nf%m_npool_to_fire
     write(10,*) 'veg_nf%m_frootn_to_fire' 
     write(10,*) veg_nf%m_frootn_to_fire
     write(10,*) 'veg_nf%m_livecrootn_xfer_to_fire' 
     write(10,*) veg_nf%m_livecrootn_xfer_to_fire
     write(10,*) 'veg_nf%m_leafn_xfer_to_fire' 
     write(10,*) veg_nf%m_leafn_xfer_to_fire
     write(10,*) 'veg_nf%m_deadstemn_to_fire' 
     write(10,*) veg_nf%m_deadstemn_to_fire
     write(10,*) 'veg_nf%m_deadstemn_xfer_to_litter_fire' 
     write(10,*) veg_nf%m_deadstemn_xfer_to_litter_fire
     write(10,*) 'veg_nf%m_deadcrootn_storage_to_litter_fire' 
     write(10,*) veg_nf%m_deadcrootn_storage_to_litter_fire
     write(10,*) 'veg_nf%m_deadstemn_to_litter_fire' 
     write(10,*) veg_nf%m_deadstemn_to_litter_fire
     write(10,*) 'veg_nf%m_deadcrootn_storage_to_fire' 
     write(10,*) veg_nf%m_deadcrootn_storage_to_fire
     write(10,*) 'veg_nf%m_frootn_storage_to_fire' 
     write(10,*) veg_nf%m_frootn_storage_to_fire
     write(10,*) 'veg_nf%m_livecrootn_to_litter_fire' 
     write(10,*) veg_nf%m_livecrootn_to_litter_fire
     write(10,*) 'veg_nf%m_livecrootn_storage_to_litter_fire' 
     write(10,*) veg_nf%m_livecrootn_storage_to_litter_fire
     write(10,*) 'veg_nf%m_deadcrootn_to_litter_fire' 
     write(10,*) veg_nf%m_deadcrootn_to_litter_fire
     write(10,*) 'veg_nf%m_frootn_to_litter_fire' 
     write(10,*) veg_nf%m_frootn_to_litter_fire
     write(10,*) 'veg_nf%m_livestemn_xfer_to_litter_fire' 
     write(10,*) veg_nf%m_livestemn_xfer_to_litter_fire
     write(10,*) 'veg_nf%m_livecrootn_to_deadcrootn_fire' 
     write(10,*) veg_nf%m_livecrootn_to_deadcrootn_fire
     write(10,*) 'veg_nf%m_frootn_storage_to_litter_fire' 
     write(10,*) veg_nf%m_frootn_storage_to_litter_fire
     write(10,*) 'veg_nf%m_deadcrootn_to_fire' 
     write(10,*) veg_nf%m_deadcrootn_to_fire
     write(10,*) 'veg_nf%m_deadstemn_storage_to_litter_fire' 
     write(10,*) veg_nf%m_deadstemn_storage_to_litter_fire
     write(10,*) 'veg_nf%m_livecrootn_storage_to_fire' 
     write(10,*) veg_nf%m_livecrootn_storage_to_fire
     write(10,*) 'veg_nf%m_deadcrootn_xfer_to_litter_fire' 
     write(10,*) veg_nf%m_deadcrootn_xfer_to_litter_fire
     write(10,*) 'veg_pf%sminp_to_plant' 
     write(10,*) veg_pf%sminp_to_plant
     write(10,*) 'veg_pf%plant_pdemand_vr' 
     write(10,*) veg_pf%plant_pdemand_vr
     write(10,*) 'veg_pf%ppool_to_grainp' 
     write(10,*) veg_pf%ppool_to_grainp
     write(10,*) 'veg_pf%ppool_to_livestemp_storage' 
     write(10,*) veg_pf%ppool_to_livestemp_storage
     write(10,*) 'veg_pf%ppool_to_leafp' 
     write(10,*) veg_pf%ppool_to_leafp
     write(10,*) 'veg_pf%plant_palloc' 
     write(10,*) veg_pf%plant_palloc
     write(10,*) 'veg_pf%ppool_to_livecrootp' 
     write(10,*) veg_pf%ppool_to_livecrootp
     write(10,*) 'veg_pf%ppool_to_deadcrootp_storage' 
     write(10,*) veg_pf%ppool_to_deadcrootp_storage
     write(10,*) 'veg_pf%avail_retransp' 
     write(10,*) veg_pf%avail_retransp
     write(10,*) 'veg_pf%ppool_to_frootp_storage' 
     write(10,*) veg_pf%ppool_to_frootp_storage
     write(10,*) 'veg_pf%ppool_to_deadstemp_storage' 
     write(10,*) veg_pf%ppool_to_deadstemp_storage
     write(10,*) 'veg_pf%retransp_to_ppool' 
     write(10,*) veg_pf%retransp_to_ppool
     write(10,*) 'veg_pf%ppool_to_deadstemp' 
     write(10,*) veg_pf%ppool_to_deadstemp
     write(10,*) 'veg_pf%ppool_to_grainp_storage' 
     write(10,*) veg_pf%ppool_to_grainp_storage
     write(10,*) 'veg_pf%ppool_to_leafp_storage' 
     write(10,*) veg_pf%ppool_to_leafp_storage
     write(10,*) 'veg_pf%supplement_to_plantp' 
     write(10,*) veg_pf%supplement_to_plantp
     write(10,*) 'veg_pf%ppool_to_deadcrootp' 
     write(10,*) veg_pf%ppool_to_deadcrootp
     write(10,*) 'veg_pf%sminp_to_ppool' 
     write(10,*) veg_pf%sminp_to_ppool
     write(10,*) 'veg_pf%ppool_to_livecrootp_storage' 
     write(10,*) veg_pf%ppool_to_livecrootp_storage
     write(10,*) 'veg_pf%ppool_to_frootp' 
     write(10,*) veg_pf%ppool_to_frootp
     write(10,*) 'veg_pf%ppool_to_livestemp' 
     write(10,*) veg_pf%ppool_to_livestemp
     write(10,*) 'veg_pf%frootp_xfer_to_frootp' 
     write(10,*) veg_pf%frootp_xfer_to_frootp
     write(10,*) 'veg_pf%deadstemp_xfer_to_deadstemp' 
     write(10,*) veg_pf%deadstemp_xfer_to_deadstemp
     write(10,*) 'veg_pf%leafp_xfer_to_leafp' 
     write(10,*) veg_pf%leafp_xfer_to_leafp
     write(10,*) 'veg_pf%deadstemp_storage_to_xfer' 
     write(10,*) veg_pf%deadstemp_storage_to_xfer
     write(10,*) 'veg_pf%livestemp_xfer_to_livestemp' 
     write(10,*) veg_pf%livestemp_xfer_to_livestemp
     write(10,*) 'veg_pf%deadcrootp_xfer_to_deadcrootp' 
     write(10,*) veg_pf%deadcrootp_xfer_to_deadcrootp
     write(10,*) 'veg_pf%livestemp_storage_to_xfer' 
     write(10,*) veg_pf%livestemp_storage_to_xfer
     write(10,*) 'veg_pf%leafp_storage_to_xfer' 
     write(10,*) veg_pf%leafp_storage_to_xfer
     write(10,*) 'veg_pf%frootp_storage_to_xfer' 
     write(10,*) veg_pf%frootp_storage_to_xfer
     write(10,*) 'veg_pf%livecrootp_xfer_to_livecrootp' 
     write(10,*) veg_pf%livecrootp_xfer_to_livecrootp
     write(10,*) 'veg_pf%deadcrootp_storage_to_xfer' 
     write(10,*) veg_pf%deadcrootp_storage_to_xfer
     write(10,*) 'veg_pf%livecrootp_storage_to_xfer' 
     write(10,*) veg_pf%livecrootp_storage_to_xfer
     write(10,*) 'veg_pf%frootp_xfer_to_frootp' 
     write(10,*) veg_pf%frootp_xfer_to_frootp
     write(10,*) 'veg_pf%deadstemp_xfer_to_deadstemp' 
     write(10,*) veg_pf%deadstemp_xfer_to_deadstemp
     write(10,*) 'veg_pf%leafp_xfer_to_leafp' 
     write(10,*) veg_pf%leafp_xfer_to_leafp
     write(10,*) 'veg_pf%deadstemp_storage_to_xfer' 
     write(10,*) veg_pf%deadstemp_storage_to_xfer
     write(10,*) 'veg_pf%livestemp_xfer_to_livestemp' 
     write(10,*) veg_pf%livestemp_xfer_to_livestemp
     write(10,*) 'veg_pf%deadcrootp_xfer_to_deadcrootp' 
     write(10,*) veg_pf%deadcrootp_xfer_to_deadcrootp
     write(10,*) 'veg_pf%livestemp_storage_to_xfer' 
     write(10,*) veg_pf%livestemp_storage_to_xfer
     write(10,*) 'veg_pf%leafp_storage_to_xfer' 
     write(10,*) veg_pf%leafp_storage_to_xfer
     write(10,*) 'veg_pf%frootp_storage_to_xfer' 
     write(10,*) veg_pf%frootp_storage_to_xfer
     write(10,*) 'veg_pf%livecrootp_xfer_to_livecrootp' 
     write(10,*) veg_pf%livecrootp_xfer_to_livecrootp
     write(10,*) 'veg_pf%deadcrootp_storage_to_xfer' 
     write(10,*) veg_pf%deadcrootp_storage_to_xfer
     write(10,*) 'veg_pf%livecrootp_storage_to_xfer' 
     write(10,*) veg_pf%livecrootp_storage_to_xfer
     write(10,*) 'veg_pf%crop_seedp_to_leaf' 
     write(10,*) veg_pf%crop_seedp_to_leaf
     write(10,*) 'veg_pf%frootp_xfer_to_frootp' 
     write(10,*) veg_pf%frootp_xfer_to_frootp
     write(10,*) 'veg_pf%deadstemp_xfer_to_deadstemp' 
     write(10,*) veg_pf%deadstemp_xfer_to_deadstemp
     write(10,*) 'veg_pf%leafp_xfer_to_leafp' 
     write(10,*) veg_pf%leafp_xfer_to_leafp
     write(10,*) 'veg_pf%livestemp_xfer_to_livestemp' 
     write(10,*) veg_pf%livestemp_xfer_to_livestemp
     write(10,*) 'veg_pf%deadcrootp_xfer_to_deadcrootp' 
     write(10,*) veg_pf%deadcrootp_xfer_to_deadcrootp
     write(10,*) 'veg_pf%livecrootp_xfer_to_livecrootp' 
     write(10,*) veg_pf%livecrootp_xfer_to_livecrootp
     write(10,*) 'veg_pf%hrv_livestemp_to_prod1p' 
     write(10,*) veg_pf%hrv_livestemp_to_prod1p
     write(10,*) 'veg_pf%hrv_grainp_to_prod1p' 
     write(10,*) veg_pf%hrv_grainp_to_prod1p
     write(10,*) 'veg_pf%hrv_leafp_to_prod1p' 
     write(10,*) veg_pf%hrv_leafp_to_prod1p
     write(10,*) 'veg_pf%hrv_cropp_to_prod1p' 
     write(10,*) veg_pf%hrv_cropp_to_prod1p
     write(10,*) 'veg_pf%leafp_to_litter' 
     write(10,*) veg_pf%leafp_to_litter
     write(10,*) 'veg_pf%frootp_to_litter' 
     write(10,*) veg_pf%frootp_to_litter
     write(10,*) 'veg_pf%livestemp_to_litter' 
     write(10,*) veg_pf%livestemp_to_litter
     write(10,*) 'veg_pf%leafp_to_retransp' 
     write(10,*) veg_pf%leafp_to_retransp
     write(10,*) 'veg_pf%leafp_to_litter' 
     write(10,*) veg_pf%leafp_to_litter
     write(10,*) 'veg_pf%frootp_to_litter' 
     write(10,*) veg_pf%frootp_to_litter
     write(10,*) 'veg_pf%leafp_to_retransp' 
     write(10,*) veg_pf%leafp_to_retransp
     write(10,*) 'veg_pf%livestemp_to_deadstemp' 
     write(10,*) veg_pf%livestemp_to_deadstemp
     write(10,*) 'veg_pf%livecrootp_to_deadcrootp' 
     write(10,*) veg_pf%livecrootp_to_deadcrootp
     write(10,*) 'veg_pf%livestemp_to_retransp' 
     write(10,*) veg_pf%livestemp_to_retransp
     write(10,*) 'veg_pf%livecrootp_to_retransp' 
     write(10,*) veg_pf%livecrootp_to_retransp
     write(10,*) 'veg_pf%ppool_to_livestemp_storage' 
     write(10,*) veg_pf%ppool_to_livestemp_storage
     write(10,*) 'veg_pf%deadstemp_storage_to_xfer' 
     write(10,*) veg_pf%deadstemp_storage_to_xfer
     write(10,*) 'veg_pf%livestemp_to_deadstemp' 
     write(10,*) veg_pf%livestemp_to_deadstemp
     write(10,*) 'veg_pf%retransp_to_ppool' 
     write(10,*) veg_pf%retransp_to_ppool
     write(10,*) 'veg_pf%livecrootp_to_retransp' 
     write(10,*) veg_pf%livecrootp_to_retransp
     write(10,*) 'veg_pf%livecrootp_xfer_to_livecrootp' 
     write(10,*) veg_pf%livecrootp_xfer_to_livecrootp
     write(10,*) 'veg_pf%ppool_to_grainp' 
     write(10,*) veg_pf%ppool_to_grainp
     write(10,*) 'veg_pf%livestemp_xfer_to_livestemp' 
     write(10,*) veg_pf%livestemp_xfer_to_livestemp
     write(10,*) 'veg_pf%ppool_to_deadcrootp_storage' 
     write(10,*) veg_pf%ppool_to_deadcrootp_storage
     write(10,*) 'veg_pf%ppool_to_grainp_storage' 
     write(10,*) veg_pf%ppool_to_grainp_storage
     write(10,*) 'veg_pf%ppool_to_leafp_storage' 
     write(10,*) veg_pf%ppool_to_leafp_storage
     write(10,*) 'veg_pf%leafp_to_retransp' 
     write(10,*) veg_pf%leafp_to_retransp
     write(10,*) 'veg_pf%ppool_to_deadcrootp' 
     write(10,*) veg_pf%ppool_to_deadcrootp
     write(10,*) 'veg_pf%livestemp_to_litter' 
     write(10,*) veg_pf%livestemp_to_litter
     write(10,*) 'veg_pf%livecrootp_to_deadcrootp' 
     write(10,*) veg_pf%livecrootp_to_deadcrootp
     write(10,*) 'veg_pf%deadstemp_xfer_to_deadstemp' 
     write(10,*) veg_pf%deadstemp_xfer_to_deadstemp
     write(10,*) 'veg_pf%frootp_to_litter' 
     write(10,*) veg_pf%frootp_to_litter
     write(10,*) 'veg_pf%ppool_to_frootp_storage' 
     write(10,*) veg_pf%ppool_to_frootp_storage
     write(10,*) 'veg_pf%grainp_xfer_to_grainp' 
     write(10,*) veg_pf%grainp_xfer_to_grainp
     write(10,*) 'veg_pf%frootp_to_retransp' 
     write(10,*) veg_pf%frootp_to_retransp
     write(10,*) 'veg_pf%livestemp_to_retransp' 
     write(10,*) veg_pf%livestemp_to_retransp
     write(10,*) 'veg_pf%supplement_to_plantp' 
     write(10,*) veg_pf%supplement_to_plantp
     write(10,*) 'veg_pf%livecrootp_storage_to_xfer' 
     write(10,*) veg_pf%livecrootp_storage_to_xfer
     write(10,*) 'veg_pf%ppool_to_livestemp' 
     write(10,*) veg_pf%ppool_to_livestemp
     write(10,*) 'veg_pf%frootp_xfer_to_frootp' 
     write(10,*) veg_pf%frootp_xfer_to_frootp
     write(10,*) 'veg_pf%leafp_xfer_to_leafp' 
     write(10,*) veg_pf%leafp_xfer_to_leafp
     write(10,*) 'veg_pf%ppool_to_leafp' 
     write(10,*) veg_pf%ppool_to_leafp
     write(10,*) 'veg_pf%ppool_to_livecrootp' 
     write(10,*) veg_pf%ppool_to_livecrootp
     write(10,*) 'veg_pf%grainp_to_food' 
     write(10,*) veg_pf%grainp_to_food
     write(10,*) 'veg_pf%deadcrootp_xfer_to_deadcrootp' 
     write(10,*) veg_pf%deadcrootp_xfer_to_deadcrootp
     write(10,*) 'veg_pf%livestemp_storage_to_xfer' 
     write(10,*) veg_pf%livestemp_storage_to_xfer
     write(10,*) 'veg_pf%ppool_to_deadstemp_storage' 
     write(10,*) veg_pf%ppool_to_deadstemp_storage
     write(10,*) 'veg_pf%leafp_to_litter' 
     write(10,*) veg_pf%leafp_to_litter
     write(10,*) 'veg_pf%ppool_to_deadstemp' 
     write(10,*) veg_pf%ppool_to_deadstemp
     write(10,*) 'veg_pf%frootp_storage_to_xfer' 
     write(10,*) veg_pf%frootp_storage_to_xfer
     write(10,*) 'veg_pf%ppool_to_livecrootp_storage' 
     write(10,*) veg_pf%ppool_to_livecrootp_storage
     write(10,*) 'veg_pf%ppool_to_frootp' 
     write(10,*) veg_pf%ppool_to_frootp
     write(10,*) 'veg_pf%deadcrootp_storage_to_xfer' 
     write(10,*) veg_pf%deadcrootp_storage_to_xfer
     write(10,*) 'veg_pf%m_frootp_storage_to_litter' 
     write(10,*) veg_pf%m_frootp_storage_to_litter
     write(10,*) 'veg_pf%m_frootp_to_litter' 
     write(10,*) veg_pf%m_frootp_to_litter
     write(10,*) 'veg_pf%m_deadstemp_to_litter' 
     write(10,*) veg_pf%m_deadstemp_to_litter
     write(10,*) 'veg_pf%m_livecrootp_storage_to_litter' 
     write(10,*) veg_pf%m_livecrootp_storage_to_litter
     write(10,*) 'veg_pf%m_livestemp_xfer_to_litter' 
     write(10,*) veg_pf%m_livestemp_xfer_to_litter
     write(10,*) 'veg_pf%m_leafp_storage_to_litter' 
     write(10,*) veg_pf%m_leafp_storage_to_litter
     write(10,*) 'veg_pf%m_deadstemp_xfer_to_litter' 
     write(10,*) veg_pf%m_deadstemp_xfer_to_litter
     write(10,*) 'veg_pf%m_livestemp_storage_to_litter' 
     write(10,*) veg_pf%m_livestemp_storage_to_litter
     write(10,*) 'veg_pf%m_livecrootp_to_litter' 
     write(10,*) veg_pf%m_livecrootp_to_litter
     write(10,*) 'veg_pf%m_deadstemp_storage_to_litter' 
     write(10,*) veg_pf%m_deadstemp_storage_to_litter
     write(10,*) 'veg_pf%m_livestemp_to_litter' 
     write(10,*) veg_pf%m_livestemp_to_litter
     write(10,*) 'veg_pf%m_livecrootp_xfer_to_litter' 
     write(10,*) veg_pf%m_livecrootp_xfer_to_litter
     write(10,*) 'veg_pf%m_deadcrootp_storage_to_litter' 
     write(10,*) veg_pf%m_deadcrootp_storage_to_litter
     write(10,*) 'veg_pf%m_leafp_to_litter' 
     write(10,*) veg_pf%m_leafp_to_litter
     write(10,*) 'veg_pf%m_retransp_to_litter' 
     write(10,*) veg_pf%m_retransp_to_litter
     write(10,*) 'veg_pf%m_deadcrootp_xfer_to_litter' 
     write(10,*) veg_pf%m_deadcrootp_xfer_to_litter
     write(10,*) 'veg_pf%m_deadcrootp_to_litter' 
     write(10,*) veg_pf%m_deadcrootp_to_litter
     write(10,*) 'veg_pf%m_leafp_xfer_to_litter' 
     write(10,*) veg_pf%m_leafp_xfer_to_litter
     write(10,*) 'veg_pf%m_ppool_to_litter' 
     write(10,*) veg_pf%m_ppool_to_litter
     write(10,*) 'veg_pf%m_frootp_xfer_to_litter' 
     write(10,*) veg_pf%m_frootp_xfer_to_litter
     write(10,*) 'veg_pf%hrv_leafp_to_litter' 
     write(10,*) veg_pf%hrv_leafp_to_litter
     write(10,*) 'veg_pf%hrv_retransp_to_litter' 
     write(10,*) veg_pf%hrv_retransp_to_litter
     write(10,*) 'veg_pf%hrv_deadcrootp_to_litter' 
     write(10,*) veg_pf%hrv_deadcrootp_to_litter
     write(10,*) 'veg_pf%hrv_leafp_xfer_to_litter' 
     write(10,*) veg_pf%hrv_leafp_xfer_to_litter
     write(10,*) 'veg_pf%hrv_ppool_to_litter' 
     write(10,*) veg_pf%hrv_ppool_to_litter
     write(10,*) 'veg_pf%hrv_deadstemp_storage_to_litter' 
     write(10,*) veg_pf%hrv_deadstemp_storage_to_litter
     write(10,*) 'veg_pf%hrv_livestemp_storage_to_litter' 
     write(10,*) veg_pf%hrv_livestemp_storage_to_litter
     write(10,*) 'veg_pf%hrv_deadcrootp_xfer_to_litter' 
     write(10,*) veg_pf%hrv_deadcrootp_xfer_to_litter
     write(10,*) 'veg_pf%hrv_frootp_storage_to_litter' 
     write(10,*) veg_pf%hrv_frootp_storage_to_litter
     write(10,*) 'veg_pf%hrv_deadcrootp_storage_to_litter' 
     write(10,*) veg_pf%hrv_deadcrootp_storage_to_litter
     write(10,*) 'veg_pf%hrv_deadstemp_to_prod10p' 
     write(10,*) veg_pf%hrv_deadstemp_to_prod10p
     write(10,*) 'veg_pf%hrv_deadstemp_xfer_to_litter' 
     write(10,*) veg_pf%hrv_deadstemp_xfer_to_litter
     write(10,*) 'veg_pf%hrv_livestemp_xfer_to_litter' 
     write(10,*) veg_pf%hrv_livestemp_xfer_to_litter
     write(10,*) 'veg_pf%hrv_livecrootp_to_litter' 
     write(10,*) veg_pf%hrv_livecrootp_to_litter
     write(10,*) 'veg_pf%hrv_frootp_to_litter' 
     write(10,*) veg_pf%hrv_frootp_to_litter
     write(10,*) 'veg_pf%hrv_livecrootp_storage_to_litter' 
     write(10,*) veg_pf%hrv_livecrootp_storage_to_litter
     write(10,*) 'veg_pf%hrv_deadstemp_to_prod100p' 
     write(10,*) veg_pf%hrv_deadstemp_to_prod100p
     write(10,*) 'veg_pf%hrv_leafp_storage_to_litter' 
     write(10,*) veg_pf%hrv_leafp_storage_to_litter
     write(10,*) 'veg_pf%hrv_livestemp_to_litter' 
     write(10,*) veg_pf%hrv_livestemp_to_litter
     write(10,*) 'veg_pf%hrv_frootp_xfer_to_litter' 
     write(10,*) veg_pf%hrv_frootp_xfer_to_litter
     write(10,*) 'veg_pf%hrv_livecrootp_xfer_to_litter' 
     write(10,*) veg_pf%hrv_livecrootp_xfer_to_litter
     write(10,*) 'veg_pf%m_leafp_to_fire' 
     write(10,*) veg_pf%m_leafp_to_fire
     write(10,*) 'veg_pf%m_deadstemp_xfer_to_litter_fire' 
     write(10,*) veg_pf%m_deadstemp_xfer_to_litter_fire
     write(10,*) 'veg_pf%m_livestemp_storage_to_litter_fire' 
     write(10,*) veg_pf%m_livestemp_storage_to_litter_fire
     write(10,*) 'veg_pf%m_livestemp_xfer_to_litter_fire' 
     write(10,*) veg_pf%m_livestemp_xfer_to_litter_fire
     write(10,*) 'veg_pf%m_frootp_xfer_to_litter_fire' 
     write(10,*) veg_pf%m_frootp_xfer_to_litter_fire
     write(10,*) 'veg_pf%m_retransp_to_litter_fire' 
     write(10,*) veg_pf%m_retransp_to_litter_fire
     write(10,*) 'veg_pf%m_deadcrootp_xfer_to_litter_fire' 
     write(10,*) veg_pf%m_deadcrootp_xfer_to_litter_fire
     write(10,*) 'veg_pf%m_deadcrootp_xfer_to_fire' 
     write(10,*) veg_pf%m_deadcrootp_xfer_to_fire
     write(10,*) 'veg_pf%m_deadcrootp_to_fire' 
     write(10,*) veg_pf%m_deadcrootp_to_fire
     write(10,*) 'veg_pf%m_livecrootp_storage_to_litter_fire' 
     write(10,*) veg_pf%m_livecrootp_storage_to_litter_fire
     write(10,*) 'veg_pf%m_livecrootp_storage_to_fire' 
     write(10,*) veg_pf%m_livecrootp_storage_to_fire
     write(10,*) 'veg_pf%m_livestemp_to_deadstemp_fire' 
     write(10,*) veg_pf%m_livestemp_to_deadstemp_fire
     write(10,*) 'veg_pf%m_frootp_to_fire' 
     write(10,*) veg_pf%m_frootp_to_fire
     write(10,*) 'veg_pf%m_leafp_to_litter_fire' 
     write(10,*) veg_pf%m_leafp_to_litter_fire
     write(10,*) 'veg_pf%m_livecrootp_to_deadcrootp_fire' 
     write(10,*) veg_pf%m_livecrootp_to_deadcrootp_fire
     write(10,*) 'veg_pf%m_deadstemp_to_fire' 
     write(10,*) veg_pf%m_deadstemp_to_fire
     write(10,*) 'veg_pf%m_deadstemp_storage_to_fire' 
     write(10,*) veg_pf%m_deadstemp_storage_to_fire
     write(10,*) 'veg_pf%m_ppool_to_fire' 
     write(10,*) veg_pf%m_ppool_to_fire
     write(10,*) 'veg_pf%m_livestemp_to_fire' 
     write(10,*) veg_pf%m_livestemp_to_fire
     write(10,*) 'veg_pf%m_frootp_storage_to_fire' 
     write(10,*) veg_pf%m_frootp_storage_to_fire
     write(10,*) 'veg_pf%m_livecrootp_to_litter_fire' 
     write(10,*) veg_pf%m_livecrootp_to_litter_fire
     write(10,*) 'veg_pf%m_frootp_storage_to_litter_fire' 
     write(10,*) veg_pf%m_frootp_storage_to_litter_fire
     write(10,*) 'veg_pf%m_livecrootp_xfer_to_fire' 
     write(10,*) veg_pf%m_livecrootp_xfer_to_fire
     write(10,*) 'veg_pf%m_ppool_to_litter_fire' 
     write(10,*) veg_pf%m_ppool_to_litter_fire
     write(10,*) 'veg_pf%m_deadstemp_storage_to_litter_fire' 
     write(10,*) veg_pf%m_deadstemp_storage_to_litter_fire
     write(10,*) 'veg_pf%m_frootp_xfer_to_fire' 
     write(10,*) veg_pf%m_frootp_xfer_to_fire
     write(10,*) 'veg_pf%m_livestemp_storage_to_fire' 
     write(10,*) veg_pf%m_livestemp_storage_to_fire
     write(10,*) 'veg_pf%m_livestemp_to_litter_fire' 
     write(10,*) veg_pf%m_livestemp_to_litter_fire
     write(10,*) 'veg_pf%m_deadstemp_xfer_to_fire' 
     write(10,*) veg_pf%m_deadstemp_xfer_to_fire
     write(10,*) 'veg_pf%m_retransp_to_fire' 
     write(10,*) veg_pf%m_retransp_to_fire
     write(10,*) 'veg_pf%m_livecrootp_to_fire' 
     write(10,*) veg_pf%m_livecrootp_to_fire
     write(10,*) 'veg_pf%m_leafp_storage_to_litter_fire' 
     write(10,*) veg_pf%m_leafp_storage_to_litter_fire
     write(10,*) 'veg_pf%m_leafp_xfer_to_fire' 
     write(10,*) veg_pf%m_leafp_xfer_to_fire
     write(10,*) 'veg_pf%m_frootp_to_litter_fire' 
     write(10,*) veg_pf%m_frootp_to_litter_fire
     write(10,*) 'veg_pf%m_leafp_storage_to_fire' 
     write(10,*) veg_pf%m_leafp_storage_to_fire
     write(10,*) 'veg_pf%m_deadcrootp_storage_to_fire' 
     write(10,*) veg_pf%m_deadcrootp_storage_to_fire
     write(10,*) 'veg_pf%m_deadcrootp_storage_to_litter_fire' 
     write(10,*) veg_pf%m_deadcrootp_storage_to_litter_fire
     write(10,*) 'veg_pf%m_livecrootp_xfer_to_litter_fire' 
     write(10,*) veg_pf%m_livecrootp_xfer_to_litter_fire
     write(10,*) 'veg_pf%m_deadstemp_to_litter_fire' 
     write(10,*) veg_pf%m_deadstemp_to_litter_fire
     write(10,*) 'veg_pf%m_deadcrootp_to_litter_fire' 
     write(10,*) veg_pf%m_deadcrootp_to_litter_fire
     write(10,*) 'veg_pf%m_leafp_xfer_to_litter_fire' 
     write(10,*) veg_pf%m_leafp_xfer_to_litter_fire
     write(10,*) 'veg_pf%m_livestemp_xfer_to_fire' 
     write(10,*) veg_pf%m_livestemp_xfer_to_fire
     write(10,*) 'veg_ns%pnup_pfrootc' 
     write(10,*) veg_ns%pnup_pfrootc
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn' 
     write(10,*) veg_ns%livestemn
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ns%grainn_storage' 
     write(10,*) veg_ns%grainn_storage
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%retransn' 
     write(10,*) veg_ns%retransn
     write(10,*) 'veg_ns%npool' 
     write(10,*) veg_ns%npool
     write(10,*) 'veg_ns%cropseedn_deficit' 
     write(10,*) veg_ns%cropseedn_deficit
     write(10,*) 'veg_ns%grainn_xfer' 
     write(10,*) veg_ns%grainn_xfer
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn_storage' 
     write(10,*) veg_ns%livestemn_storage
     write(10,*) 'veg_ns%frootn' 
     write(10,*) veg_ns%frootn
     write(10,*) 'veg_ns%livecrootn' 
     write(10,*) veg_ns%livecrootn
     write(10,*) 'veg_ns%leafn' 
     write(10,*) veg_ns%leafn
     write(10,*) 'veg_ns%livecrootn_storage' 
     write(10,*) veg_ns%livecrootn_storage
     write(10,*) 'veg_ns%frootn_storage' 
     write(10,*) veg_ns%frootn_storage
     write(10,*) 'veg_ns%leafn_storage' 
     write(10,*) veg_ns%leafn_storage
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%deadstemn_storage' 
     write(10,*) veg_ns%deadstemn_storage
     write(10,*) 'veg_ns%deadcrootn_storage' 
     write(10,*) veg_ns%deadcrootn_storage
     write(10,*) 'veg_ns%deadcrootn' 
     write(10,*) veg_ns%deadcrootn
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%grainn' 
     write(10,*) veg_ns%grainn
     write(10,*) 'veg_ns%livestemn' 
     write(10,*) veg_ns%livestemn
     write(10,*) 'veg_ns%deadcrootn_storage' 
     write(10,*) veg_ns%deadcrootn_storage
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ns%npool' 
     write(10,*) veg_ns%npool
     write(10,*) 'veg_ns%frootn' 
     write(10,*) veg_ns%frootn
     write(10,*) 'veg_ns%deadstemn_storage' 
     write(10,*) veg_ns%deadstemn_storage
     write(10,*) 'veg_ns%livecrootn' 
     write(10,*) veg_ns%livecrootn
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%livestemn_storage' 
     write(10,*) veg_ns%livestemn_storage
     write(10,*) 'veg_ns%deadcrootn' 
     write(10,*) veg_ns%deadcrootn
     write(10,*) 'veg_ns%livecrootn_storage' 
     write(10,*) veg_ns%livecrootn_storage
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%frootn_storage' 
     write(10,*) veg_ns%frootn_storage
     write(10,*) 'veg_ns%leafn_storage' 
     write(10,*) veg_ns%leafn_storage
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%retransn' 
     write(10,*) veg_ns%retransn
     write(10,*) 'veg_ns%leafn' 
     write(10,*) veg_ns%leafn
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn' 
     write(10,*) veg_ns%livestemn
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%retransn' 
     write(10,*) veg_ns%retransn
     write(10,*) 'veg_ns%npool' 
     write(10,*) veg_ns%npool
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn_storage' 
     write(10,*) veg_ns%livestemn_storage
     write(10,*) 'veg_ns%frootn' 
     write(10,*) veg_ns%frootn
     write(10,*) 'veg_ns%leafn' 
     write(10,*) veg_ns%leafn
     write(10,*) 'veg_ns%livecrootn_storage' 
     write(10,*) veg_ns%livecrootn_storage
     write(10,*) 'veg_ns%frootn_storage' 
     write(10,*) veg_ns%frootn_storage
     write(10,*) 'veg_ns%leafn_storage' 
     write(10,*) veg_ns%leafn_storage
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%deadstemn_storage' 
     write(10,*) veg_ns%deadstemn_storage
     write(10,*) 'veg_ns%grainn' 
     write(10,*) veg_ns%grainn
     write(10,*) 'veg_ns%deadcrootn_storage' 
     write(10,*) veg_ns%deadcrootn_storage
     write(10,*) 'veg_ns%deadcrootn' 
     write(10,*) veg_ns%deadcrootn
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%livecrootn' 
     write(10,*) veg_ns%livecrootn
     write(10,*) 'col_ns%smin_no3_vr' 
     write(10,*) col_ns%smin_no3_vr
     write(10,*) 'col_ns%smin_nh4_vr' 
     write(10,*) col_ns%smin_nh4_vr
     write(10,*) 'col_ns%smin_nh4_vr' 
     write(10,*) col_ns%smin_nh4_vr
     write(10,*) 'col_ns%smin_no3_vr' 
     write(10,*) col_ns%smin_no3_vr
     write(10,*) 'col_ns%sminn_vr' 
     write(10,*) col_ns%sminn_vr
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'col_ns%prod100n' 
     write(10,*) col_ns%prod100n
     write(10,*) 'col_ns%prod10n' 
     write(10,*) col_ns%prod10n
     write(10,*) 'col_ns%prod1n' 
     write(10,*) col_ns%prod1n
     write(10,*) 'veg_pp%itype' 
     write(10,*) veg_pp%itype
     write(10,*) 'veg_pp%itype' 
     write(10,*) veg_pp%itype
     write(10,*) 'veg_vp%woody' 
     write(10,*) veg_vp%woody
     write(10,*) 'veg_cf%cpool_to_leafc_storage' 
     write(10,*) veg_cf%cpool_to_leafc_storage
     write(10,*) 'veg_cf%xsmrpool_turnover' 
     write(10,*) veg_cf%xsmrpool_turnover
     write(10,*) 'veg_cf%grain_curmr' 
     write(10,*) veg_cf%grain_curmr
     write(10,*) 'veg_cf%froot_mr' 
     write(10,*) veg_cf%froot_mr
     write(10,*) 'veg_cf%cpool_to_deadcrootc_storage' 
     write(10,*) veg_cf%cpool_to_deadcrootc_storage
     write(10,*) 'veg_cf%livestem_curmr' 
     write(10,*) veg_cf%livestem_curmr
     write(10,*) 'veg_cf%grain_mr' 
     write(10,*) veg_cf%grain_mr
     write(10,*) 'veg_cf%livestem_xsmr' 
     write(10,*) veg_cf%livestem_xsmr
     write(10,*) 'veg_cf%allocation_stem' 
     write(10,*) veg_cf%allocation_stem
     write(10,*) 'veg_cf%cpool_to_grainc_storage' 
     write(10,*) veg_cf%cpool_to_grainc_storage
     write(10,*) 'veg_cf%froot_curmr' 
     write(10,*) veg_cf%froot_curmr
     write(10,*) 'veg_cf%availc' 
     write(10,*) veg_cf%availc
     write(10,*) 'veg_cf%livecroot_curmr' 
     write(10,*) veg_cf%livecroot_curmr
     write(10,*) 'veg_cf%livecroot_mr' 
     write(10,*) veg_cf%livecroot_mr
     write(10,*) 'veg_cf%plant_calloc' 
     write(10,*) veg_cf%plant_calloc
     write(10,*) 'veg_cf%cpool_to_grainc' 
     write(10,*) veg_cf%cpool_to_grainc
     write(10,*) 'veg_cf%livestem_mr' 
     write(10,*) veg_cf%livestem_mr
     write(10,*) 'veg_cf%leaf_curmr' 
     write(10,*) veg_cf%leaf_curmr
     write(10,*) 'veg_cf%cpool_to_gresp_storage' 
     write(10,*) veg_cf%cpool_to_gresp_storage
     write(10,*) 'veg_cf%cpool_to_livestemc_storage' 
     write(10,*) veg_cf%cpool_to_livestemc_storage
     write(10,*) 'veg_cf%grain_xsmr' 
     write(10,*) veg_cf%grain_xsmr
     write(10,*) 'veg_cf%cpool_to_deadstemc_storage' 
     write(10,*) veg_cf%cpool_to_deadstemc_storage
     write(10,*) 'veg_cf%froot_xsmr' 
     write(10,*) veg_cf%froot_xsmr
     write(10,*) 'veg_cf%cpool_to_frootc_storage' 
     write(10,*) veg_cf%cpool_to_frootc_storage
     write(10,*) 'veg_cf%excess_cflux' 
     write(10,*) veg_cf%excess_cflux
     write(10,*) 'veg_cf%psnshade_to_cpool' 
     write(10,*) veg_cf%psnshade_to_cpool
     write(10,*) 'veg_cf%psnsun_to_cpool' 
     write(10,*) veg_cf%psnsun_to_cpool
     write(10,*) 'veg_cf%cpool_to_leafc' 
     write(10,*) veg_cf%cpool_to_leafc
     write(10,*) 'veg_cf%cpool_to_livecrootc' 
     write(10,*) veg_cf%cpool_to_livecrootc
     write(10,*) 'veg_cf%livecroot_xsmr' 
     write(10,*) veg_cf%livecroot_xsmr
     write(10,*) 'veg_cf%allocation_froot' 
     write(10,*) veg_cf%allocation_froot
     write(10,*) 'veg_cf%leaf_mr' 
     write(10,*) veg_cf%leaf_mr
     write(10,*) 'veg_cf%cpool_to_frootc' 
     write(10,*) veg_cf%cpool_to_frootc
     write(10,*) 'veg_cf%xsmrpool_recover' 
     write(10,*) veg_cf%xsmrpool_recover
     write(10,*) 'veg_cf%allocation_leaf' 
     write(10,*) veg_cf%allocation_leaf
     write(10,*) 'veg_cf%cpool_to_deadstemc' 
     write(10,*) veg_cf%cpool_to_deadstemc
     write(10,*) 'veg_cf%cpool_to_xsmrpool' 
     write(10,*) veg_cf%cpool_to_xsmrpool
     write(10,*) 'veg_cf%cpool_to_deadcrootc' 
     write(10,*) veg_cf%cpool_to_deadcrootc
     write(10,*) 'veg_cf%leaf_xsmr' 
     write(10,*) veg_cf%leaf_xsmr
     write(10,*) 'veg_cf%cpool_to_livecrootc_storage' 
     write(10,*) veg_cf%cpool_to_livecrootc_storage
     write(10,*) 'veg_cf%cpool_to_livestemc' 
     write(10,*) veg_cf%cpool_to_livestemc
     write(10,*) 'veg_cf%prev_frootc_to_litter' 
     write(10,*) veg_cf%prev_frootc_to_litter
     write(10,*) 'veg_cf%frootc_storage_to_xfer' 
     write(10,*) veg_cf%frootc_storage_to_xfer
     write(10,*) 'veg_cf%leafc_xfer_to_leafc' 
     write(10,*) veg_cf%leafc_xfer_to_leafc
     write(10,*) 'veg_cf%livecrootc_storage_to_xfer' 
     write(10,*) veg_cf%livecrootc_storage_to_xfer
     write(10,*) 'veg_cf%livestemc_xfer_to_livestemc' 
     write(10,*) veg_cf%livestemc_xfer_to_livestemc
     write(10,*) 'veg_cf%livestemc_storage_to_xfer' 
     write(10,*) veg_cf%livestemc_storage_to_xfer
     write(10,*) 'veg_cf%deadstemc_xfer_to_deadstemc' 
     write(10,*) veg_cf%deadstemc_xfer_to_deadstemc
     write(10,*) 'veg_cf%deadstemc_storage_to_xfer' 
     write(10,*) veg_cf%deadstemc_storage_to_xfer
     write(10,*) 'veg_cf%leafc_storage_to_xfer' 
     write(10,*) veg_cf%leafc_storage_to_xfer
     write(10,*) 'veg_cf%gresp_storage_to_xfer' 
     write(10,*) veg_cf%gresp_storage_to_xfer
     write(10,*) 'veg_cf%prev_leafc_to_litter' 
     write(10,*) veg_cf%prev_leafc_to_litter
     write(10,*) 'veg_cf%deadcrootc_xfer_to_deadcrootc' 
     write(10,*) veg_cf%deadcrootc_xfer_to_deadcrootc
     write(10,*) 'veg_cf%livecrootc_xfer_to_livecrootc' 
     write(10,*) veg_cf%livecrootc_xfer_to_livecrootc
     write(10,*) 'veg_cf%frootc_xfer_to_frootc' 
     write(10,*) veg_cf%frootc_xfer_to_frootc
     write(10,*) 'veg_cf%deadcrootc_storage_to_xfer' 
     write(10,*) veg_cf%deadcrootc_storage_to_xfer
     write(10,*) 'veg_cf%prev_frootc_to_litter' 
     write(10,*) veg_cf%prev_frootc_to_litter
     write(10,*) 'veg_cf%frootc_storage_to_xfer' 
     write(10,*) veg_cf%frootc_storage_to_xfer
     write(10,*) 'veg_cf%leafc_xfer_to_leafc' 
     write(10,*) veg_cf%leafc_xfer_to_leafc
     write(10,*) 'veg_cf%livecrootc_storage_to_xfer' 
     write(10,*) veg_cf%livecrootc_storage_to_xfer
     write(10,*) 'veg_cf%livestemc_xfer_to_livestemc' 
     write(10,*) veg_cf%livestemc_xfer_to_livestemc
     write(10,*) 'veg_cf%livestemc_storage_to_xfer' 
     write(10,*) veg_cf%livestemc_storage_to_xfer
     write(10,*) 'veg_cf%deadstemc_xfer_to_deadstemc' 
     write(10,*) veg_cf%deadstemc_xfer_to_deadstemc
     write(10,*) 'veg_cf%deadstemc_storage_to_xfer' 
     write(10,*) veg_cf%deadstemc_storage_to_xfer
     write(10,*) 'veg_cf%leafc_storage_to_xfer' 
     write(10,*) veg_cf%leafc_storage_to_xfer
     write(10,*) 'veg_cf%gresp_storage_to_xfer' 
     write(10,*) veg_cf%gresp_storage_to_xfer
     write(10,*) 'veg_cf%prev_leafc_to_litter' 
     write(10,*) veg_cf%prev_leafc_to_litter
     write(10,*) 'veg_cf%deadcrootc_xfer_to_deadcrootc' 
     write(10,*) veg_cf%deadcrootc_xfer_to_deadcrootc
     write(10,*) 'veg_cf%livecrootc_xfer_to_livecrootc' 
     write(10,*) veg_cf%livecrootc_xfer_to_livecrootc
     write(10,*) 'veg_cf%frootc_xfer_to_frootc' 
     write(10,*) veg_cf%frootc_xfer_to_frootc
     write(10,*) 'veg_cf%deadcrootc_storage_to_xfer' 
     write(10,*) veg_cf%deadcrootc_storage_to_xfer
     write(10,*) 'veg_cf%crop_seedc_to_leaf' 
     write(10,*) veg_cf%crop_seedc_to_leaf
     write(10,*) 'veg_cf%leafc_xfer_to_leafc' 
     write(10,*) veg_cf%leafc_xfer_to_leafc
     write(10,*) 'veg_cf%livestemc_xfer_to_livestemc' 
     write(10,*) veg_cf%livestemc_xfer_to_livestemc
     write(10,*) 'veg_cf%deadstemc_xfer_to_deadstemc' 
     write(10,*) veg_cf%deadstemc_xfer_to_deadstemc
     write(10,*) 'veg_cf%deadcrootc_xfer_to_deadcrootc' 
     write(10,*) veg_cf%deadcrootc_xfer_to_deadcrootc
     write(10,*) 'veg_cf%livecrootc_xfer_to_livecrootc' 
     write(10,*) veg_cf%livecrootc_xfer_to_livecrootc
     write(10,*) 'veg_cf%frootc_xfer_to_frootc' 
     write(10,*) veg_cf%frootc_xfer_to_frootc
     write(10,*) 'veg_cf%hrv_livestemc_to_prod1c' 
     write(10,*) veg_cf%hrv_livestemc_to_prod1c
     write(10,*) 'veg_cf%hrv_leafc_to_prod1c' 
     write(10,*) veg_cf%hrv_leafc_to_prod1c
     write(10,*) 'veg_cf%hrv_grainc_to_prod1c' 
     write(10,*) veg_cf%hrv_grainc_to_prod1c
     write(10,*) 'veg_cf%hrv_cropc_to_prod1c' 
     write(10,*) veg_cf%hrv_cropc_to_prod1c
     write(10,*) 'veg_cf%prev_frootc_to_litter' 
     write(10,*) veg_cf%prev_frootc_to_litter
     write(10,*) 'veg_cf%frootc_to_litter' 
     write(10,*) veg_cf%frootc_to_litter
     write(10,*) 'veg_cf%leafc_to_litter' 
     write(10,*) veg_cf%leafc_to_litter
     write(10,*) 'veg_cf%livestemc_to_litter' 
     write(10,*) veg_cf%livestemc_to_litter
     write(10,*) 'veg_cf%prev_leafc_to_litter' 
     write(10,*) veg_cf%prev_leafc_to_litter
     write(10,*) 'veg_cf%frootc_to_litter' 
     write(10,*) veg_cf%frootc_to_litter
     write(10,*) 'veg_cf%leafc_to_litter' 
     write(10,*) veg_cf%leafc_to_litter
     write(10,*) 'veg_cf%livestemc_to_deadstemc' 
     write(10,*) veg_cf%livestemc_to_deadstemc
     write(10,*) 'veg_cf%livecrootc_to_deadcrootc' 
     write(10,*) veg_cf%livecrootc_to_deadcrootc
     write(10,*) 'veg_cf%cpool_livecroot_storage_gr' 
     write(10,*) veg_cf%cpool_livecroot_storage_gr
     write(10,*) 'veg_cf%transfer_deadstem_gr' 
     write(10,*) veg_cf%transfer_deadstem_gr
     write(10,*) 'veg_cf%transfer_livecroot_gr' 
     write(10,*) veg_cf%transfer_livecroot_gr
     write(10,*) 'veg_cf%cpool_grain_storage_gr' 
     write(10,*) veg_cf%cpool_grain_storage_gr
     write(10,*) 'veg_cf%cpool_leaf_storage_gr' 
     write(10,*) veg_cf%cpool_leaf_storage_gr
     write(10,*) 'veg_cf%transfer_livestem_gr' 
     write(10,*) veg_cf%transfer_livestem_gr
     write(10,*) 'veg_cf%transfer_grain_gr' 
     write(10,*) veg_cf%transfer_grain_gr
     write(10,*) 'veg_cf%cpool_froot_storage_gr' 
     write(10,*) veg_cf%cpool_froot_storage_gr
     write(10,*) 'veg_cf%cpool_deadstem_storage_gr' 
     write(10,*) veg_cf%cpool_deadstem_storage_gr
     write(10,*) 'veg_cf%cpool_livestem_gr' 
     write(10,*) veg_cf%cpool_livestem_gr
     write(10,*) 'veg_cf%cpool_livestem_storage_gr' 
     write(10,*) veg_cf%cpool_livestem_storage_gr
     write(10,*) 'veg_cf%cpool_livecroot_gr' 
     write(10,*) veg_cf%cpool_livecroot_gr
     write(10,*) 'veg_cf%cpool_deadcroot_gr' 
     write(10,*) veg_cf%cpool_deadcroot_gr
     write(10,*) 'veg_cf%cpool_froot_gr' 
     write(10,*) veg_cf%cpool_froot_gr
     write(10,*) 'veg_cf%cpool_grain_gr' 
     write(10,*) veg_cf%cpool_grain_gr
     write(10,*) 'veg_cf%transfer_froot_gr' 
     write(10,*) veg_cf%transfer_froot_gr
     write(10,*) 'veg_cf%transfer_leaf_gr' 
     write(10,*) veg_cf%transfer_leaf_gr
     write(10,*) 'veg_cf%cpool_leaf_gr' 
     write(10,*) veg_cf%cpool_leaf_gr
     write(10,*) 'veg_cf%cpool_deadcroot_storage_gr' 
     write(10,*) veg_cf%cpool_deadcroot_storage_gr
     write(10,*) 'veg_cf%transfer_deadcroot_gr' 
     write(10,*) veg_cf%transfer_deadcroot_gr
     write(10,*) 'veg_cf%cpool_deadstem_gr' 
     write(10,*) veg_cf%cpool_deadstem_gr
     write(10,*) 'veg_cf%cpool_livecroot_storage_gr' 
     write(10,*) veg_cf%cpool_livecroot_storage_gr
     write(10,*) 'veg_cf%grainc_xfer_to_grainc' 
     write(10,*) veg_cf%grainc_xfer_to_grainc
     write(10,*) 'veg_cf%cpool_to_leafc_storage' 
     write(10,*) veg_cf%cpool_to_leafc_storage
     write(10,*) 'veg_cf%livestemc_xfer_to_livestemc' 
     write(10,*) veg_cf%livestemc_xfer_to_livestemc
     write(10,*) 'veg_cf%grain_curmr' 
     write(10,*) veg_cf%grain_curmr
     write(10,*) 'veg_cf%deadstemc_xfer_to_deadstemc' 
     write(10,*) veg_cf%deadstemc_xfer_to_deadstemc
     write(10,*) 'veg_cf%cpool_to_deadcrootc_storage' 
     write(10,*) veg_cf%cpool_to_deadcrootc_storage
     write(10,*) 'veg_cf%cpool_grain_storage_gr' 
     write(10,*) veg_cf%cpool_grain_storage_gr
     write(10,*) 'veg_cf%cpool_leaf_storage_gr' 
     write(10,*) veg_cf%cpool_leaf_storage_gr
     write(10,*) 'veg_cf%livestem_curmr' 
     write(10,*) veg_cf%livestem_curmr
     write(10,*) 'veg_cf%livestemc_to_deadstemc' 
     write(10,*) veg_cf%livestemc_to_deadstemc
     write(10,*) 'veg_cf%deadcrootc_storage_to_xfer' 
     write(10,*) veg_cf%deadcrootc_storage_to_xfer
     write(10,*) 'veg_cf%cpool_to_grainc_storage' 
     write(10,*) veg_cf%cpool_to_grainc_storage
     write(10,*) 'veg_cf%froot_curmr' 
     write(10,*) veg_cf%froot_curmr
     write(10,*) 'veg_cf%cpool_froot_storage_gr' 
     write(10,*) veg_cf%cpool_froot_storage_gr
     write(10,*) 'veg_cf%cpool_deadstem_storage_gr' 
     write(10,*) veg_cf%cpool_deadstem_storage_gr
     write(10,*) 'veg_cf%cpool_livestem_gr' 
     write(10,*) veg_cf%cpool_livestem_gr
     write(10,*) 'veg_cf%livecroot_curmr' 
     write(10,*) veg_cf%livecroot_curmr
     write(10,*) 'veg_cf%cpool_livestem_storage_gr' 
     write(10,*) veg_cf%cpool_livestem_storage_gr
     write(10,*) 'veg_cf%livestemc_storage_to_xfer' 
     write(10,*) veg_cf%livestemc_storage_to_xfer
     write(10,*) 'veg_cf%cpool_deadcroot_gr' 
     write(10,*) veg_cf%cpool_deadcroot_gr
     write(10,*) 'veg_cf%cpool_livecroot_gr' 
     write(10,*) veg_cf%cpool_livecroot_gr
     write(10,*) 'veg_cf%cpool_to_grainc' 
     write(10,*) veg_cf%cpool_to_grainc
     write(10,*) 'veg_cf%gresp_storage_to_xfer' 
     write(10,*) veg_cf%gresp_storage_to_xfer
     write(10,*) 'veg_cf%leaf_curmr' 
     write(10,*) veg_cf%leaf_curmr
     write(10,*) 'veg_cf%grainc_storage_to_xfer' 
     write(10,*) veg_cf%grainc_storage_to_xfer
     write(10,*) 'veg_cf%livecrootc_to_deadcrootc' 
     write(10,*) veg_cf%livecrootc_to_deadcrootc
     write(10,*) 'veg_cf%deadcrootc_xfer_to_deadcrootc' 
     write(10,*) veg_cf%deadcrootc_xfer_to_deadcrootc
     write(10,*) 'veg_cf%cpool_to_gresp_storage' 
     write(10,*) veg_cf%cpool_to_gresp_storage
     write(10,*) 'veg_cf%cpool_to_livestemc_storage' 
     write(10,*) veg_cf%cpool_to_livestemc_storage
     write(10,*) 'veg_cf%frootc_xfer_to_frootc' 
     write(10,*) veg_cf%frootc_xfer_to_frootc
     write(10,*) 'veg_cf%frootc_storage_to_xfer' 
     write(10,*) veg_cf%frootc_storage_to_xfer
     write(10,*) 'veg_cf%leafc_xfer_to_leafc' 
     write(10,*) veg_cf%leafc_xfer_to_leafc
     write(10,*) 'veg_cf%cpool_to_deadstemc_storage' 
     write(10,*) veg_cf%cpool_to_deadstemc_storage
     write(10,*) 'veg_cf%cpool_to_frootc_storage' 
     write(10,*) veg_cf%cpool_to_frootc_storage
     write(10,*) 'veg_cf%cpool_froot_gr' 
     write(10,*) veg_cf%cpool_froot_gr
     write(10,*) 'veg_cf%grainc_to_food' 
     write(10,*) veg_cf%grainc_to_food
     write(10,*) 'veg_cf%frootc_to_litter' 
     write(10,*) veg_cf%frootc_to_litter
     write(10,*) 'veg_cf%cpool_to_leafc' 
     write(10,*) veg_cf%cpool_to_leafc
     write(10,*) 'veg_cf%cpool_to_livecrootc' 
     write(10,*) veg_cf%cpool_to_livecrootc
     write(10,*) 'veg_cf%cpool_grain_gr' 
     write(10,*) veg_cf%cpool_grain_gr
     write(10,*) 'veg_cf%livecrootc_xfer_to_livecrootc' 
     write(10,*) veg_cf%livecrootc_xfer_to_livecrootc
     write(10,*) 'veg_cf%cpool_to_frootc' 
     write(10,*) veg_cf%cpool_to_frootc
     write(10,*) 'veg_cf%livecrootc_storage_to_xfer' 
     write(10,*) veg_cf%livecrootc_storage_to_xfer
     write(10,*) 'veg_cf%cpool_leaf_gr' 
     write(10,*) veg_cf%cpool_leaf_gr
     write(10,*) 'veg_cf%cpool_deadcroot_storage_gr' 
     write(10,*) veg_cf%cpool_deadcroot_storage_gr
     write(10,*) 'veg_cf%deadstemc_storage_to_xfer' 
     write(10,*) veg_cf%deadstemc_storage_to_xfer
     write(10,*) 'veg_cf%leafc_storage_to_xfer' 
     write(10,*) veg_cf%leafc_storage_to_xfer
     write(10,*) 'veg_cf%cpool_to_xsmrpool' 
     write(10,*) veg_cf%cpool_to_xsmrpool
     write(10,*) 'veg_cf%leafc_to_litter' 
     write(10,*) veg_cf%leafc_to_litter
     write(10,*) 'veg_cf%cpool_to_deadcrootc' 
     write(10,*) veg_cf%cpool_to_deadcrootc
     write(10,*) 'veg_cf%cpool_to_deadstemc' 
     write(10,*) veg_cf%cpool_to_deadstemc
     write(10,*) 'veg_cf%livestemc_to_litter' 
     write(10,*) veg_cf%livestemc_to_litter
     write(10,*) 'veg_cf%cpool_to_livecrootc_storage' 
     write(10,*) veg_cf%cpool_to_livecrootc_storage
     write(10,*) 'veg_cf%cpool_deadstem_gr' 
     write(10,*) veg_cf%cpool_deadstem_gr
     write(10,*) 'veg_cf%cpool_to_livestemc' 
     write(10,*) veg_cf%cpool_to_livestemc
     write(10,*) 'veg_cf%xsmrpool_to_atm' 
     write(10,*) veg_cf%xsmrpool_to_atm
     write(10,*) 'veg_cf%m_livecrootc_storage_to_litter' 
     write(10,*) veg_cf%m_livecrootc_storage_to_litter
     write(10,*) 'veg_cf%m_gresp_xfer_to_litter' 
     write(10,*) veg_cf%m_gresp_xfer_to_litter
     write(10,*) 'veg_cf%m_deadcrootc_to_litter' 
     write(10,*) veg_cf%m_deadcrootc_to_litter
     write(10,*) 'veg_cf%m_deadcrootc_xfer_to_litter' 
     write(10,*) veg_cf%m_deadcrootc_xfer_to_litter
     write(10,*) 'veg_cf%m_livecrootc_to_litter' 
     write(10,*) veg_cf%m_livecrootc_to_litter
     write(10,*) 'veg_cf%m_deadstemc_storage_to_litter' 
     write(10,*) veg_cf%m_deadstemc_storage_to_litter
     write(10,*) 'veg_cf%m_livestemc_to_litter' 
     write(10,*) veg_cf%m_livestemc_to_litter
     write(10,*) 'veg_cf%m_gresp_storage_to_litter' 
     write(10,*) veg_cf%m_gresp_storage_to_litter
     write(10,*) 'veg_cf%m_leafc_to_litter' 
     write(10,*) veg_cf%m_leafc_to_litter
     write(10,*) 'veg_cf%m_frootc_storage_to_litter' 
     write(10,*) veg_cf%m_frootc_storage_to_litter
     write(10,*) 'veg_cf%m_cpool_to_litter' 
     write(10,*) veg_cf%m_cpool_to_litter
     write(10,*) 'veg_cf%m_livestemc_storage_to_litter' 
     write(10,*) veg_cf%m_livestemc_storage_to_litter
     write(10,*) 'veg_cf%m_frootc_xfer_to_litter' 
     write(10,*) veg_cf%m_frootc_xfer_to_litter
     write(10,*) 'veg_cf%m_deadstemc_to_litter' 
     write(10,*) veg_cf%m_deadstemc_to_litter
     write(10,*) 'veg_cf%m_leafc_xfer_to_litter' 
     write(10,*) veg_cf%m_leafc_xfer_to_litter
     write(10,*) 'veg_cf%m_frootc_to_litter' 
     write(10,*) veg_cf%m_frootc_to_litter
     write(10,*) 'veg_cf%m_deadstemc_xfer_to_litter' 
     write(10,*) veg_cf%m_deadstemc_xfer_to_litter
     write(10,*) 'veg_cf%m_leafc_storage_to_litter' 
     write(10,*) veg_cf%m_leafc_storage_to_litter
     write(10,*) 'veg_cf%m_livestemc_xfer_to_litter' 
     write(10,*) veg_cf%m_livestemc_xfer_to_litter
     write(10,*) 'veg_cf%m_deadcrootc_storage_to_litter' 
     write(10,*) veg_cf%m_deadcrootc_storage_to_litter
     write(10,*) 'veg_cf%m_livecrootc_xfer_to_litter' 
     write(10,*) veg_cf%m_livecrootc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_deadstemc_storage_to_litter' 
     write(10,*) veg_cf%hrv_deadstemc_storage_to_litter
     write(10,*) 'veg_cf%hrv_cpool_to_litter' 
     write(10,*) veg_cf%hrv_cpool_to_litter
     write(10,*) 'veg_cf%hrv_livestemc_storage_to_litter' 
     write(10,*) veg_cf%hrv_livestemc_storage_to_litter
     write(10,*) 'veg_cf%hrv_leafc_storage_to_litter' 
     write(10,*) veg_cf%hrv_leafc_storage_to_litter
     write(10,*) 'veg_cf%hrv_deadcrootc_to_litter' 
     write(10,*) veg_cf%hrv_deadcrootc_to_litter
     write(10,*) 'veg_cf%hrv_deadstemc_to_prod10c' 
     write(10,*) veg_cf%hrv_deadstemc_to_prod10c
     write(10,*) 'veg_cf%hrv_deadstemc_to_prod100c' 
     write(10,*) veg_cf%hrv_deadstemc_to_prod100c
     write(10,*) 'veg_cf%hrv_leafc_xfer_to_litter' 
     write(10,*) veg_cf%hrv_leafc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_leafc_to_litter' 
     write(10,*) veg_cf%hrv_leafc_to_litter
     write(10,*) 'veg_cf%hrv_livecrootc_xfer_to_litter' 
     write(10,*) veg_cf%hrv_livecrootc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_frootc_to_litter' 
     write(10,*) veg_cf%hrv_frootc_to_litter
     write(10,*) 'veg_cf%hrv_deadcrootc_xfer_to_litter' 
     write(10,*) veg_cf%hrv_deadcrootc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_xsmrpool_to_atm' 
     write(10,*) veg_cf%hrv_xsmrpool_to_atm
     write(10,*) 'veg_cf%hrv_livestemc_to_litter' 
     write(10,*) veg_cf%hrv_livestemc_to_litter
     write(10,*) 'veg_cf%hrv_deadcrootc_storage_to_litter' 
     write(10,*) veg_cf%hrv_deadcrootc_storage_to_litter
     write(10,*) 'veg_cf%hrv_deadstemc_xfer_to_litter' 
     write(10,*) veg_cf%hrv_deadstemc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_gresp_storage_to_litter' 
     write(10,*) veg_cf%hrv_gresp_storage_to_litter
     write(10,*) 'veg_cf%hrv_livecrootc_storage_to_litter' 
     write(10,*) veg_cf%hrv_livecrootc_storage_to_litter
     write(10,*) 'veg_cf%hrv_livestemc_xfer_to_litter' 
     write(10,*) veg_cf%hrv_livestemc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_frootc_xfer_to_litter' 
     write(10,*) veg_cf%hrv_frootc_xfer_to_litter
     write(10,*) 'veg_cf%hrv_livecrootc_to_litter' 
     write(10,*) veg_cf%hrv_livecrootc_to_litter
     write(10,*) 'veg_cf%hrv_frootc_storage_to_litter' 
     write(10,*) veg_cf%hrv_frootc_storage_to_litter
     write(10,*) 'veg_cf%hrv_gresp_xfer_to_litter' 
     write(10,*) veg_cf%hrv_gresp_xfer_to_litter
     write(10,*) 'veg_cf%m_livecrootc_xfer_to_fire' 
     write(10,*) veg_cf%m_livecrootc_xfer_to_fire
     write(10,*) 'veg_cf%m_leafc_xfer_to_fire' 
     write(10,*) veg_cf%m_leafc_xfer_to_fire
     write(10,*) 'veg_cf%m_livestemc_to_fire' 
     write(10,*) veg_cf%m_livestemc_to_fire
     write(10,*) 'veg_cf%m_livestemc_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_livestemc_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_frootc_storage_to_fire' 
     write(10,*) veg_cf%m_frootc_storage_to_fire
     write(10,*) 'veg_cf%m_frootc_to_litter_fire' 
     write(10,*) veg_cf%m_frootc_to_litter_fire
     write(10,*) 'veg_cf%m_deadcrootc_to_litter_fire' 
     write(10,*) veg_cf%m_deadcrootc_to_litter_fire
     write(10,*) 'veg_cf%m_deadstemc_storage_to_litter_fire' 
     write(10,*) veg_cf%m_deadstemc_storage_to_litter_fire
     write(10,*) 'veg_cf%m_livestemc_storage_to_litter_fire' 
     write(10,*) veg_cf%m_livestemc_storage_to_litter_fire
     write(10,*) 'veg_cf%m_livestemc_to_litter_fire' 
     write(10,*) veg_cf%m_livestemc_to_litter_fire
     write(10,*) 'veg_cf%m_deadstemc_xfer_to_fire' 
     write(10,*) veg_cf%m_deadstemc_xfer_to_fire
     write(10,*) 'veg_cf%m_gresp_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_gresp_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_deadstemc_to_litter_fire' 
     write(10,*) veg_cf%m_deadstemc_to_litter_fire
     write(10,*) 'veg_cf%m_leafc_to_fire' 
     write(10,*) veg_cf%m_leafc_to_fire
     write(10,*) 'veg_cf%m_deadstemc_storage_to_fire' 
     write(10,*) veg_cf%m_deadstemc_storage_to_fire
     write(10,*) 'veg_cf%m_deadcrootc_storage_to_fire' 
     write(10,*) veg_cf%m_deadcrootc_storage_to_fire
     write(10,*) 'veg_cf%m_deadcrootc_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_deadcrootc_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_livestemc_to_deadstemc_fire' 
     write(10,*) veg_cf%m_livestemc_to_deadstemc_fire
     write(10,*) 'veg_cf%m_frootc_storage_to_litter_fire' 
     write(10,*) veg_cf%m_frootc_storage_to_litter_fire
     write(10,*) 'veg_cf%m_frootc_to_fire' 
     write(10,*) veg_cf%m_frootc_to_fire
     write(10,*) 'veg_cf%m_deadcrootc_xfer_to_fire' 
     write(10,*) veg_cf%m_deadcrootc_xfer_to_fire
     write(10,*) 'veg_cf%m_gresp_storage_to_litter_fire' 
     write(10,*) veg_cf%m_gresp_storage_to_litter_fire
     write(10,*) 'veg_cf%m_livestemc_xfer_to_fire' 
     write(10,*) veg_cf%m_livestemc_xfer_to_fire
     write(10,*) 'veg_cf%m_frootc_xfer_to_fire' 
     write(10,*) veg_cf%m_frootc_xfer_to_fire
     write(10,*) 'veg_cf%m_cpool_to_litter_fire' 
     write(10,*) veg_cf%m_cpool_to_litter_fire
     write(10,*) 'veg_cf%m_livecrootc_storage_to_litter_fire' 
     write(10,*) veg_cf%m_livecrootc_storage_to_litter_fire
     write(10,*) 'veg_cf%m_deadcrootc_storage_to_litter_fire' 
     write(10,*) veg_cf%m_deadcrootc_storage_to_litter_fire
     write(10,*) 'veg_cf%m_deadstemc_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_deadstemc_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_deadcrootc_to_fire' 
     write(10,*) veg_cf%m_deadcrootc_to_fire
     write(10,*) 'veg_cf%m_livecrootc_storage_to_fire' 
     write(10,*) veg_cf%m_livecrootc_storage_to_fire
     write(10,*) 'veg_cf%m_gresp_xfer_to_fire' 
     write(10,*) veg_cf%m_gresp_xfer_to_fire
     write(10,*) 'veg_cf%m_leafc_to_litter_fire' 
     write(10,*) veg_cf%m_leafc_to_litter_fire
     write(10,*) 'veg_cf%m_frootc_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_frootc_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_livestemc_storage_to_fire' 
     write(10,*) veg_cf%m_livestemc_storage_to_fire
     write(10,*) 'veg_cf%m_livecrootc_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_livecrootc_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_livecrootc_to_litter_fire' 
     write(10,*) veg_cf%m_livecrootc_to_litter_fire
     write(10,*) 'veg_cf%m_leafc_xfer_to_litter_fire' 
     write(10,*) veg_cf%m_leafc_xfer_to_litter_fire
     write(10,*) 'veg_cf%m_deadstemc_to_fire' 
     write(10,*) veg_cf%m_deadstemc_to_fire
     write(10,*) 'veg_cf%m_livecrootc_to_fire' 
     write(10,*) veg_cf%m_livecrootc_to_fire
     write(10,*) 'veg_cf%m_gresp_storage_to_fire' 
     write(10,*) veg_cf%m_gresp_storage_to_fire
     write(10,*) 'veg_cf%m_leafc_storage_to_litter_fire' 
     write(10,*) veg_cf%m_leafc_storage_to_litter_fire
     write(10,*) 'veg_cf%m_leafc_storage_to_fire' 
     write(10,*) veg_cf%m_leafc_storage_to_fire
     write(10,*) 'veg_cf%m_cpool_to_fire' 
     write(10,*) veg_cf%m_cpool_to_fire
     write(10,*) 'veg_cf%m_livecrootc_to_deadcrootc_fire' 
     write(10,*) veg_cf%m_livecrootc_to_deadcrootc_fire
     write(10,*) 'canopystate_vars%laisun_patch' 
     write(10,*) canopystate_vars%laisun_patch
     write(10,*) 'canopystate_vars%laisha_patch' 
     write(10,*) canopystate_vars%laisha_patch
     if ( use_c13 ) then
     write(10,*) 'c13_veg_cf%psnshade_to_cpool' 
     write(10,*) c13_veg_cf%psnshade_to_cpool
     write(10,*) 'c13_veg_cf%psnsun_to_cpool' 
     write(10,*) c13_veg_cf%psnsun_to_cpool
     endif
     if ( use_c14 ) then
     write(10,*) 'c14_veg_cf%psnshade_to_cpool' 
     write(10,*) c14_veg_cf%psnshade_to_cpool
     write(10,*) 'c14_veg_cf%psnsun_to_cpool' 
     write(10,*) c14_veg_cf%psnsun_to_cpool
     endif
     write(10,*) 'veg_es%gdd820' 
     write(10,*) veg_es%gdd820
     write(10,*) 'veg_es%gdd1020' 
     write(10,*) veg_es%gdd1020
     write(10,*) 'veg_es%gdd020' 
     write(10,*) veg_es%gdd020
     write(10,*) 'veg_es%t_ref2m' 
     write(10,*) veg_es%t_ref2m
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%cpool' 
     write(10,*) veg_cs%cpool
     write(10,*) 'veg_cs%cpool' 
     write(10,*) veg_cs%cpool
     write(10,*) 'veg_cs%gresp_xfer' 
     write(10,*) veg_cs%gresp_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%grainc' 
     write(10,*) veg_cs%grainc
     write(10,*) 'veg_cs%leafc_storage' 
     write(10,*) veg_cs%leafc_storage
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     write(10,*) 'veg_cs%deadcrootc' 
     write(10,*) veg_cs%deadcrootc
     write(10,*) 'veg_cs%livecrootc_storage' 
     write(10,*) veg_cs%livecrootc_storage
     write(10,*) 'veg_cs%livestemc_storage' 
     write(10,*) veg_cs%livestemc_storage
     write(10,*) 'veg_cs%frootc_storage' 
     write(10,*) veg_cs%frootc_storage
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%cropseedc_deficit' 
     write(10,*) veg_cs%cropseedc_deficit
     write(10,*) 'veg_cs%xsmrpool' 
     write(10,*) veg_cs%xsmrpool
     write(10,*) 'veg_cs%grainc_storage' 
     write(10,*) veg_cs%grainc_storage
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%gresp_storage' 
     write(10,*) veg_cs%gresp_storage
     write(10,*) 'veg_cs%grainc_xfer' 
     write(10,*) veg_cs%grainc_xfer
     write(10,*) 'veg_cs%frootc' 
     write(10,*) veg_cs%frootc
     write(10,*) 'veg_cs%livestemc' 
     write(10,*) veg_cs%livestemc
     write(10,*) 'veg_cs%leafc' 
     write(10,*) veg_cs%leafc
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%deadstemc_storage' 
     write(10,*) veg_cs%deadstemc_storage
     write(10,*) 'veg_cs%livecrootc' 
     write(10,*) veg_cs%livecrootc
     write(10,*) 'veg_cs%deadcrootc_storage' 
     write(10,*) veg_cs%deadcrootc_storage
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     write(10,*) 'veg_cs%cpool' 
     write(10,*) veg_cs%cpool
     write(10,*) 'veg_cs%gresp_xfer' 
     write(10,*) veg_cs%gresp_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%leafc_storage' 
     write(10,*) veg_cs%leafc_storage
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     write(10,*) 'veg_cs%deadcrootc' 
     write(10,*) veg_cs%deadcrootc
     write(10,*) 'veg_cs%livecrootc_storage' 
     write(10,*) veg_cs%livecrootc_storage
     write(10,*) 'veg_cs%livestemc_storage' 
     write(10,*) veg_cs%livestemc_storage
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%frootc_storage' 
     write(10,*) veg_cs%frootc_storage
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%gresp_storage' 
     write(10,*) veg_cs%gresp_storage
     write(10,*) 'veg_cs%frootc' 
     write(10,*) veg_cs%frootc
     write(10,*) 'veg_cs%livestemc' 
     write(10,*) veg_cs%livestemc
     write(10,*) 'veg_cs%leafc' 
     write(10,*) veg_cs%leafc
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%deadstemc_storage' 
     write(10,*) veg_cs%deadstemc_storage
     write(10,*) 'veg_cs%livecrootc' 
     write(10,*) veg_cs%livecrootc
     write(10,*) 'veg_cs%deadcrootc_storage' 
     write(10,*) veg_cs%deadcrootc_storage
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     write(10,*) 'veg_cs%cpool' 
     write(10,*) veg_cs%cpool
     write(10,*) 'veg_cs%gresp_xfer' 
     write(10,*) veg_cs%gresp_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%grainc' 
     write(10,*) veg_cs%grainc
     write(10,*) 'veg_cs%leafc_storage' 
     write(10,*) veg_cs%leafc_storage
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     write(10,*) 'veg_cs%deadcrootc' 
     write(10,*) veg_cs%deadcrootc
     write(10,*) 'veg_cs%livecrootc_storage' 
     write(10,*) veg_cs%livecrootc_storage
     write(10,*) 'veg_cs%livestemc_storage' 
     write(10,*) veg_cs%livestemc_storage
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%frootc_storage' 
     write(10,*) veg_cs%frootc_storage
     write(10,*) 'veg_cs%xsmrpool' 
     write(10,*) veg_cs%xsmrpool
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%gresp_storage' 
     write(10,*) veg_cs%gresp_storage
     write(10,*) 'veg_cs%frootc' 
     write(10,*) veg_cs%frootc
     write(10,*) 'veg_cs%livestemc' 
     write(10,*) veg_cs%livestemc
     write(10,*) 'veg_cs%leafc' 
     write(10,*) veg_cs%leafc
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%deadstemc_storage' 
     write(10,*) veg_cs%deadstemc_storage
     write(10,*) 'veg_cs%livecrootc' 
     write(10,*) veg_cs%livecrootc
     write(10,*) 'veg_cs%deadcrootc_storage' 
     write(10,*) veg_cs%deadcrootc_storage
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     write(10,*) 'veg_cs%totvegc' 
     write(10,*) veg_cs%totvegc
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     write(10,*) 'veg_cs%leafc' 
     write(10,*) veg_cs%leafc
     write(10,*) 'veg_cs%cpool' 
     write(10,*) veg_cs%cpool
     write(10,*) 'veg_cs%gresp_xfer' 
     write(10,*) veg_cs%gresp_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%leafc_storage' 
     write(10,*) veg_cs%leafc_storage
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     write(10,*) 'veg_cs%deadcrootc' 
     write(10,*) veg_cs%deadcrootc
     write(10,*) 'veg_cs%livecrootc_storage' 
     write(10,*) veg_cs%livecrootc_storage
     write(10,*) 'veg_cs%livestemc_storage' 
     write(10,*) veg_cs%livestemc_storage
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%frootc_storage' 
     write(10,*) veg_cs%frootc_storage
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%gresp_storage' 
     write(10,*) veg_cs%gresp_storage
     write(10,*) 'veg_cs%frootc' 
     write(10,*) veg_cs%frootc
     write(10,*) 'veg_cs%livestemc' 
     write(10,*) veg_cs%livestemc
     write(10,*) 'veg_cs%leafc' 
     write(10,*) veg_cs%leafc
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%deadstemc_storage' 
     write(10,*) veg_cs%deadstemc_storage
     write(10,*) 'veg_cs%livecrootc' 
     write(10,*) veg_cs%livecrootc
     write(10,*) 'veg_cs%deadcrootc_storage' 
     write(10,*) veg_cs%deadcrootc_storage
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%deadcrootp_storage' 
     write(10,*) veg_ps%deadcrootp_storage
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%livecrootp' 
     write(10,*) veg_ps%livecrootp
     write(10,*) 'veg_ps%cropseedp_deficit' 
     write(10,*) veg_ps%cropseedp_deficit
     write(10,*) 'veg_ps%retransp' 
     write(10,*) veg_ps%retransp
     write(10,*) 'veg_ps%leafp_storage' 
     write(10,*) veg_ps%leafp_storage
     write(10,*) 'veg_ps%deadstemp_storage' 
     write(10,*) veg_ps%deadstemp_storage
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%ppool' 
     write(10,*) veg_ps%ppool
     write(10,*) 'veg_ps%livecrootp_storage' 
     write(10,*) veg_ps%livecrootp_storage
     write(10,*) 'veg_ps%leafp' 
     write(10,*) veg_ps%leafp
     write(10,*) 'veg_ps%grainp' 
     write(10,*) veg_ps%grainp
     write(10,*) 'veg_ps%livestemp_storage' 
     write(10,*) veg_ps%livestemp_storage
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%grainp_storage' 
     write(10,*) veg_ps%grainp_storage
     write(10,*) 'veg_ps%livestemp' 
     write(10,*) veg_ps%livestemp
     write(10,*) 'veg_ps%grainp_xfer' 
     write(10,*) veg_ps%grainp_xfer
     write(10,*) 'veg_ps%deadcrootp' 
     write(10,*) veg_ps%deadcrootp
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%frootp_storage' 
     write(10,*) veg_ps%frootp_storage
     write(10,*) 'veg_ps%frootp' 
     write(10,*) veg_ps%frootp
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'veg_ps%leafp' 
     write(10,*) veg_ps%leafp
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%deadcrootp_storage' 
     write(10,*) veg_ps%deadcrootp_storage
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%livecrootp' 
     write(10,*) veg_ps%livecrootp
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%retransp' 
     write(10,*) veg_ps%retransp
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%deadstemp_storage' 
     write(10,*) veg_ps%deadstemp_storage
     write(10,*) 'veg_ps%livestemp' 
     write(10,*) veg_ps%livestemp
     write(10,*) 'veg_ps%leafp_storage' 
     write(10,*) veg_ps%leafp_storage
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'veg_ps%deadcrootp' 
     write(10,*) veg_ps%deadcrootp
     write(10,*) 'veg_ps%frootp_storage' 
     write(10,*) veg_ps%frootp_storage
     write(10,*) 'veg_ps%frootp' 
     write(10,*) veg_ps%frootp
     write(10,*) 'veg_ps%livecrootp_storage' 
     write(10,*) veg_ps%livecrootp_storage
     write(10,*) 'veg_ps%livestemp_storage' 
     write(10,*) veg_ps%livestemp_storage
     write(10,*) 'veg_ps%ppool' 
     write(10,*) veg_ps%ppool
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%deadcrootp_storage' 
     write(10,*) veg_ps%deadcrootp_storage
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%livecrootp' 
     write(10,*) veg_ps%livecrootp
     write(10,*) 'veg_ps%retransp' 
     write(10,*) veg_ps%retransp
     write(10,*) 'veg_ps%leafp_storage' 
     write(10,*) veg_ps%leafp_storage
     write(10,*) 'veg_ps%deadstemp_storage' 
     write(10,*) veg_ps%deadstemp_storage
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%livecrootp_storage' 
     write(10,*) veg_ps%livecrootp_storage
     write(10,*) 'veg_ps%ppool' 
     write(10,*) veg_ps%ppool
     write(10,*) 'veg_ps%leafp' 
     write(10,*) veg_ps%leafp
     write(10,*) 'veg_ps%grainp' 
     write(10,*) veg_ps%grainp
     write(10,*) 'veg_ps%livestemp_storage' 
     write(10,*) veg_ps%livestemp_storage
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%livestemp' 
     write(10,*) veg_ps%livestemp
     write(10,*) 'veg_ps%deadcrootp' 
     write(10,*) veg_ps%deadcrootp
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%frootp_storage' 
     write(10,*) veg_ps%frootp_storage
     write(10,*) 'veg_ps%frootp' 
     write(10,*) veg_ps%frootp
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'top_as%windbot' 
     write(10,*) top_as%windbot
     write(10,*) 'top_as%pbot' 
     write(10,*) top_as%pbot
     write(10,*) 'top_as%rhbot' 
     write(10,*) top_as%rhbot
     write(10,*) 'veg_ef%eflx_soil_grnd' 
     write(10,*) veg_ef%eflx_soil_grnd
     write(10,*) 'veg_ef%netrad' 
     write(10,*) veg_ef%netrad
     write(10,*) 'crop_vars%prev_xt_bar_patch' 
     write(10,*) crop_vars%prev_xt_bar_patch
     write(10,*) 'crop_vars%xt_patch' 
     write(10,*) crop_vars%xt_patch
     write(10,*) 'crop_vars%p2ETo_patch' 
     write(10,*) crop_vars%p2ETo_patch
     write(10,*) 'crop_vars%cvp_patch' 
     write(10,*) crop_vars%cvp_patch
     write(10,*) 'crop_vars%ETo_patch' 
     write(10,*) crop_vars%ETo_patch
     write(10,*) 'crop_vars%xp_patch' 
     write(10,*) crop_vars%xp_patch
     write(10,*) 'crop_vars%cvt_patch' 
     write(10,*) crop_vars%cvt_patch
     write(10,*) 'crop_vars%plantmonth_patch' 
     write(10,*) crop_vars%plantmonth_patch
     write(10,*) 'crop_vars%xt_bar_patch' 
     write(10,*) crop_vars%xt_bar_patch
     write(10,*) 'crop_vars%p2ETo_bar_patch' 
     write(10,*) crop_vars%p2ETo_bar_patch
     write(10,*) 'crop_vars%prev_p2ETo_bar_patch' 
     write(10,*) crop_vars%prev_p2ETo_bar_patch
     write(10,*) 'crop_vars%prev_xp_bar_patch' 
     write(10,*) crop_vars%prev_xp_bar_patch
     write(10,*) 'crop_vars%P2E_rm_patch' 
     write(10,*) crop_vars%P2E_rm_patch
     write(10,*) 'crop_vars%xp_bar_patch' 
     write(10,*) crop_vars%xp_bar_patch
     write(10,*) 'crop_vars%plantday_patch' 
     write(10,*) crop_vars%plantday_patch
     write(10,*) 'crop_vars%vf_patch' 
     write(10,*) crop_vars%vf_patch
     write(10,*) 'crop_vars%crpyld_patch' 
     write(10,*) crop_vars%crpyld_patch
     write(10,*) 'crop_vars%cropplant_patch' 
     write(10,*) crop_vars%cropplant_patch
     write(10,*) 'crop_vars%croplive_patch' 
     write(10,*) crop_vars%croplive_patch
     write(10,*) 'crop_vars%harvdate_patch' 
     write(10,*) crop_vars%harvdate_patch
     write(10,*) 'crop_vars%dmyield_patch' 
     write(10,*) crop_vars%dmyield_patch
     write(10,*) 'crop_vars%harvday_patch' 
     write(10,*) crop_vars%harvday_patch
     write(10,*) 'crop_vars%vf_patch' 
     write(10,*) crop_vars%vf_patch
     write(10,*) 'crop_vars%crpyld_patch' 
     write(10,*) crop_vars%crpyld_patch
     write(10,*) 'crop_vars%dmyield_patch' 
     write(10,*) crop_vars%dmyield_patch
     write(10,*) 'soilstate_vars%rootfr_patch' 
     write(10,*) soilstate_vars%rootfr_patch
     write(10,*) 'soilstate_vars%root_depth_patch' 
     write(10,*) soilstate_vars%root_depth_patch
     write(10,*) 'col_pp%nlevbed' 
     write(10,*) col_pp%nlevbed
     write(10,*) 'decomp_cascade_con%cascade_donor_pool' 
     write(10,*) decomp_cascade_con%cascade_donor_pool
     write(10,*) 'decomp_cascade_con%cascade_receiver_pool' 
     write(10,*) decomp_cascade_con%cascade_receiver_pool
     write(10,*) 'decomp_cascade_con%cascade_donor_pool' 
     write(10,*) decomp_cascade_con%cascade_donor_pool
     write(10,*) 'decomp_cascade_con%cascade_receiver_pool' 
     write(10,*) decomp_cascade_con%cascade_receiver_pool
     write(10,*) 'decomp_cascade_con%cascade_donor_pool' 
     write(10,*) decomp_cascade_con%cascade_donor_pool
     write(10,*) 'decomp_cascade_con%cascade_receiver_pool' 
     write(10,*) decomp_cascade_con%cascade_receiver_pool
     write(10,*) 'col_cs%decomp_som2c_vr' 
     write(10,*) col_cs%decomp_som2c_vr
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     write(10,*) 'col_cs%prod100c' 
     write(10,*) col_cs%prod100c
     write(10,*) 'col_cs%prod10c' 
     write(10,*) col_cs%prod10c
     write(10,*) 'col_cs%prod1c' 
     write(10,*) col_cs%prod1c
     write(10,*) 'col_cs%leafc' 
     write(10,*) col_cs%leafc
     write(10,*) 'col_cs%fuelc' 
     write(10,*) col_cs%fuelc
     write(10,*) 'col_cs%deadstemc' 
     write(10,*) col_cs%deadstemc
     write(10,*) 'col_cs%totvegc' 
     write(10,*) col_cs%totvegc
     write(10,*) 'col_cs%rootc' 
     write(10,*) col_cs%rootc
     write(10,*) 'col_cs%fuelc_crop' 
     write(10,*) col_cs%fuelc_crop
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     write(10,*) 'col_ps%solutionp_vr' 
     write(10,*) col_ps%solutionp_vr
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_ps%prod100p' 
     write(10,*) col_ps%prod100p
     write(10,*) 'col_ps%prod10p' 
     write(10,*) col_ps%prod10p
     write(10,*) 'col_ps%prod1p' 
     write(10,*) col_ps%prod1p
     if ( use_c13 ) then
     write(10,*) 'c13_col_cs%decomp_cpools_vr' 
     write(10,*) c13_col_cs%decomp_cpools_vr
     write(10,*) 'c13_col_cs%prod100c' 
     write(10,*) c13_col_cs%prod100c
     write(10,*) 'c13_col_cs%prod10c' 
     write(10,*) c13_col_cs%prod10c
     write(10,*) 'c13_col_cs%prod1c' 
     write(10,*) c13_col_cs%prod1c
     endif
     if ( use_c14 ) then
     write(10,*) 'c14_col_cs%decomp_cpools_vr' 
     write(10,*) c14_col_cs%decomp_cpools_vr
     write(10,*) 'c14_col_cs%prod100c' 
     write(10,*) c14_col_cs%prod100c
     write(10,*) 'c14_col_cs%prod10c' 
     write(10,*) c14_col_cs%prod10c
     write(10,*) 'c14_col_cs%prod1c' 
     write(10,*) c14_col_cs%prod1c
     write(10,*) 'c14_col_cs%seedc' 
     write(10,*) c14_col_cs%seedc
     write(10,*) 'c14_col_cs%decomp_cpools_vr' 
     write(10,*) c14_col_cs%decomp_cpools_vr
     endif
     if ( use_c13 ) then
     write(10,*) 'c13_col_cf%decomp_cpools_transport_tendency' 
     write(10,*) c13_col_cf%decomp_cpools_transport_tendency
     write(10,*) 'c13_col_cf%prod10c_loss' 
     write(10,*) c13_col_cf%prod10c_loss
     write(10,*) 'c13_col_cf%prod100c_loss' 
     write(10,*) c13_col_cf%prod100c_loss
     write(10,*) 'c13_col_cf%prod1c_loss' 
     write(10,*) c13_col_cf%prod1c_loss
     endif
     if ( use_c14 ) then
     write(10,*) 'c14_col_cf%decomp_cpools_transport_tendency' 
     write(10,*) c14_col_cf%decomp_cpools_transport_tendency
     write(10,*) 'c14_col_cf%prod10c_loss' 
     write(10,*) c14_col_cf%prod10c_loss
     write(10,*) 'c14_col_cf%prod100c_loss' 
     write(10,*) c14_col_cf%prod100c_loss
     write(10,*) 'c14_col_cf%prod1c_loss' 
     write(10,*) c14_col_cf%prod1c_loss
     write(10,*) 'c14_veg_cs%cpool' 
     write(10,*) c14_veg_cs%cpool
     write(10,*) 'c14_veg_cs%gresp_xfer' 
     write(10,*) c14_veg_cs%gresp_xfer
     write(10,*) 'c14_veg_cs%deadstemc_xfer' 
     write(10,*) c14_veg_cs%deadstemc_xfer
     write(10,*) 'c14_veg_cs%leafc_storage' 
     write(10,*) c14_veg_cs%leafc_storage
     write(10,*) 'c14_veg_cs%deadstemc' 
     write(10,*) c14_veg_cs%deadstemc
     write(10,*) 'c14_veg_cs%deadcrootc' 
     write(10,*) c14_veg_cs%deadcrootc
     write(10,*) 'c14_veg_cs%livecrootc_storage' 
     write(10,*) c14_veg_cs%livecrootc_storage
     write(10,*) 'c14_veg_cs%livestemc_storage' 
     write(10,*) c14_veg_cs%livestemc_storage
     write(10,*) 'c14_veg_cs%frootc_storage' 
     write(10,*) c14_veg_cs%frootc_storage
     write(10,*) 'c14_veg_cs%livestemc_xfer' 
     write(10,*) c14_veg_cs%livestemc_xfer
     write(10,*) 'c14_veg_cs%xsmrpool' 
     write(10,*) c14_veg_cs%xsmrpool
     write(10,*) 'c14_veg_cs%leafc_xfer' 
     write(10,*) c14_veg_cs%leafc_xfer
     write(10,*) 'c14_veg_cs%gresp_storage' 
     write(10,*) c14_veg_cs%gresp_storage
     write(10,*) 'c14_veg_cs%ctrunc' 
     write(10,*) c14_veg_cs%ctrunc
     write(10,*) 'c14_veg_cs%frootc' 
     write(10,*) c14_veg_cs%frootc
     write(10,*) 'c14_veg_cs%livestemc' 
     write(10,*) c14_veg_cs%livestemc
     write(10,*) 'c14_veg_cs%leafc' 
     write(10,*) c14_veg_cs%leafc
     write(10,*) 'c14_veg_cs%livecrootc_xfer' 
     write(10,*) c14_veg_cs%livecrootc_xfer
     write(10,*) 'c14_veg_cs%deadstemc_storage' 
     write(10,*) c14_veg_cs%deadstemc_storage
     write(10,*) 'c14_veg_cs%livecrootc' 
     write(10,*) c14_veg_cs%livecrootc
     write(10,*) 'c14_veg_cs%deadcrootc_storage' 
     write(10,*) c14_veg_cs%deadcrootc_storage
     write(10,*) 'c14_veg_cs%frootc_xfer' 
     write(10,*) c14_veg_cs%frootc_xfer
     write(10,*) 'c14_veg_cs%deadcrootc_xfer' 
     write(10,*) c14_veg_cs%deadcrootc_xfer
     endif
     close(10)
end subroutine 
subroutine update_vars_SurfaceAlbedo(gpu)
     use clm_instMod, only : surfalb_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SurfaceAlbedo.txt"
     else
          file='cpu_SurfaceAlbedo.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc surfalb_vars%albgrd_dst_col, & 
     !$acc surfalb_vars%fabd_sun_patch, & 
     !$acc surfalb_vars%fabd_sun_z_patch, & 
     !$acc surfalb_vars%ftdd_patch, & 
     !$acc surfalb_vars%albd_patch, & 
     !$acc surfalb_vars%ftid_patch, & 
     !$acc surfalb_vars%fabi_sha_z_patch, & 
     !$acc surfalb_vars%coszen_col, & 
     !$acc surfalb_vars%fabi_patch, & 
     !$acc surfalb_vars%fabi_sun_z_patch, & 
     !$acc surfalb_vars%tsai_z_patch, & 
     !$acc surfalb_vars%albgrd_pur_col, & 
     !$acc surfalb_vars%flx_absin_col, & 
     !$acc surfalb_vars%nrad_patch, & 
     !$acc surfalb_vars%fabi_sun_patch, & 
     !$acc surfalb_vars%ncan_patch, & 
     !$acc surfalb_vars%vcmaxcintsun_patch, & 
     !$acc surfalb_vars%albgrd_col, & 
     !$acc surfalb_vars%albgri_bc_col, & 
     !$acc surfalb_vars%ftii_patch, & 
     !$acc surfalb_vars%tlai_z_patch, & 
     !$acc surfalb_vars%fabd_sha_patch, & 
     !$acc surfalb_vars%flx_absdn_col, & 
     !$acc surfalb_vars%fabd_sha_z_patch, & 
     !$acc surfalb_vars%albgri_oc_col, & 
     !$acc surfalb_vars%albgrd_oc_col, & 
     !$acc surfalb_vars%albsni_hst_col, & 
     !$acc surfalb_vars%albgri_pur_col, & 
     !$acc surfalb_vars%albgri_dst_col, & 
     !$acc surfalb_vars%fsun_z_patch, & 
     !$acc surfalb_vars%vcmaxcintsha_patch, & 
     !$acc surfalb_vars%albgri_col, & 
     !$acc surfalb_vars%albi_patch, & 
     !$acc surfalb_vars%flx_absdv_col, & 
     !$acc surfalb_vars%fabd_patch, & 
     !$acc surfalb_vars%albsod_col, & 
     !$acc surfalb_vars%fabi_sha_patch, & 
     !$acc surfalb_vars%albgrd_bc_col, & 
     !$acc surfalb_vars%flx_absiv_col, & 
     !$acc surfalb_vars%albsnd_hst_col, & 
     !$acc surfalb_vars%albsoi_col, & 
     !$acc surfalb_vars%albsod_col, & 
     !$acc surfalb_vars%albsoi_col, & 
     !$acc surfalb_vars%fabd_sun_patch, & 
     !$acc surfalb_vars%fsun_z_patch, & 
     !$acc surfalb_vars%fabi_sun_z_patch, & 
     !$acc surfalb_vars%vcmaxcintsha_patch, & 
     !$acc surfalb_vars%albi_patch, & 
     !$acc surfalb_vars%fabd_sun_z_patch, & 
     !$acc surfalb_vars%ftii_patch, & 
     !$acc surfalb_vars%fabd_patch, & 
     !$acc surfalb_vars%fabi_sha_patch, & 
     !$acc surfalb_vars%ftdd_patch, & 
     !$acc surfalb_vars%fabd_sha_z_patch, & 
     !$acc surfalb_vars%albd_patch, & 
     !$acc surfalb_vars%fabd_sha_patch, & 
     !$acc surfalb_vars%fabi_sun_patch, & 
     !$acc surfalb_vars%fabi_sha_z_patch, & 
     !$acc surfalb_vars%ftid_patch, & 
     !$acc surfalb_vars%vcmaxcintsun_patch, & 
     !$acc surfalb_vars%fabi_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'surfalb_vars%albgrd_dst_col' 
     write(10,*) surfalb_vars%albgrd_dst_col
     write(10,*) 'surfalb_vars%fabd_sun_patch' 
     write(10,*) surfalb_vars%fabd_sun_patch
     write(10,*) 'surfalb_vars%fabd_sun_z_patch' 
     write(10,*) surfalb_vars%fabd_sun_z_patch
     write(10,*) 'surfalb_vars%ftdd_patch' 
     write(10,*) surfalb_vars%ftdd_patch
     write(10,*) 'surfalb_vars%albd_patch' 
     write(10,*) surfalb_vars%albd_patch
     write(10,*) 'surfalb_vars%ftid_patch' 
     write(10,*) surfalb_vars%ftid_patch
     write(10,*) 'surfalb_vars%fabi_sha_z_patch' 
     write(10,*) surfalb_vars%fabi_sha_z_patch
     write(10,*) 'surfalb_vars%coszen_col' 
     write(10,*) surfalb_vars%coszen_col
     write(10,*) 'surfalb_vars%fabi_patch' 
     write(10,*) surfalb_vars%fabi_patch
     write(10,*) 'surfalb_vars%fabi_sun_z_patch' 
     write(10,*) surfalb_vars%fabi_sun_z_patch
     write(10,*) 'surfalb_vars%tsai_z_patch' 
     write(10,*) surfalb_vars%tsai_z_patch
     write(10,*) 'surfalb_vars%albgrd_pur_col' 
     write(10,*) surfalb_vars%albgrd_pur_col
     write(10,*) 'surfalb_vars%flx_absin_col' 
     write(10,*) surfalb_vars%flx_absin_col
     write(10,*) 'surfalb_vars%nrad_patch' 
     write(10,*) surfalb_vars%nrad_patch
     write(10,*) 'surfalb_vars%fabi_sun_patch' 
     write(10,*) surfalb_vars%fabi_sun_patch
     write(10,*) 'surfalb_vars%ncan_patch' 
     write(10,*) surfalb_vars%ncan_patch
     write(10,*) 'surfalb_vars%vcmaxcintsun_patch' 
     write(10,*) surfalb_vars%vcmaxcintsun_patch
     write(10,*) 'surfalb_vars%albgrd_col' 
     write(10,*) surfalb_vars%albgrd_col
     write(10,*) 'surfalb_vars%albgri_bc_col' 
     write(10,*) surfalb_vars%albgri_bc_col
     write(10,*) 'surfalb_vars%ftii_patch' 
     write(10,*) surfalb_vars%ftii_patch
     write(10,*) 'surfalb_vars%tlai_z_patch' 
     write(10,*) surfalb_vars%tlai_z_patch
     write(10,*) 'surfalb_vars%fabd_sha_patch' 
     write(10,*) surfalb_vars%fabd_sha_patch
     write(10,*) 'surfalb_vars%flx_absdn_col' 
     write(10,*) surfalb_vars%flx_absdn_col
     write(10,*) 'surfalb_vars%fabd_sha_z_patch' 
     write(10,*) surfalb_vars%fabd_sha_z_patch
     write(10,*) 'surfalb_vars%albgri_oc_col' 
     write(10,*) surfalb_vars%albgri_oc_col
     write(10,*) 'surfalb_vars%albgrd_oc_col' 
     write(10,*) surfalb_vars%albgrd_oc_col
     write(10,*) 'surfalb_vars%albsni_hst_col' 
     write(10,*) surfalb_vars%albsni_hst_col
     write(10,*) 'surfalb_vars%albgri_pur_col' 
     write(10,*) surfalb_vars%albgri_pur_col
     write(10,*) 'surfalb_vars%albgri_dst_col' 
     write(10,*) surfalb_vars%albgri_dst_col
     write(10,*) 'surfalb_vars%fsun_z_patch' 
     write(10,*) surfalb_vars%fsun_z_patch
     write(10,*) 'surfalb_vars%vcmaxcintsha_patch' 
     write(10,*) surfalb_vars%vcmaxcintsha_patch
     write(10,*) 'surfalb_vars%albgri_col' 
     write(10,*) surfalb_vars%albgri_col
     write(10,*) 'surfalb_vars%albi_patch' 
     write(10,*) surfalb_vars%albi_patch
     write(10,*) 'surfalb_vars%flx_absdv_col' 
     write(10,*) surfalb_vars%flx_absdv_col
     write(10,*) 'surfalb_vars%fabd_patch' 
     write(10,*) surfalb_vars%fabd_patch
     write(10,*) 'surfalb_vars%albsod_col' 
     write(10,*) surfalb_vars%albsod_col
     write(10,*) 'surfalb_vars%fabi_sha_patch' 
     write(10,*) surfalb_vars%fabi_sha_patch
     write(10,*) 'surfalb_vars%albgrd_bc_col' 
     write(10,*) surfalb_vars%albgrd_bc_col
     write(10,*) 'surfalb_vars%flx_absiv_col' 
     write(10,*) surfalb_vars%flx_absiv_col
     write(10,*) 'surfalb_vars%albsnd_hst_col' 
     write(10,*) surfalb_vars%albsnd_hst_col
     write(10,*) 'surfalb_vars%albsoi_col' 
     write(10,*) surfalb_vars%albsoi_col
     write(10,*) 'surfalb_vars%albsod_col' 
     write(10,*) surfalb_vars%albsod_col
     write(10,*) 'surfalb_vars%albsoi_col' 
     write(10,*) surfalb_vars%albsoi_col
     write(10,*) 'surfalb_vars%fabd_sun_patch' 
     write(10,*) surfalb_vars%fabd_sun_patch
     write(10,*) 'surfalb_vars%fsun_z_patch' 
     write(10,*) surfalb_vars%fsun_z_patch
     write(10,*) 'surfalb_vars%fabi_sun_z_patch' 
     write(10,*) surfalb_vars%fabi_sun_z_patch
     write(10,*) 'surfalb_vars%vcmaxcintsha_patch' 
     write(10,*) surfalb_vars%vcmaxcintsha_patch
     write(10,*) 'surfalb_vars%albi_patch' 
     write(10,*) surfalb_vars%albi_patch
     write(10,*) 'surfalb_vars%fabd_sun_z_patch' 
     write(10,*) surfalb_vars%fabd_sun_z_patch
     write(10,*) 'surfalb_vars%ftii_patch' 
     write(10,*) surfalb_vars%ftii_patch
     write(10,*) 'surfalb_vars%fabd_patch' 
     write(10,*) surfalb_vars%fabd_patch
     write(10,*) 'surfalb_vars%fabi_sha_patch' 
     write(10,*) surfalb_vars%fabi_sha_patch
     write(10,*) 'surfalb_vars%ftdd_patch' 
     write(10,*) surfalb_vars%ftdd_patch
     write(10,*) 'surfalb_vars%fabd_sha_z_patch' 
     write(10,*) surfalb_vars%fabd_sha_z_patch
     write(10,*) 'surfalb_vars%albd_patch' 
     write(10,*) surfalb_vars%albd_patch
     write(10,*) 'surfalb_vars%fabd_sha_patch' 
     write(10,*) surfalb_vars%fabd_sha_patch
     write(10,*) 'surfalb_vars%fabi_sun_patch' 
     write(10,*) surfalb_vars%fabi_sun_patch
     write(10,*) 'surfalb_vars%fabi_sha_z_patch' 
     write(10,*) surfalb_vars%fabi_sha_z_patch
     write(10,*) 'surfalb_vars%ftid_patch' 
     write(10,*) surfalb_vars%ftid_patch
     write(10,*) 'surfalb_vars%vcmaxcintsun_patch' 
     write(10,*) surfalb_vars%vcmaxcintsun_patch
     write(10,*) 'surfalb_vars%fabi_patch' 
     write(10,*) surfalb_vars%fabi_patch
     close(10)
end subroutine 
subroutine update_vars_UrbanAlbedo(gpu)
     use LandunitType, only : lun_pp 
     use UrbanParamsType, only : urbanparams_vars 
     use clm_instMod, only : solarabs_vars 
     use clm_instMod, only : surfalb_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_UrbanAlbedo.txt"
     else
          file='cpu_UrbanAlbedo.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc lun_pp%wtroad_perv, & 
     !$acc lun_pp%canyon_hwr )
     !$acc update self(& 
     !$acc urbanparams_vars%alb_wall_dif, & 
     !$acc urbanparams_vars%alb_wall_dir )
     !$acc update self(& 
     !$acc solarabs_vars%sabs_shadewall_dir_lun, & 
     !$acc solarabs_vars%sabs_sunwall_dir_lun, & 
     !$acc solarabs_vars%sabs_shadewall_dif_lun, & 
     !$acc solarabs_vars%sabs_improad_dir_lun, & 
     !$acc solarabs_vars%sabs_perroad_dif_lun, & 
     !$acc solarabs_vars%sabs_sunwall_dif_lun, & 
     !$acc solarabs_vars%sabs_perroad_dir_lun, & 
     !$acc solarabs_vars%sabs_improad_dif_lun, & 
     !$acc solarabs_vars%sabs_roof_dif_lun, & 
     !$acc solarabs_vars%sabs_roof_dir_lun, & 
     !$acc solarabs_vars%sabs_shadewall_dir_lun, & 
     !$acc solarabs_vars%sabs_sunwall_dir_lun, & 
     !$acc solarabs_vars%sabs_shadewall_dif_lun, & 
     !$acc solarabs_vars%sabs_improad_dir_lun, & 
     !$acc solarabs_vars%sabs_perroad_dif_lun, & 
     !$acc solarabs_vars%sabs_sunwall_dif_lun, & 
     !$acc solarabs_vars%sabs_perroad_dir_lun, & 
     !$acc solarabs_vars%sabs_improad_dif_lun, & 
     !$acc solarabs_vars%sabs_roof_dif_lun, & 
     !$acc solarabs_vars%sabs_roof_dir_lun )
     !$acc update self(& 
     !$acc surfalb_vars%albgrd_col, & 
     !$acc surfalb_vars%fabd_sun_patch, & 
     !$acc surfalb_vars%albgri_col, & 
     !$acc surfalb_vars%albi_patch, & 
     !$acc surfalb_vars%ftii_patch, & 
     !$acc surfalb_vars%fabd_patch, & 
     !$acc surfalb_vars%fabi_sha_patch, & 
     !$acc surfalb_vars%ftdd_patch, & 
     !$acc surfalb_vars%fabd_sha_patch, & 
     !$acc surfalb_vars%fabi_sun_patch, & 
     !$acc surfalb_vars%albd_patch, & 
     !$acc surfalb_vars%ftid_patch, & 
     !$acc surfalb_vars%fabi_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'lun_pp%wtroad_perv' 
     write(10,*) lun_pp%wtroad_perv
     write(10,*) 'lun_pp%canyon_hwr' 
     write(10,*) lun_pp%canyon_hwr
     write(10,*) 'urbanparams_vars%alb_wall_dif' 
     write(10,*) urbanparams_vars%alb_wall_dif
     write(10,*) 'urbanparams_vars%alb_wall_dir' 
     write(10,*) urbanparams_vars%alb_wall_dir
     write(10,*) 'solarabs_vars%sabs_shadewall_dir_lun' 
     write(10,*) solarabs_vars%sabs_shadewall_dir_lun
     write(10,*) 'solarabs_vars%sabs_sunwall_dir_lun' 
     write(10,*) solarabs_vars%sabs_sunwall_dir_lun
     write(10,*) 'solarabs_vars%sabs_shadewall_dif_lun' 
     write(10,*) solarabs_vars%sabs_shadewall_dif_lun
     write(10,*) 'solarabs_vars%sabs_improad_dir_lun' 
     write(10,*) solarabs_vars%sabs_improad_dir_lun
     write(10,*) 'solarabs_vars%sabs_perroad_dif_lun' 
     write(10,*) solarabs_vars%sabs_perroad_dif_lun
     write(10,*) 'solarabs_vars%sabs_sunwall_dif_lun' 
     write(10,*) solarabs_vars%sabs_sunwall_dif_lun
     write(10,*) 'solarabs_vars%sabs_perroad_dir_lun' 
     write(10,*) solarabs_vars%sabs_perroad_dir_lun
     write(10,*) 'solarabs_vars%sabs_improad_dif_lun' 
     write(10,*) solarabs_vars%sabs_improad_dif_lun
     write(10,*) 'solarabs_vars%sabs_roof_dif_lun' 
     write(10,*) solarabs_vars%sabs_roof_dif_lun
     write(10,*) 'solarabs_vars%sabs_roof_dir_lun' 
     write(10,*) solarabs_vars%sabs_roof_dir_lun
     write(10,*) 'solarabs_vars%sabs_shadewall_dir_lun' 
     write(10,*) solarabs_vars%sabs_shadewall_dir_lun
     write(10,*) 'solarabs_vars%sabs_sunwall_dir_lun' 
     write(10,*) solarabs_vars%sabs_sunwall_dir_lun
     write(10,*) 'solarabs_vars%sabs_shadewall_dif_lun' 
     write(10,*) solarabs_vars%sabs_shadewall_dif_lun
     write(10,*) 'solarabs_vars%sabs_improad_dir_lun' 
     write(10,*) solarabs_vars%sabs_improad_dir_lun
     write(10,*) 'solarabs_vars%sabs_perroad_dif_lun' 
     write(10,*) solarabs_vars%sabs_perroad_dif_lun
     write(10,*) 'solarabs_vars%sabs_sunwall_dif_lun' 
     write(10,*) solarabs_vars%sabs_sunwall_dif_lun
     write(10,*) 'solarabs_vars%sabs_perroad_dir_lun' 
     write(10,*) solarabs_vars%sabs_perroad_dir_lun
     write(10,*) 'solarabs_vars%sabs_improad_dif_lun' 
     write(10,*) solarabs_vars%sabs_improad_dif_lun
     write(10,*) 'solarabs_vars%sabs_roof_dif_lun' 
     write(10,*) solarabs_vars%sabs_roof_dif_lun
     write(10,*) 'solarabs_vars%sabs_roof_dir_lun' 
     write(10,*) solarabs_vars%sabs_roof_dir_lun
     write(10,*) 'surfalb_vars%albgrd_col' 
     write(10,*) surfalb_vars%albgrd_col
     write(10,*) 'surfalb_vars%fabd_sun_patch' 
     write(10,*) surfalb_vars%fabd_sun_patch
     write(10,*) 'surfalb_vars%albgri_col' 
     write(10,*) surfalb_vars%albgri_col
     write(10,*) 'surfalb_vars%albi_patch' 
     write(10,*) surfalb_vars%albi_patch
     write(10,*) 'surfalb_vars%ftii_patch' 
     write(10,*) surfalb_vars%ftii_patch
     write(10,*) 'surfalb_vars%fabd_patch' 
     write(10,*) surfalb_vars%fabd_patch
     write(10,*) 'surfalb_vars%fabi_sha_patch' 
     write(10,*) surfalb_vars%fabi_sha_patch
     write(10,*) 'surfalb_vars%ftdd_patch' 
     write(10,*) surfalb_vars%ftdd_patch
     write(10,*) 'surfalb_vars%fabd_sha_patch' 
     write(10,*) surfalb_vars%fabd_sha_patch
     write(10,*) 'surfalb_vars%fabi_sun_patch' 
     write(10,*) surfalb_vars%fabi_sun_patch
     write(10,*) 'surfalb_vars%albd_patch' 
     write(10,*) surfalb_vars%albd_patch
     write(10,*) 'surfalb_vars%ftid_patch' 
     write(10,*) surfalb_vars%ftid_patch
     write(10,*) 'surfalb_vars%fabi_patch' 
     write(10,*) surfalb_vars%fabi_patch
     close(10)
end subroutine 
subroutine update_vars_AnnualUpdate(gpu)
     use VegetationDataType, only : veg_cf 
     use clm_instMod, only : cnstate_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_AnnualUpdate.txt"
     else
          file='cpu_AnnualUpdate.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc veg_cf%annsum_npp, & 
     !$acc veg_cf%tempsum_npp )
     !$acc update self(& 
     !$acc cnstate_vars%tempavg_t2m_patch, & 
     !$acc cnstate_vars%annsum_counter_col, & 
     !$acc cnstate_vars%tempmax_retransp_patch, & 
     !$acc cnstate_vars%tempmax_retransn_patch, & 
     !$acc cnstate_vars%annsum_potential_gpp_patch, & 
     !$acc cnstate_vars%tempsum_potential_gpp_patch, & 
     !$acc cnstate_vars%annmax_retransn_patch, & 
     !$acc cnstate_vars%annmax_retransp_patch, & 
     !$acc cnstate_vars%annavg_t2m_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'veg_cf%annsum_npp' 
     write(10,*) veg_cf%annsum_npp
     write(10,*) 'veg_cf%tempsum_npp' 
     write(10,*) veg_cf%tempsum_npp
     write(10,*) 'cnstate_vars%tempavg_t2m_patch' 
     write(10,*) cnstate_vars%tempavg_t2m_patch
     write(10,*) 'cnstate_vars%annsum_counter_col' 
     write(10,*) cnstate_vars%annsum_counter_col
     write(10,*) 'cnstate_vars%tempmax_retransp_patch' 
     write(10,*) cnstate_vars%tempmax_retransp_patch
     write(10,*) 'cnstate_vars%tempmax_retransn_patch' 
     write(10,*) cnstate_vars%tempmax_retransn_patch
     write(10,*) 'cnstate_vars%annsum_potential_gpp_patch' 
     write(10,*) cnstate_vars%annsum_potential_gpp_patch
     write(10,*) 'cnstate_vars%tempsum_potential_gpp_patch' 
     write(10,*) cnstate_vars%tempsum_potential_gpp_patch
     write(10,*) 'cnstate_vars%annmax_retransn_patch' 
     write(10,*) cnstate_vars%annmax_retransn_patch
     write(10,*) 'cnstate_vars%annmax_retransp_patch' 
     write(10,*) cnstate_vars%annmax_retransp_patch
     write(10,*) 'cnstate_vars%annavg_t2m_patch' 
     write(10,*) cnstate_vars%annavg_t2m_patch
     close(10)
end subroutine 
subroutine update_vars_SatellitePhenology(gpu)
     use clm_instMod, only : canopystate_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SatellitePhenology.txt"
     else
          file='cpu_SatellitePhenology.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc canopystate_vars%esai_patch, & 
     !$acc canopystate_vars%hbot_patch, & 
     !$acc canopystate_vars%htop_patch, & 
     !$acc canopystate_vars%frac_veg_nosno_alb_patch, & 
     !$acc canopystate_vars%tlai_patch, & 
     !$acc canopystate_vars%tsai_patch, & 
     !$acc canopystate_vars%elai_patch, & 
     !$acc canopystate_vars%tlai_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'canopystate_vars%esai_patch' 
     write(10,*) canopystate_vars%esai_patch
     write(10,*) 'canopystate_vars%hbot_patch' 
     write(10,*) canopystate_vars%hbot_patch
     write(10,*) 'canopystate_vars%htop_patch' 
     write(10,*) canopystate_vars%htop_patch
     write(10,*) 'canopystate_vars%frac_veg_nosno_alb_patch' 
     write(10,*) canopystate_vars%frac_veg_nosno_alb_patch
     write(10,*) 'canopystate_vars%tlai_patch' 
     write(10,*) canopystate_vars%tlai_patch
     write(10,*) 'canopystate_vars%tsai_patch' 
     write(10,*) canopystate_vars%tsai_patch
     write(10,*) 'canopystate_vars%elai_patch' 
     write(10,*) canopystate_vars%elai_patch
     write(10,*) 'canopystate_vars%tlai_patch' 
     write(10,*) canopystate_vars%tlai_patch
     close(10)
end subroutine 
subroutine update_vars_depvel_compute(gpu)
     use clm_instMod, only : drydepvel_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_depvel_compute.txt"
     else
          file='cpu_depvel_compute.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc drydepvel_vars%velocity_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'drydepvel_vars%velocity_patch' 
#if 0
     write(10,*) drydepvel_vars%velocity_patch
#endif
     close(10)
end subroutine 
subroutine update_vars_CH4(gpu)
     use TopounitDataType, only : top_as 
     use clm_instMod, only : soilstate_vars 
     use clm_instMod, only : ch4_vars 
     use clm_instMod, only : lnd2atm_vars 
     use VegetationDataType, only : veg_cf 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_CH4.txt"
     else
          file='cpu_CH4.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc top_as%pch4bot )
     !$acc update self(& 
     !$acc soilstate_vars%rootfr_col )
     !$acc update self(& 
     !$acc ch4_vars%ch4prodg_grc, & 
     !$acc ch4_vars%finundated_lag_col, & 
     !$acc ch4_vars%fsat_bef_col, & 
     !$acc ch4_vars%ch4_surf_ebul_lake_col, & 
     !$acc ch4_vars%ch4_surf_flux_tot_col, & 
     !$acc ch4_vars%conc_o2_lake_col, & 
     !$acc ch4_vars%zwt_ch4_unsat_col, & 
     !$acc ch4_vars%c_atm_grc, & 
     !$acc ch4_vars%layer_sat_lag_col, & 
     !$acc ch4_vars%finundated_col, & 
     !$acc ch4_vars%grnd_ch4_cond_patch, & 
     !$acc ch4_vars%totcolch4_col, & 
     !$acc ch4_vars%qflx_surf_lag_col, & 
     !$acc ch4_vars%conc_ch4_sat_col, & 
     !$acc ch4_vars%ch4_surf_diff_lake_col, & 
     !$acc ch4_vars%ch4_dfsat_flux_col, & 
     !$acc ch4_vars%conc_ch4_lake_col, & 
     !$acc ch4_vars%ch4_oxid_depth_lake_col, & 
     !$acc ch4_vars%grnd_ch4_cond_col, & 
     !$acc ch4_vars%lake_soilc_col, & 
     !$acc ch4_vars%ch4co2f_grc, & 
     !$acc ch4_vars%ch4_prod_depth_lake_col, & 
     !$acc ch4_vars%annsum_counter_col, & 
     !$acc ch4_vars%tempavg_finrw_col, & 
     !$acc ch4_vars%tempavg_somhr_col, & 
     !$acc ch4_vars%annavg_somhr_col, & 
     !$acc ch4_vars%annavg_finrw_col, & 
     !$acc ch4_vars%ch4_prod_depth_sat_col, & 
     !$acc ch4_vars%sif_col, & 
     !$acc ch4_vars%ch4_prod_depth_unsat_col, & 
     !$acc ch4_vars%o2_decomp_depth_sat_col, & 
     !$acc ch4_vars%o2_decomp_depth_unsat_col, & 
     !$acc ch4_vars%o2_oxid_depth_sat_col, & 
     !$acc ch4_vars%o2_oxid_depth_unsat_col, & 
     !$acc ch4_vars%ch4_oxid_depth_sat_col, & 
     !$acc ch4_vars%ch4_oxid_depth_unsat_col, & 
     !$acc ch4_vars%ch4_aere_depth_sat_col, & 
     !$acc ch4_vars%ch4_tran_depth_unsat_col, & 
     !$acc ch4_vars%o2_aere_depth_unsat_col, & 
     !$acc ch4_vars%o2_aere_depth_sat_col, & 
     !$acc ch4_vars%ch4_tran_depth_sat_col, & 
     !$acc ch4_vars%ch4_aere_depth_unsat_col, & 
     !$acc ch4_vars%ch4_ebul_depth_sat_col, & 
     !$acc ch4_vars%ch4_ebul_depth_unsat_col, & 
     !$acc ch4_vars%ch4_surf_diff_unsat_col, & 
     !$acc ch4_vars%ch4_ebul_total_sat_col, & 
     !$acc ch4_vars%ch4_ebul_total_unsat_col, & 
     !$acc ch4_vars%ch4_surf_aere_unsat_col, & 
     !$acc ch4_vars%ch4_surf_ebul_unsat_col, & 
     !$acc ch4_vars%o2_decomp_depth_sat_col, & 
     !$acc ch4_vars%conc_o2_unsat_col, & 
     !$acc ch4_vars%ch4_oxid_depth_sat_col, & 
     !$acc ch4_vars%o2stress_sat_col, & 
     !$acc ch4_vars%ch4_surf_aere_sat_col, & 
     !$acc ch4_vars%o2stress_unsat_col, & 
     !$acc ch4_vars%ch4_aere_depth_unsat_col, & 
     !$acc ch4_vars%conc_ch4_unsat_col, & 
     !$acc ch4_vars%ch4_aere_depth_sat_col, & 
     !$acc ch4_vars%o2_oxid_depth_sat_col, & 
     !$acc ch4_vars%conc_o2_sat_col, & 
     !$acc ch4_vars%o2_oxid_depth_unsat_col, & 
     !$acc ch4_vars%ch4_surf_diff_sat_col, & 
     !$acc ch4_vars%conc_ch4_sat_col, & 
     !$acc ch4_vars%ch4_oxid_depth_unsat_col, & 
     !$acc ch4_vars%ch4_surf_ebul_sat_col, & 
     !$acc ch4_vars%grnd_ch4_cond_col, & 
     !$acc ch4_vars%ch4_ebul_depth_sat_col, & 
     !$acc ch4_vars%ch4stress_unsat_col, & 
     !$acc ch4_vars%ch4stress_sat_col, & 
     !$acc ch4_vars%o2_decomp_depth_unsat_col, & 
     !$acc ch4_vars%ch4_ebul_depth_unsat_col )
     !$acc update self(& 
     !$acc lnd2atm_vars%nem_grc )
     !$acc update self(& 
     !$acc veg_cf%annavg_agnpp, & 
     !$acc veg_cf%annavg_bgnpp, & 
     !$acc veg_cf%tempavg_agnpp, & 
     !$acc veg_cf%tempavg_bgnpp )
     end if 
     !! CPU print statements !! 
     write(10,*) 'top_as%pch4bot' 
     write(10,*) top_as%pch4bot
     write(10,*) 'soilstate_vars%rootfr_col' 
     write(10,*) soilstate_vars%rootfr_col
     write(10,*) 'ch4_vars%ch4prodg_grc' 
     write(10,*) ch4_vars%ch4prodg_grc
     write(10,*) 'ch4_vars%finundated_lag_col' 
     write(10,*) ch4_vars%finundated_lag_col
     write(10,*) 'ch4_vars%fsat_bef_col' 
     write(10,*) ch4_vars%fsat_bef_col
     write(10,*) 'ch4_vars%ch4_surf_ebul_lake_col' 
     write(10,*) ch4_vars%ch4_surf_ebul_lake_col
     write(10,*) 'ch4_vars%ch4_surf_flux_tot_col' 
     write(10,*) ch4_vars%ch4_surf_flux_tot_col
     write(10,*) 'ch4_vars%conc_o2_lake_col' 
     write(10,*) ch4_vars%conc_o2_lake_col
     write(10,*) 'ch4_vars%zwt_ch4_unsat_col' 
     write(10,*) ch4_vars%zwt_ch4_unsat_col
     write(10,*) 'ch4_vars%c_atm_grc' 
     write(10,*) ch4_vars%c_atm_grc
     write(10,*) 'ch4_vars%layer_sat_lag_col' 
     write(10,*) ch4_vars%layer_sat_lag_col
     write(10,*) 'ch4_vars%finundated_col' 
     write(10,*) ch4_vars%finundated_col
     write(10,*) 'ch4_vars%grnd_ch4_cond_patch' 
     write(10,*) ch4_vars%grnd_ch4_cond_patch
     write(10,*) 'ch4_vars%totcolch4_col' 
     write(10,*) ch4_vars%totcolch4_col
     write(10,*) 'ch4_vars%qflx_surf_lag_col' 
     write(10,*) ch4_vars%qflx_surf_lag_col
     write(10,*) 'ch4_vars%conc_ch4_sat_col' 
     write(10,*) ch4_vars%conc_ch4_sat_col
     write(10,*) 'ch4_vars%ch4_surf_diff_lake_col' 
     write(10,*) ch4_vars%ch4_surf_diff_lake_col
     write(10,*) 'ch4_vars%ch4_dfsat_flux_col' 
     write(10,*) ch4_vars%ch4_dfsat_flux_col
     write(10,*) 'ch4_vars%conc_ch4_lake_col' 
     write(10,*) ch4_vars%conc_ch4_lake_col
     write(10,*) 'ch4_vars%ch4_oxid_depth_lake_col' 
     write(10,*) ch4_vars%ch4_oxid_depth_lake_col
     write(10,*) 'ch4_vars%grnd_ch4_cond_col' 
     write(10,*) ch4_vars%grnd_ch4_cond_col
     write(10,*) 'ch4_vars%lake_soilc_col' 
     write(10,*) ch4_vars%lake_soilc_col
     write(10,*) 'ch4_vars%ch4co2f_grc' 
     write(10,*) ch4_vars%ch4co2f_grc
     write(10,*) 'ch4_vars%ch4_prod_depth_lake_col' 
     write(10,*) ch4_vars%ch4_prod_depth_lake_col
     write(10,*) 'ch4_vars%annsum_counter_col' 
     write(10,*) ch4_vars%annsum_counter_col
     write(10,*) 'ch4_vars%tempavg_finrw_col' 
     write(10,*) ch4_vars%tempavg_finrw_col
     write(10,*) 'ch4_vars%tempavg_somhr_col' 
     write(10,*) ch4_vars%tempavg_somhr_col
     write(10,*) 'ch4_vars%annavg_somhr_col' 
     write(10,*) ch4_vars%annavg_somhr_col
     write(10,*) 'ch4_vars%annavg_finrw_col' 
     write(10,*) ch4_vars%annavg_finrw_col
     write(10,*) 'ch4_vars%ch4_prod_depth_sat_col' 
     write(10,*) ch4_vars%ch4_prod_depth_sat_col
     write(10,*) 'ch4_vars%sif_col' 
     write(10,*) ch4_vars%sif_col
     write(10,*) 'ch4_vars%ch4_prod_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_prod_depth_unsat_col
     write(10,*) 'ch4_vars%o2_decomp_depth_sat_col' 
     write(10,*) ch4_vars%o2_decomp_depth_sat_col
     write(10,*) 'ch4_vars%o2_decomp_depth_unsat_col' 
     write(10,*) ch4_vars%o2_decomp_depth_unsat_col
     write(10,*) 'ch4_vars%o2_oxid_depth_sat_col' 
     write(10,*) ch4_vars%o2_oxid_depth_sat_col
     write(10,*) 'ch4_vars%o2_oxid_depth_unsat_col' 
     write(10,*) ch4_vars%o2_oxid_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_oxid_depth_sat_col' 
     write(10,*) ch4_vars%ch4_oxid_depth_sat_col
     write(10,*) 'ch4_vars%ch4_oxid_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_oxid_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_aere_depth_sat_col' 
     write(10,*) ch4_vars%ch4_aere_depth_sat_col
     write(10,*) 'ch4_vars%ch4_tran_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_tran_depth_unsat_col
     write(10,*) 'ch4_vars%o2_aere_depth_unsat_col' 
     write(10,*) ch4_vars%o2_aere_depth_unsat_col
     write(10,*) 'ch4_vars%o2_aere_depth_sat_col' 
     write(10,*) ch4_vars%o2_aere_depth_sat_col
     write(10,*) 'ch4_vars%ch4_tran_depth_sat_col' 
     write(10,*) ch4_vars%ch4_tran_depth_sat_col
     write(10,*) 'ch4_vars%ch4_aere_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_aere_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_ebul_depth_sat_col' 
     write(10,*) ch4_vars%ch4_ebul_depth_sat_col
     write(10,*) 'ch4_vars%ch4_ebul_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_ebul_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_surf_diff_unsat_col' 
     write(10,*) ch4_vars%ch4_surf_diff_unsat_col
     write(10,*) 'ch4_vars%ch4_ebul_total_sat_col' 
     write(10,*) ch4_vars%ch4_ebul_total_sat_col
     write(10,*) 'ch4_vars%ch4_ebul_total_unsat_col' 
     write(10,*) ch4_vars%ch4_ebul_total_unsat_col
     write(10,*) 'ch4_vars%ch4_surf_aere_unsat_col' 
     write(10,*) ch4_vars%ch4_surf_aere_unsat_col
     write(10,*) 'ch4_vars%ch4_surf_ebul_unsat_col' 
     write(10,*) ch4_vars%ch4_surf_ebul_unsat_col
     write(10,*) 'ch4_vars%o2_decomp_depth_sat_col' 
     write(10,*) ch4_vars%o2_decomp_depth_sat_col
     write(10,*) 'ch4_vars%conc_o2_unsat_col' 
     write(10,*) ch4_vars%conc_o2_unsat_col
     write(10,*) 'ch4_vars%ch4_oxid_depth_sat_col' 
     write(10,*) ch4_vars%ch4_oxid_depth_sat_col
     write(10,*) 'ch4_vars%o2stress_sat_col' 
     write(10,*) ch4_vars%o2stress_sat_col
     write(10,*) 'ch4_vars%ch4_surf_aere_sat_col' 
     write(10,*) ch4_vars%ch4_surf_aere_sat_col
     write(10,*) 'ch4_vars%o2stress_unsat_col' 
     write(10,*) ch4_vars%o2stress_unsat_col
     write(10,*) 'ch4_vars%ch4_aere_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_aere_depth_unsat_col
     write(10,*) 'ch4_vars%conc_ch4_unsat_col' 
     write(10,*) ch4_vars%conc_ch4_unsat_col
     write(10,*) 'ch4_vars%ch4_aere_depth_sat_col' 
     write(10,*) ch4_vars%ch4_aere_depth_sat_col
     write(10,*) 'ch4_vars%o2_oxid_depth_sat_col' 
     write(10,*) ch4_vars%o2_oxid_depth_sat_col
     write(10,*) 'ch4_vars%conc_o2_sat_col' 
     write(10,*) ch4_vars%conc_o2_sat_col
     write(10,*) 'ch4_vars%o2_oxid_depth_unsat_col' 
     write(10,*) ch4_vars%o2_oxid_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_surf_diff_sat_col' 
     write(10,*) ch4_vars%ch4_surf_diff_sat_col
     write(10,*) 'ch4_vars%conc_ch4_sat_col' 
     write(10,*) ch4_vars%conc_ch4_sat_col
     write(10,*) 'ch4_vars%ch4_oxid_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_oxid_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_surf_ebul_sat_col' 
     write(10,*) ch4_vars%ch4_surf_ebul_sat_col
     write(10,*) 'ch4_vars%grnd_ch4_cond_col' 
     write(10,*) ch4_vars%grnd_ch4_cond_col
     write(10,*) 'ch4_vars%ch4_ebul_depth_sat_col' 
     write(10,*) ch4_vars%ch4_ebul_depth_sat_col
     write(10,*) 'ch4_vars%ch4stress_unsat_col' 
     write(10,*) ch4_vars%ch4stress_unsat_col
     write(10,*) 'ch4_vars%ch4stress_sat_col' 
     write(10,*) ch4_vars%ch4stress_sat_col
     write(10,*) 'ch4_vars%o2_decomp_depth_unsat_col' 
     write(10,*) ch4_vars%o2_decomp_depth_unsat_col
     write(10,*) 'ch4_vars%ch4_ebul_depth_unsat_col' 
     write(10,*) ch4_vars%ch4_ebul_depth_unsat_col
     write(10,*) 'lnd2atm_vars%nem_grc' 
     write(10,*) lnd2atm_vars%nem_grc
     write(10,*) 'veg_cf%annavg_agnpp' 
     write(10,*) veg_cf%annavg_agnpp
     write(10,*) 'veg_cf%annavg_bgnpp' 
     write(10,*) veg_cf%annavg_bgnpp
     write(10,*) 'veg_cf%tempavg_agnpp' 
     write(10,*) veg_cf%tempavg_agnpp
     write(10,*) 'veg_cf%tempavg_bgnpp' 
     write(10,*) veg_cf%tempavg_bgnpp
     close(10)
end subroutine 
subroutine update_vars_HydrologyDrainage(gpu)
     use ColumnDataType, only : col_ws 
     use ColumnDataType, only : col_wf 
     use clm_instMod, only : soilhydrology_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_HydrologyDrainage.txt"
     else
          file='cpu_HydrologyDrainage.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%endwb, & 
     !$acc col_ws%h2osoi_liq_depth_intg, & 
     !$acc col_ws%h2osoi_vol, & 
     !$acc col_ws%h2osoi_ice_depth_intg, & 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osfc )
     !$acc update self(& 
     !$acc col_wf%qflx_infl, & 
     !$acc col_wf%qflx_glcice_frz, & 
     !$acc col_wf%qflx_irr_demand, & 
     !$acc col_wf%qflx_h2osfc_surf, & 
     !$acc col_wf%qflx_runoff, & 
     !$acc col_wf%qflx_drain_perched, & 
     !$acc col_wf%qflx_surf, & 
     !$acc col_wf%qflx_qrgwl, & 
     !$acc col_wf%qflx_runoff_u, & 
     !$acc col_wf%qflx_glcice, & 
     !$acc col_wf%qflx_drain, & 
     !$acc col_wf%qflx_rsub_sat, & 
     !$acc col_wf%qflx_runoff_r, & 
     !$acc col_wf%qflx_snwcp_ice, & 
     !$acc col_wf%qflx_drain_perched, & 
     !$acc col_wf%qflx_qrgwl, & 
     !$acc col_wf%qflx_drain, & 
     !$acc col_wf%qflx_rsub_sat )
     !$acc update self(& 
     !$acc soilhydrology_vars%moist_col, & 
     !$acc soilhydrology_vars%moist_vol_col, & 
     !$acc soilhydrology_vars%ice_col, & 
     !$acc soilhydrology_vars%frost_table_col, & 
     !$acc soilhydrology_vars%icefrac_col, & 
     !$acc soilhydrology_vars%zwt_perched_col, & 
     !$acc soilhydrology_vars%wa_col, & 
     !$acc soilhydrology_vars%zwt_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%endwb' 
     write(10,*) col_ws%endwb
     write(10,*) 'col_ws%h2osoi_liq_depth_intg' 
     write(10,*) col_ws%h2osoi_liq_depth_intg
     write(10,*) 'col_ws%h2osoi_vol' 
     write(10,*) col_ws%h2osoi_vol
     write(10,*) 'col_ws%h2osoi_ice_depth_intg' 
     write(10,*) col_ws%h2osoi_ice_depth_intg
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_wf%qflx_infl' 
     write(10,*) col_wf%qflx_infl
     write(10,*) 'col_wf%qflx_glcice_frz' 
     write(10,*) col_wf%qflx_glcice_frz
     write(10,*) 'col_wf%qflx_irr_demand' 
     write(10,*) col_wf%qflx_irr_demand
     write(10,*) 'col_wf%qflx_h2osfc_surf' 
     write(10,*) col_wf%qflx_h2osfc_surf
     write(10,*) 'col_wf%qflx_runoff' 
     write(10,*) col_wf%qflx_runoff
     write(10,*) 'col_wf%qflx_drain_perched' 
     write(10,*) col_wf%qflx_drain_perched
     write(10,*) 'col_wf%qflx_surf' 
     write(10,*) col_wf%qflx_surf
     write(10,*) 'col_wf%qflx_qrgwl' 
     write(10,*) col_wf%qflx_qrgwl
     write(10,*) 'col_wf%qflx_runoff_u' 
     write(10,*) col_wf%qflx_runoff_u
     write(10,*) 'col_wf%qflx_glcice' 
     write(10,*) col_wf%qflx_glcice
     write(10,*) 'col_wf%qflx_drain' 
     write(10,*) col_wf%qflx_drain
     write(10,*) 'col_wf%qflx_rsub_sat' 
     write(10,*) col_wf%qflx_rsub_sat
     write(10,*) 'col_wf%qflx_runoff_r' 
     write(10,*) col_wf%qflx_runoff_r
     write(10,*) 'col_wf%qflx_snwcp_ice' 
     write(10,*) col_wf%qflx_snwcp_ice
     write(10,*) 'col_wf%qflx_drain_perched' 
     write(10,*) col_wf%qflx_drain_perched
     write(10,*) 'col_wf%qflx_qrgwl' 
     write(10,*) col_wf%qflx_qrgwl
     write(10,*) 'col_wf%qflx_drain' 
     write(10,*) col_wf%qflx_drain
     write(10,*) 'col_wf%qflx_rsub_sat' 
     write(10,*) col_wf%qflx_rsub_sat
     write(10,*) 'soilhydrology_vars%moist_col' 
     write(10,*) soilhydrology_vars%moist_col
     write(10,*) 'soilhydrology_vars%moist_vol_col' 
     write(10,*) soilhydrology_vars%moist_vol_col
     write(10,*) 'soilhydrology_vars%ice_col' 
     write(10,*) soilhydrology_vars%ice_col
     write(10,*) 'soilhydrology_vars%frost_table_col' 
     write(10,*) soilhydrology_vars%frost_table_col
     write(10,*) 'soilhydrology_vars%icefrac_col' 
     write(10,*) soilhydrology_vars%icefrac_col
     write(10,*) 'soilhydrology_vars%zwt_perched_col' 
     write(10,*) soilhydrology_vars%zwt_perched_col
     write(10,*) 'soilhydrology_vars%wa_col' 
     write(10,*) soilhydrology_vars%wa_col
     write(10,*) 'soilhydrology_vars%zwt_col' 
     write(10,*) soilhydrology_vars%zwt_col
     close(10)
end subroutine 
subroutine update_vars_VegStructUpdate(gpu)
     use clm_instMod, only : cnstate_vars 
     use clm_instMod, only : canopystate_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_VegStructUpdate.txt"
     else
          file='cpu_VegStructUpdate.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc cnstate_vars%peaklai_patch, & 
     !$acc cnstate_vars%htmx_patch )
     !$acc update self(& 
     !$acc canopystate_vars%hbot_patch, & 
     !$acc canopystate_vars%esai_patch, & 
     !$acc canopystate_vars%htop_patch, & 
     !$acc canopystate_vars%frac_veg_nosno_alb_patch, & 
     !$acc canopystate_vars%tlai_patch, & 
     !$acc canopystate_vars%tsai_patch, & 
     !$acc canopystate_vars%elai_patch )
     end if 
     !! CPU print statements !! 
     write(10,*) 'cnstate_vars%peaklai_patch' 
     write(10,*) cnstate_vars%peaklai_patch
     write(10,*) 'cnstate_vars%htmx_patch' 
     write(10,*) cnstate_vars%htmx_patch
     write(10,*) 'canopystate_vars%hbot_patch' 
     write(10,*) canopystate_vars%hbot_patch
     write(10,*) 'canopystate_vars%esai_patch' 
     write(10,*) canopystate_vars%esai_patch
     write(10,*) 'canopystate_vars%htop_patch' 
     write(10,*) canopystate_vars%htop_patch
     write(10,*) 'canopystate_vars%frac_veg_nosno_alb_patch' 
     write(10,*) canopystate_vars%frac_veg_nosno_alb_patch
     write(10,*) 'canopystate_vars%tlai_patch' 
     write(10,*) canopystate_vars%tlai_patch
     write(10,*) 'canopystate_vars%tsai_patch' 
     write(10,*) canopystate_vars%tsai_patch
     write(10,*) 'canopystate_vars%elai_patch' 
     write(10,*) canopystate_vars%elai_patch
     close(10)
end subroutine 
subroutine update_vars_ColWaterBalanceCheck(gpu)
     use TopounitDataType, only : top_af 
     use ColumnDataType, only : col_ws 
     use ColumnDataType, only : col_wf 
     use ColumnDataType, only : col_ef 
     use VegetationDataType, only : veg_ef 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_ColWaterBalanceCheck.txt"
     else
          file='cpu_ColWaterBalanceCheck.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc top_af%rain, & 
     !$acc top_af%snow )
     !$acc update self(& 
     !$acc col_ws%errh2osno, & 
     !$acc col_ws%begwb, & 
     !$acc col_ws%total_plant_stored_h2o, & 
     !$acc col_ws%errh2o, & 
     !$acc col_ws%endwb )
     !$acc update self(& 
     !$acc col_wf%qflx_glcice_frz, & 
     !$acc col_wf%qflx_irrig, & 
     !$acc col_wf%qflx_over_supply, & 
     !$acc col_wf%qflx_h2osfc_surf, & 
     !$acc col_wf%qflx_snwcp_ice, & 
     !$acc col_wf%qflx_drain_perched, & 
     !$acc col_wf%qflx_surf, & 
     !$acc col_wf%snow_sources, & 
     !$acc col_wf%qflx_surf_irrig, & 
     !$acc col_wf%qflx_qrgwl, & 
     !$acc col_wf%qflx_lateral, & 
     !$acc col_wf%qflx_evap_tot, & 
     !$acc col_wf%dwb, & 
     !$acc col_wf%qflx_drain, & 
     !$acc col_wf%qflx_glcice_melt, & 
     !$acc col_wf%snow_sinks )
     !$acc update self(& 
     !$acc col_ef%errsoi )
     !$acc update self(& 
     !$acc veg_ef%errseb, & 
     !$acc veg_ef%errlon, & 
     !$acc veg_ef%errsol, & 
     !$acc veg_ef%netrad )
     end if 
     !! CPU print statements !! 
     write(10,*) 'top_af%rain' 
     write(10,*) top_af%rain
     write(10,*) 'top_af%snow' 
     write(10,*) top_af%snow
     write(10,*) 'col_ws%errh2osno' 
     write(10,*) col_ws%errh2osno
     write(10,*) 'col_ws%begwb' 
     write(10,*) col_ws%begwb
     write(10,*) 'col_ws%total_plant_stored_h2o' 
     write(10,*) col_ws%total_plant_stored_h2o
     write(10,*) 'col_ws%errh2o' 
     write(10,*) col_ws%errh2o
     write(10,*) 'col_ws%endwb' 
     write(10,*) col_ws%endwb
     write(10,*) 'col_wf%qflx_glcice_frz' 
     write(10,*) col_wf%qflx_glcice_frz
     write(10,*) 'col_wf%qflx_irrig' 
     write(10,*) col_wf%qflx_irrig
     write(10,*) 'col_wf%qflx_over_supply' 
     write(10,*) col_wf%qflx_over_supply
     write(10,*) 'col_wf%qflx_h2osfc_surf' 
     write(10,*) col_wf%qflx_h2osfc_surf
     write(10,*) 'col_wf%qflx_snwcp_ice' 
     write(10,*) col_wf%qflx_snwcp_ice
     write(10,*) 'col_wf%qflx_drain_perched' 
     write(10,*) col_wf%qflx_drain_perched
     write(10,*) 'col_wf%qflx_surf' 
     write(10,*) col_wf%qflx_surf
     write(10,*) 'col_wf%snow_sources' 
     write(10,*) col_wf%snow_sources
     write(10,*) 'col_wf%qflx_surf_irrig' 
     write(10,*) col_wf%qflx_surf_irrig
     write(10,*) 'col_wf%qflx_qrgwl' 
     write(10,*) col_wf%qflx_qrgwl
     write(10,*) 'col_wf%qflx_lateral' 
     write(10,*) col_wf%qflx_lateral
     write(10,*) 'col_wf%qflx_evap_tot' 
     write(10,*) col_wf%qflx_evap_tot
     write(10,*) 'col_wf%dwb' 
     write(10,*) col_wf%dwb
     write(10,*) 'col_wf%qflx_drain' 
     write(10,*) col_wf%qflx_drain
     write(10,*) 'col_wf%qflx_glcice_melt' 
     write(10,*) col_wf%qflx_glcice_melt
     write(10,*) 'col_wf%snow_sinks' 
     write(10,*) col_wf%snow_sinks
     write(10,*) 'col_ef%errsoi' 
     write(10,*) col_ef%errsoi
     write(10,*) 'veg_ef%errseb' 
     write(10,*) veg_ef%errseb
     write(10,*) 'veg_ef%errlon' 
     write(10,*) veg_ef%errlon
     write(10,*) 'veg_ef%errsol' 
     write(10,*) veg_ef%errsol
     write(10,*) 'veg_ef%netrad' 
     write(10,*) veg_ef%netrad
     close(10)
end subroutine 
subroutine update_vars_GridBalanceCheck(gpu)
     use ColumnDataType, only : col_ws 
     use GridcellDataType, only : grc_ws 
     use clm_instMod, only : soilhydrology_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_GridBalanceCheck.txt"
     else
          file='cpu_GridBalanceCheck.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%h2osoi_ice_depth_intg, & 
     !$acc col_ws%h2ocan, & 
     !$acc col_ws%h2osoi_liq_depth_intg, & 
     !$acc col_ws%h2osfc, & 
     !$acc col_ws%errh2o, & 
     !$acc col_ws%endwb, & 
     !$acc col_ws%h2osno )
     !$acc update self(& 
     !$acc grc_ws%end_h2ocan, & 
     !$acc grc_ws%end_h2osfc, & 
     !$acc grc_ws%end_h2osoi_ice, & 
     !$acc grc_ws%errh2o, & 
     !$acc grc_ws%end_h2osoi_liq, & 
     !$acc grc_ws%end_h2osno, & 
     !$acc grc_ws%endwb )
     !$acc update self(& 
     !$acc soilhydrology_vars%end_wa_grc )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%h2osoi_ice_depth_intg' 
     write(10,*) col_ws%h2osoi_ice_depth_intg
     write(10,*) 'col_ws%h2ocan' 
     write(10,*) col_ws%h2ocan
     write(10,*) 'col_ws%h2osoi_liq_depth_intg' 
     write(10,*) col_ws%h2osoi_liq_depth_intg
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_ws%errh2o' 
     write(10,*) col_ws%errh2o
     write(10,*) 'col_ws%endwb' 
     write(10,*) col_ws%endwb
     write(10,*) 'col_ws%h2osno' 
     write(10,*) col_ws%h2osno
     write(10,*) 'grc_ws%end_h2ocan' 
     write(10,*) grc_ws%end_h2ocan
     write(10,*) 'grc_ws%end_h2osfc' 
     write(10,*) grc_ws%end_h2osfc
     write(10,*) 'grc_ws%end_h2osoi_ice' 
     write(10,*) grc_ws%end_h2osoi_ice
     write(10,*) 'grc_ws%errh2o' 
     write(10,*) grc_ws%errh2o
     write(10,*) 'grc_ws%end_h2osoi_liq' 
     write(10,*) grc_ws%end_h2osoi_liq
     write(10,*) 'grc_ws%end_h2osno' 
     write(10,*) grc_ws%end_h2osno
     write(10,*) 'grc_ws%endwb' 
     write(10,*) grc_ws%endwb
     write(10,*) 'soilhydrology_vars%end_wa_grc' 
     write(10,*) soilhydrology_vars%end_wa_grc
     close(10)
end subroutine 
subroutine update_vars_ColCBalanceCheck(gpu)
     use ColumnDataType, only : col_cf 
     use ColumnDataType, only : col_cs 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_ColCBalanceCheck.txt"
     else
          file='cpu_ColCBalanceCheck.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_cf%er )
     !$acc update self(& 
     !$acc col_cs%endcb, & 
     !$acc col_cs%errcb )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_cf%er' 
     write(10,*) col_cf%er
     write(10,*) 'col_cs%endcb' 
     write(10,*) col_cs%endcb
     write(10,*) 'col_cs%errcb' 
     write(10,*) col_cs%errcb
     close(10)
end subroutine 
subroutine update_vars_ColNBalanceCheck(gpu)
     use ColumnDataType, only : col_nf 
     use ColumnDataType, only : col_ns 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_ColNBalanceCheck.txt"
     else
          file='cpu_ColNBalanceCheck.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_nf%ninputs, & 
     !$acc col_nf%noutputs )
     !$acc update self(& 
     !$acc col_ns%endnb, & 
     !$acc col_ns%errnb )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_nf%ninputs' 
     write(10,*) col_nf%ninputs
     write(10,*) 'col_nf%noutputs' 
     write(10,*) col_nf%noutputs
     write(10,*) 'col_ns%endnb' 
     write(10,*) col_ns%endnb
     write(10,*) 'col_ns%errnb' 
     write(10,*) col_ns%errnb
     close(10)
end subroutine 
subroutine update_vars_ColPBalanceCheck(gpu)
     use ColumnDataType, only : col_pf 
     use ColumnDataType, only : col_ps 
     use VegetationDataType, only : veg_pf 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_ColPBalanceCheck.txt"
     else
          file='cpu_ColPBalanceCheck.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pf%poutputs, & 
     !$acc col_pf%pinputs )
     !$acc update self(& 
     !$acc col_ps%errpb, & 
     !$acc col_ps%endpb )
     !$acc update self(& 
     !$acc veg_pf%leafp_to_litter, & 
     !$acc veg_pf%frootp_to_litter )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pf%poutputs' 
     write(10,*) col_pf%poutputs
     write(10,*) 'col_pf%pinputs' 
     write(10,*) col_pf%pinputs
     write(10,*) 'col_ps%errpb' 
     write(10,*) col_ps%errpb
     write(10,*) 'col_ps%endpb' 
     write(10,*) col_ps%endpb
     write(10,*) 'veg_pf%leafp_to_litter' 
     write(10,*) veg_pf%leafp_to_litter
     write(10,*) 'veg_pf%frootp_to_litter' 
     write(10,*) veg_pf%frootp_to_litter
     close(10)
end subroutine 
subroutine update_vars_WaterBudget_SetEndingMonthlyStates(gpu)
     use ColumnDataType, only : col_ws 
     use GridcellDataType, only : grc_ws 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_WaterBudget_SetEndingMonthlyStates.txt"
     else
          file='cpu_WaterBudget_SetEndingMonthlyStates.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%endwb )
     !$acc update self(& 
     !$acc grc_ws%tws_month_end )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%endwb' 
     write(10,*) col_ws%endwb
     write(10,*) 'grc_ws%tws_month_end' 
     write(10,*) grc_ws%tws_month_end
     close(10)
end subroutine 
subroutine update_vars_SoilErosion(gpu)
     use clm_instMod, only : sedflux_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_SoilErosion.txt"
     else
          file='cpu_SoilErosion.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc sedflux_vars%sed_p_ero_col, & 
     !$acc sedflux_vars%sed_ero_col, & 
     !$acc sedflux_vars%sed_q_ero_col, & 
     !$acc sedflux_vars%sed_yld_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'sedflux_vars%sed_p_ero_col' 
     write(10,*) sedflux_vars%sed_p_ero_col
     write(10,*) 'sedflux_vars%sed_ero_col' 
     write(10,*) sedflux_vars%sed_ero_col
     write(10,*) 'sedflux_vars%sed_q_ero_col' 
     write(10,*) sedflux_vars%sed_q_ero_col
     write(10,*) 'sedflux_vars%sed_yld_col' 
     write(10,*) sedflux_vars%sed_yld_col
     close(10)
end subroutine 
subroutine update_vars_EcosystemDynLeaching(gpu)
     use ColumnDataType, only : col_pf 
     use ColumnDataType, only : col_nf 
     use VegetationDataType, only : veg_ns 
     use ColumnDataType, only : col_ns 
     use VegetationDataType, only : veg_ps 
     use ColumnDataType, only : col_ps 
     use VegetationDataType, only : veg_cs 
     use clm_varctl        , only : use_c13, use_c14
     use VegetationDataType, only : c13_veg_cs 
     use ColumnDataType, only : c14_col_cs 
     use VegetationDataType, only : c14_veg_cs 
     use ColumnDataType, only : col_cs 
     use ColumnDataType, only : c13_col_cs 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_EcosystemDynLeaching.txt"
     else
          file='cpu_EcosystemDynLeaching.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pf%primp_to_labilep_vr, & 
     !$acc col_pf%labilep_to_secondp_vr, & 
     !$acc col_pf%secondp_to_labilep_vr, & 
     !$acc col_pf%secondp_to_occlp_vr, & 
     !$acc col_pf%biochem_pmin_vr, & 
     !$acc col_pf%biochem_pmin_ppools_vr, & 
     !$acc col_pf%sminp_leached_vr, & 
     !$acc col_pf%desorb_to_solutionp_vr, & 
     !$acc col_pf%adsorb_to_labilep_vr, & 
     !$acc col_pf%labilep_to_secondp_vr, & 
     !$acc col_pf%sminp_leached_vr )
     !$acc update self(& 
     !$acc col_nf%sminn_leached_vr, & 
     !$acc col_nf%smin_no3_leached_vr, & 
     !$acc col_nf%smin_no3_runoff_vr )
     !$acc update self(& 
     !$acc veg_ns%livestemn, & 
     !$acc veg_ns%deadcrootn_storage, & 
     !$acc veg_ns%deadstemn, & 
     !$acc veg_ns%npool, & 
     !$acc veg_ns%frootn, & 
     !$acc veg_ns%deadstemn_storage, & 
     !$acc veg_ns%livecrootn, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%livestemn_storage, & 
     !$acc veg_ns%deadcrootn, & 
     !$acc veg_ns%livecrootn_storage, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%frootn_storage, & 
     !$acc veg_ns%leafn_storage, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%retransn, & 
     !$acc veg_ns%leafn, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn, & 
     !$acc veg_ns%deadstemn, & 
     !$acc veg_ns%grainn_storage, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%retransn, & 
     !$acc veg_ns%npool, & 
     !$acc veg_ns%grainn_xfer, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn_storage, & 
     !$acc veg_ns%frootn, & 
     !$acc veg_ns%livecrootn, & 
     !$acc veg_ns%leafn, & 
     !$acc veg_ns%ntrunc, & 
     !$acc veg_ns%livecrootn_storage, & 
     !$acc veg_ns%frootn_storage, & 
     !$acc veg_ns%leafn_storage, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%deadstemn_storage, & 
     !$acc veg_ns%deadcrootn_storage, & 
     !$acc veg_ns%deadcrootn, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%grainn )
     !$acc update self(& 
     !$acc col_ns%decomp_npools_vr, & 
     !$acc col_ns%smin_no3_vr, & 
     !$acc col_ns%sminn_vr, & 
     !$acc col_ns%ntrunc_vr, & 
     !$acc col_ns%decomp_npools_vr, & 
     !$acc col_ns%smin_nh4_vr, & 
     !$acc col_ns%smin_no3_vr )
     !$acc update self(& 
     !$acc veg_ps%leafp, & 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%livecrootp, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%deadstemp_storage, & 
     !$acc veg_ps%deadcrootp, & 
     !$acc veg_ps%retransp, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%livestemp, & 
     !$acc veg_ps%leafp_storage, & 
     !$acc veg_ps%deadcrootp_storage, & 
     !$acc veg_ps%ppool, & 
     !$acc veg_ps%frootp_storage, & 
     !$acc veg_ps%frootp, & 
     !$acc veg_ps%deadstemp, & 
     !$acc veg_ps%livecrootp_storage, & 
     !$acc veg_ps%livestemp_storage, & 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%livecrootp, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%deadcrootp_storage, & 
     !$acc veg_ps%retransp, & 
     !$acc veg_ps%leafp_storage, & 
     !$acc veg_ps%deadstemp_storage, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%livecrootp_storage, & 
     !$acc veg_ps%ppool, & 
     !$acc veg_ps%leafp, & 
     !$acc veg_ps%grainp, & 
     !$acc veg_ps%livestemp_storage, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%grainp_storage, & 
     !$acc veg_ps%livestemp, & 
     !$acc veg_ps%grainp_xfer, & 
     !$acc veg_ps%ptrunc, & 
     !$acc veg_ps%deadcrootp, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%frootp_storage, & 
     !$acc veg_ps%frootp, & 
     !$acc veg_ps%deadstemp )
     !$acc update self(& 
     !$acc col_ps%labilep_vr, & 
     !$acc col_ps%decomp_ppools_vr, & 
     !$acc col_ps%primp_vr_prev, & 
     !$acc col_ps%primp_vr_cur, & 
     !$acc col_ps%secondp_vr_cur, & 
     !$acc col_ps%solutionp_vr_prev, & 
     !$acc col_ps%labilep_vr_prev, & 
     !$acc col_ps%occlp_vr_prev, & 
     !$acc col_ps%secondp_vr, & 
     !$acc col_ps%occlp_vr, & 
     !$acc col_ps%secondp_vr_prev, & 
     !$acc col_ps%occlp_vr_cur, & 
     !$acc col_ps%solutionp_vr, & 
     !$acc col_ps%labilep_vr_cur, & 
     !$acc col_ps%primp_vr, & 
     !$acc col_ps%solutionp_vr_cur, & 
     !$acc col_ps%decomp_ppools_vr, & 
     !$acc col_ps%ptrunc_vr )
     !$acc update self(& 
     !$acc veg_cs%cpool, & 
     !$acc veg_cs%gresp_xfer, & 
     !$acc veg_cs%deadstemc_xfer, & 
     !$acc veg_cs%grainc, & 
     !$acc veg_cs%leafc_storage, & 
     !$acc veg_cs%deadstemc, & 
     !$acc veg_cs%deadcrootc, & 
     !$acc veg_cs%livecrootc_storage, & 
     !$acc veg_cs%livestemc_storage, & 
     !$acc veg_cs%frootc_storage, & 
     !$acc veg_cs%livestemc_xfer, & 
     !$acc veg_cs%xsmrpool, & 
     !$acc veg_cs%grainc_storage, & 
     !$acc veg_cs%leafc_xfer, & 
     !$acc veg_cs%gresp_storage, & 
     !$acc veg_cs%ctrunc, & 
     !$acc veg_cs%grainc_xfer, & 
     !$acc veg_cs%frootc, & 
     !$acc veg_cs%livestemc, & 
     !$acc veg_cs%leafc, & 
     !$acc veg_cs%livecrootc_xfer, & 
     !$acc veg_cs%deadstemc_storage, & 
     !$acc veg_cs%livecrootc, & 
     !$acc veg_cs%deadcrootc_storage, & 
     !$acc veg_cs%frootc_xfer, & 
     !$acc veg_cs%deadcrootc_xfer )
     !$acc update self(& 
     !$acc c13_veg_cs%cpool, & 
     !$acc c13_veg_cs%gresp_xfer, & 
     !$acc c13_veg_cs%deadstemc_xfer, & 
     !$acc c13_veg_cs%leafc_storage, & 
     !$acc c13_veg_cs%deadstemc, & 
     !$acc c13_veg_cs%deadcrootc, & 
     !$acc c13_veg_cs%livecrootc_storage, & 
     !$acc c13_veg_cs%livestemc_storage, & 
     !$acc c13_veg_cs%frootc_storage, & 
     !$acc c13_veg_cs%livestemc_xfer, & 
     !$acc c13_veg_cs%leafc_xfer, & 
     !$acc c13_veg_cs%gresp_storage, & 
     !$acc c13_veg_cs%ctrunc, & 
     !$acc c13_veg_cs%frootc, & 
     !$acc c13_veg_cs%livestemc, & 
     !$acc c13_veg_cs%livecrootc_xfer, & 
     !$acc c13_veg_cs%leafc, & 
     !$acc c13_veg_cs%deadstemc_storage, & 
     !$acc c13_veg_cs%livecrootc, & 
     !$acc c13_veg_cs%deadcrootc_storage, & 
     !$acc c13_veg_cs%frootc_xfer, & 
     !$acc c13_veg_cs%deadcrootc_xfer )
     !$acc update self(& 
     !$acc c14_col_cs%ctrunc_vr, & 
     !$acc c14_col_cs%decomp_cpools_vr )
     !$acc update self(& 
     !$acc c14_veg_cs%cpool, & 
     !$acc c14_veg_cs%gresp_xfer, & 
     !$acc c14_veg_cs%deadstemc_xfer, & 
     !$acc c14_veg_cs%leafc_storage, & 
     !$acc c14_veg_cs%deadstemc, & 
     !$acc c14_veg_cs%deadcrootc, & 
     !$acc c14_veg_cs%livecrootc_storage, & 
     !$acc c14_veg_cs%livestemc_storage, & 
     !$acc c14_veg_cs%frootc_storage, & 
     !$acc c14_veg_cs%livestemc_xfer, & 
     !$acc c14_veg_cs%leafc_xfer, & 
     !$acc c14_veg_cs%gresp_storage, & 
     !$acc c14_veg_cs%ctrunc, & 
     !$acc c14_veg_cs%frootc, & 
     !$acc c14_veg_cs%livestemc, & 
     !$acc c14_veg_cs%livecrootc_xfer, & 
     !$acc c14_veg_cs%leafc, & 
     !$acc c14_veg_cs%deadstemc_storage, & 
     !$acc c14_veg_cs%livecrootc, & 
     !$acc c14_veg_cs%deadcrootc_storage, & 
     !$acc c14_veg_cs%frootc_xfer, & 
     !$acc c14_veg_cs%deadcrootc_xfer )
     !$acc update self(& 
     !$acc col_cs%ctrunc_vr, & 
     !$acc col_cs%decomp_cpools_vr )
     !$acc update self(& 
     !$acc c13_col_cs%ctrunc_vr, & 
     !$acc c13_col_cs%decomp_cpools_vr )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pf%primp_to_labilep_vr' 
     write(10,*) col_pf%primp_to_labilep_vr
     write(10,*) 'col_pf%labilep_to_secondp_vr' 
     write(10,*) col_pf%labilep_to_secondp_vr
     write(10,*) 'col_pf%secondp_to_labilep_vr' 
     write(10,*) col_pf%secondp_to_labilep_vr
     write(10,*) 'col_pf%secondp_to_occlp_vr' 
     write(10,*) col_pf%secondp_to_occlp_vr
     write(10,*) 'col_pf%biochem_pmin_vr' 
     write(10,*) col_pf%biochem_pmin_vr
     write(10,*) 'col_pf%biochem_pmin_ppools_vr' 
     write(10,*) col_pf%biochem_pmin_ppools_vr
     write(10,*) 'col_pf%sminp_leached_vr' 
     write(10,*) col_pf%sminp_leached_vr
     write(10,*) 'col_pf%desorb_to_solutionp_vr' 
     write(10,*) col_pf%desorb_to_solutionp_vr
     write(10,*) 'col_pf%adsorb_to_labilep_vr' 
     write(10,*) col_pf%adsorb_to_labilep_vr
     write(10,*) 'col_pf%labilep_to_secondp_vr' 
     write(10,*) col_pf%labilep_to_secondp_vr
     write(10,*) 'col_pf%sminp_leached_vr' 
     write(10,*) col_pf%sminp_leached_vr
     write(10,*) 'col_nf%sminn_leached_vr' 
     write(10,*) col_nf%sminn_leached_vr
     write(10,*) 'col_nf%smin_no3_leached_vr' 
     write(10,*) col_nf%smin_no3_leached_vr
     write(10,*) 'col_nf%smin_no3_runoff_vr' 
     write(10,*) col_nf%smin_no3_runoff_vr
     write(10,*) 'veg_ns%livestemn' 
     write(10,*) veg_ns%livestemn
     write(10,*) 'veg_ns%deadcrootn_storage' 
     write(10,*) veg_ns%deadcrootn_storage
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ns%npool' 
     write(10,*) veg_ns%npool
     write(10,*) 'veg_ns%frootn' 
     write(10,*) veg_ns%frootn
     write(10,*) 'veg_ns%deadstemn_storage' 
     write(10,*) veg_ns%deadstemn_storage
     write(10,*) 'veg_ns%livecrootn' 
     write(10,*) veg_ns%livecrootn
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%livestemn_storage' 
     write(10,*) veg_ns%livestemn_storage
     write(10,*) 'veg_ns%deadcrootn' 
     write(10,*) veg_ns%deadcrootn
     write(10,*) 'veg_ns%livecrootn_storage' 
     write(10,*) veg_ns%livecrootn_storage
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%frootn_storage' 
     write(10,*) veg_ns%frootn_storage
     write(10,*) 'veg_ns%leafn_storage' 
     write(10,*) veg_ns%leafn_storage
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%retransn' 
     write(10,*) veg_ns%retransn
     write(10,*) 'veg_ns%leafn' 
     write(10,*) veg_ns%leafn
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn' 
     write(10,*) veg_ns%livestemn
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ns%grainn_storage' 
     write(10,*) veg_ns%grainn_storage
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%retransn' 
     write(10,*) veg_ns%retransn
     write(10,*) 'veg_ns%npool' 
     write(10,*) veg_ns%npool
     write(10,*) 'veg_ns%grainn_xfer' 
     write(10,*) veg_ns%grainn_xfer
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn_storage' 
     write(10,*) veg_ns%livestemn_storage
     write(10,*) 'veg_ns%frootn' 
     write(10,*) veg_ns%frootn
     write(10,*) 'veg_ns%livecrootn' 
     write(10,*) veg_ns%livecrootn
     write(10,*) 'veg_ns%leafn' 
     write(10,*) veg_ns%leafn
     write(10,*) 'veg_ns%ntrunc' 
     write(10,*) veg_ns%ntrunc
     write(10,*) 'veg_ns%livecrootn_storage' 
     write(10,*) veg_ns%livecrootn_storage
     write(10,*) 'veg_ns%frootn_storage' 
     write(10,*) veg_ns%frootn_storage
     write(10,*) 'veg_ns%leafn_storage' 
     write(10,*) veg_ns%leafn_storage
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%deadstemn_storage' 
     write(10,*) veg_ns%deadstemn_storage
     write(10,*) 'veg_ns%deadcrootn_storage' 
     write(10,*) veg_ns%deadcrootn_storage
     write(10,*) 'veg_ns%deadcrootn' 
     write(10,*) veg_ns%deadcrootn
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%grainn' 
     write(10,*) veg_ns%grainn
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'col_ns%smin_no3_vr' 
     write(10,*) col_ns%smin_no3_vr
     write(10,*) 'col_ns%sminn_vr' 
     write(10,*) col_ns%sminn_vr
     write(10,*) 'col_ns%ntrunc_vr' 
     write(10,*) col_ns%ntrunc_vr
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'col_ns%smin_nh4_vr' 
     write(10,*) col_ns%smin_nh4_vr
     write(10,*) 'col_ns%smin_no3_vr' 
     write(10,*) col_ns%smin_no3_vr
     write(10,*) 'veg_ps%leafp' 
     write(10,*) veg_ps%leafp
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%livecrootp' 
     write(10,*) veg_ps%livecrootp
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%deadstemp_storage' 
     write(10,*) veg_ps%deadstemp_storage
     write(10,*) 'veg_ps%deadcrootp' 
     write(10,*) veg_ps%deadcrootp
     write(10,*) 'veg_ps%retransp' 
     write(10,*) veg_ps%retransp
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%livestemp' 
     write(10,*) veg_ps%livestemp
     write(10,*) 'veg_ps%leafp_storage' 
     write(10,*) veg_ps%leafp_storage
     write(10,*) 'veg_ps%deadcrootp_storage' 
     write(10,*) veg_ps%deadcrootp_storage
     write(10,*) 'veg_ps%ppool' 
     write(10,*) veg_ps%ppool
     write(10,*) 'veg_ps%frootp_storage' 
     write(10,*) veg_ps%frootp_storage
     write(10,*) 'veg_ps%frootp' 
     write(10,*) veg_ps%frootp
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'veg_ps%livecrootp_storage' 
     write(10,*) veg_ps%livecrootp_storage
     write(10,*) 'veg_ps%livestemp_storage' 
     write(10,*) veg_ps%livestemp_storage
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%livecrootp' 
     write(10,*) veg_ps%livecrootp
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%deadcrootp_storage' 
     write(10,*) veg_ps%deadcrootp_storage
     write(10,*) 'veg_ps%retransp' 
     write(10,*) veg_ps%retransp
     write(10,*) 'veg_ps%leafp_storage' 
     write(10,*) veg_ps%leafp_storage
     write(10,*) 'veg_ps%deadstemp_storage' 
     write(10,*) veg_ps%deadstemp_storage
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%livecrootp_storage' 
     write(10,*) veg_ps%livecrootp_storage
     write(10,*) 'veg_ps%ppool' 
     write(10,*) veg_ps%ppool
     write(10,*) 'veg_ps%leafp' 
     write(10,*) veg_ps%leafp
     write(10,*) 'veg_ps%grainp' 
     write(10,*) veg_ps%grainp
     write(10,*) 'veg_ps%livestemp_storage' 
     write(10,*) veg_ps%livestemp_storage
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%grainp_storage' 
     write(10,*) veg_ps%grainp_storage
     write(10,*) 'veg_ps%livestemp' 
     write(10,*) veg_ps%livestemp
     write(10,*) 'veg_ps%grainp_xfer' 
     write(10,*) veg_ps%grainp_xfer
     write(10,*) 'veg_ps%ptrunc' 
     write(10,*) veg_ps%ptrunc
     write(10,*) 'veg_ps%deadcrootp' 
     write(10,*) veg_ps%deadcrootp
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%frootp_storage' 
     write(10,*) veg_ps%frootp_storage
     write(10,*) 'veg_ps%frootp' 
     write(10,*) veg_ps%frootp
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'col_ps%labilep_vr' 
     write(10,*) col_ps%labilep_vr
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_ps%primp_vr_prev' 
     write(10,*) col_ps%primp_vr_prev
     write(10,*) 'col_ps%primp_vr_cur' 
     write(10,*) col_ps%primp_vr_cur
     write(10,*) 'col_ps%secondp_vr_cur' 
     write(10,*) col_ps%secondp_vr_cur
     write(10,*) 'col_ps%solutionp_vr_prev' 
     write(10,*) col_ps%solutionp_vr_prev
     write(10,*) 'col_ps%labilep_vr_prev' 
     write(10,*) col_ps%labilep_vr_prev
     write(10,*) 'col_ps%occlp_vr_prev' 
     write(10,*) col_ps%occlp_vr_prev
     write(10,*) 'col_ps%secondp_vr' 
     write(10,*) col_ps%secondp_vr
     write(10,*) 'col_ps%occlp_vr' 
     write(10,*) col_ps%occlp_vr
     write(10,*) 'col_ps%secondp_vr_prev' 
     write(10,*) col_ps%secondp_vr_prev
     write(10,*) 'col_ps%occlp_vr_cur' 
     write(10,*) col_ps%occlp_vr_cur
     write(10,*) 'col_ps%solutionp_vr' 
     write(10,*) col_ps%solutionp_vr
     write(10,*) 'col_ps%labilep_vr_cur' 
     write(10,*) col_ps%labilep_vr_cur
     write(10,*) 'col_ps%primp_vr' 
     write(10,*) col_ps%primp_vr
     write(10,*) 'col_ps%solutionp_vr_cur' 
     write(10,*) col_ps%solutionp_vr_cur
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_ps%ptrunc_vr' 
     write(10,*) col_ps%ptrunc_vr
     write(10,*) 'veg_cs%cpool' 
     write(10,*) veg_cs%cpool
     write(10,*) 'veg_cs%gresp_xfer' 
     write(10,*) veg_cs%gresp_xfer
     write(10,*) 'veg_cs%deadstemc_xfer' 
     write(10,*) veg_cs%deadstemc_xfer
     write(10,*) 'veg_cs%grainc' 
     write(10,*) veg_cs%grainc
     write(10,*) 'veg_cs%leafc_storage' 
     write(10,*) veg_cs%leafc_storage
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     write(10,*) 'veg_cs%deadcrootc' 
     write(10,*) veg_cs%deadcrootc
     write(10,*) 'veg_cs%livecrootc_storage' 
     write(10,*) veg_cs%livecrootc_storage
     write(10,*) 'veg_cs%livestemc_storage' 
     write(10,*) veg_cs%livestemc_storage
     write(10,*) 'veg_cs%frootc_storage' 
     write(10,*) veg_cs%frootc_storage
     write(10,*) 'veg_cs%livestemc_xfer' 
     write(10,*) veg_cs%livestemc_xfer
     write(10,*) 'veg_cs%xsmrpool' 
     write(10,*) veg_cs%xsmrpool
     write(10,*) 'veg_cs%grainc_storage' 
     write(10,*) veg_cs%grainc_storage
     write(10,*) 'veg_cs%leafc_xfer' 
     write(10,*) veg_cs%leafc_xfer
     write(10,*) 'veg_cs%gresp_storage' 
     write(10,*) veg_cs%gresp_storage
     write(10,*) 'veg_cs%ctrunc' 
     write(10,*) veg_cs%ctrunc
     write(10,*) 'veg_cs%grainc_xfer' 
     write(10,*) veg_cs%grainc_xfer
     write(10,*) 'veg_cs%frootc' 
     write(10,*) veg_cs%frootc
     write(10,*) 'veg_cs%livestemc' 
     write(10,*) veg_cs%livestemc
     write(10,*) 'veg_cs%leafc' 
     write(10,*) veg_cs%leafc
     write(10,*) 'veg_cs%livecrootc_xfer' 
     write(10,*) veg_cs%livecrootc_xfer
     write(10,*) 'veg_cs%deadstemc_storage' 
     write(10,*) veg_cs%deadstemc_storage
     write(10,*) 'veg_cs%livecrootc' 
     write(10,*) veg_cs%livecrootc
     write(10,*) 'veg_cs%deadcrootc_storage' 
     write(10,*) veg_cs%deadcrootc_storage
     write(10,*) 'veg_cs%frootc_xfer' 
     write(10,*) veg_cs%frootc_xfer
     write(10,*) 'veg_cs%deadcrootc_xfer' 
     write(10,*) veg_cs%deadcrootc_xfer
     if ( use_c13 ) then
     write(10,*) 'c13_veg_cs%cpool' 
     write(10,*) c13_veg_cs%cpool
     write(10,*) 'c13_veg_cs%gresp_xfer' 
     write(10,*) c13_veg_cs%gresp_xfer
     write(10,*) 'c13_veg_cs%deadstemc_xfer' 
     write(10,*) c13_veg_cs%deadstemc_xfer
     write(10,*) 'c13_veg_cs%leafc_storage' 
     write(10,*) c13_veg_cs%leafc_storage
     write(10,*) 'c13_veg_cs%deadstemc' 
     write(10,*) c13_veg_cs%deadstemc
     write(10,*) 'c13_veg_cs%deadcrootc' 
     write(10,*) c13_veg_cs%deadcrootc
     write(10,*) 'c13_veg_cs%livecrootc_storage' 
     write(10,*) c13_veg_cs%livecrootc_storage
     write(10,*) 'c13_veg_cs%livestemc_storage' 
     write(10,*) c13_veg_cs%livestemc_storage
     write(10,*) 'c13_veg_cs%frootc_storage' 
     write(10,*) c13_veg_cs%frootc_storage
     write(10,*) 'c13_veg_cs%livestemc_xfer' 
     write(10,*) c13_veg_cs%livestemc_xfer
     write(10,*) 'c13_veg_cs%leafc_xfer' 
     write(10,*) c13_veg_cs%leafc_xfer
     write(10,*) 'c13_veg_cs%gresp_storage' 
     write(10,*) c13_veg_cs%gresp_storage
     write(10,*) 'c13_veg_cs%ctrunc' 
     write(10,*) c13_veg_cs%ctrunc
     write(10,*) 'c13_veg_cs%frootc' 
     write(10,*) c13_veg_cs%frootc
     write(10,*) 'c13_veg_cs%livestemc' 
     write(10,*) c13_veg_cs%livestemc
     write(10,*) 'c13_veg_cs%livecrootc_xfer' 
     write(10,*) c13_veg_cs%livecrootc_xfer
     write(10,*) 'c13_veg_cs%leafc' 
     write(10,*) c13_veg_cs%leafc
     write(10,*) 'c13_veg_cs%deadstemc_storage' 
     write(10,*) c13_veg_cs%deadstemc_storage
     write(10,*) 'c13_veg_cs%livecrootc' 
     write(10,*) c13_veg_cs%livecrootc
     write(10,*) 'c13_veg_cs%deadcrootc_storage' 
     write(10,*) c13_veg_cs%deadcrootc_storage
     write(10,*) 'c13_veg_cs%frootc_xfer' 
     write(10,*) c13_veg_cs%frootc_xfer
     write(10,*) 'c13_veg_cs%deadcrootc_xfer' 
     write(10,*) c13_veg_cs%deadcrootc_xfer
     endif
     if ( use_c14 ) then
     write(10,*) 'c14_col_cs%ctrunc_vr' 
     write(10,*) c14_col_cs%ctrunc_vr
     write(10,*) 'c14_col_cs%decomp_cpools_vr' 
     write(10,*) c14_col_cs%decomp_cpools_vr
     write(10,*) 'c14_veg_cs%cpool' 
     write(10,*) c14_veg_cs%cpool
     write(10,*) 'c14_veg_cs%gresp_xfer' 
     write(10,*) c14_veg_cs%gresp_xfer
     write(10,*) 'c14_veg_cs%deadstemc_xfer' 
     write(10,*) c14_veg_cs%deadstemc_xfer
     write(10,*) 'c14_veg_cs%leafc_storage' 
     write(10,*) c14_veg_cs%leafc_storage
     write(10,*) 'c14_veg_cs%deadstemc' 
     write(10,*) c14_veg_cs%deadstemc
     write(10,*) 'c14_veg_cs%deadcrootc' 
     write(10,*) c14_veg_cs%deadcrootc
     write(10,*) 'c14_veg_cs%livecrootc_storage' 
     write(10,*) c14_veg_cs%livecrootc_storage
     write(10,*) 'c14_veg_cs%livestemc_storage' 
     write(10,*) c14_veg_cs%livestemc_storage
     write(10,*) 'c14_veg_cs%frootc_storage' 
     write(10,*) c14_veg_cs%frootc_storage
     write(10,*) 'c14_veg_cs%livestemc_xfer' 
     write(10,*) c14_veg_cs%livestemc_xfer
     write(10,*) 'c14_veg_cs%leafc_xfer' 
     write(10,*) c14_veg_cs%leafc_xfer
     write(10,*) 'c14_veg_cs%gresp_storage' 
     write(10,*) c14_veg_cs%gresp_storage
     write(10,*) 'c14_veg_cs%ctrunc' 
     write(10,*) c14_veg_cs%ctrunc
     write(10,*) 'c14_veg_cs%frootc' 
     write(10,*) c14_veg_cs%frootc
     write(10,*) 'c14_veg_cs%livestemc' 
     write(10,*) c14_veg_cs%livestemc
     write(10,*) 'c14_veg_cs%livecrootc_xfer' 
     write(10,*) c14_veg_cs%livecrootc_xfer
     write(10,*) 'c14_veg_cs%leafc' 
     write(10,*) c14_veg_cs%leafc
     write(10,*) 'c14_veg_cs%deadstemc_storage' 
     write(10,*) c14_veg_cs%deadstemc_storage
     write(10,*) 'c14_veg_cs%livecrootc' 
     write(10,*) c14_veg_cs%livecrootc
     write(10,*) 'c14_veg_cs%deadcrootc_storage' 
     write(10,*) c14_veg_cs%deadcrootc_storage
     write(10,*) 'c14_veg_cs%frootc_xfer' 
     write(10,*) c14_veg_cs%frootc_xfer
     write(10,*) 'c14_veg_cs%deadcrootc_xfer' 
     write(10,*) c14_veg_cs%deadcrootc_xfer
     endif
     write(10,*) 'col_cs%ctrunc_vr' 
     write(10,*) col_cs%ctrunc_vr
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     if ( use_c13 ) then
     write(10,*) 'c13_col_cs%ctrunc_vr' 
     write(10,*) c13_col_cs%ctrunc_vr
     write(10,*) 'c13_col_cs%decomp_cpools_vr' 
     write(10,*) c13_col_cs%decomp_cpools_vr
     endif
     close(10)
end subroutine 
subroutine update_vars_dyn_hwcontent_init(gpu)
     use VegetationDataType, only : veg_ws 
     use ColumnDataType, only : col_es 
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : soilhydrology_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_dyn_hwcontent_init.txt"
     else
          file='cpu_dyn_hwcontent_init.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc veg_ws%h2ocan, & 
     !$acc veg_ws%h2ocan )
     !$acc update self(& 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_h2osfc, & 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osfc, & 
     !$acc col_ws%h2osoi_liq )
     !$acc update self(& 
     !$acc soilhydrology_vars%wa_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'veg_ws%h2ocan' 
     write(10,*) veg_ws%h2ocan
     write(10,*) 'veg_ws%h2ocan' 
     write(10,*) veg_ws%h2ocan
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_h2osfc' 
     write(10,*) col_es%t_h2osfc
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'soilhydrology_vars%wa_col' 
     write(10,*) soilhydrology_vars%wa_col
     close(10)
end subroutine 
subroutine update_vars_dynSubgrid_wrapup_weight_changes(gpu)
     use LandunitType, only : lun_pp 
     use ColumnType, only : col_pp 
     use VegetationType, only : veg_pp 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_dynSubgrid_wrapup_weight_changes.txt"
     else
          file='cpu_dynSubgrid_wrapup_weight_changes.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc lun_pp%wttopounit, & 
     !$acc lun_pp%wtgcell, & 
     !$acc lun_pp%active, & 
     !$acc lun_pp%active, & 
     !$acc lun_pp%wtgcell )
     !$acc update self(& 
     !$acc col_pp%wttopounit, & 
     !$acc col_pp%wtgcell, & 
     !$acc col_pp%active, & 
     !$acc col_pp%active, & 
     !$acc col_pp%wtgcell )
     !$acc update self(& 
     !$acc veg_pp%wtlunit, & 
     !$acc veg_pp%wttopounit, & 
     !$acc veg_pp%wtgcell, & 
     !$acc veg_pp%active, & 
     !$acc veg_pp%active, & 
     !$acc veg_pp%wtlunit, & 
     !$acc veg_pp%wtgcell )
     end if 
     !! CPU print statements !! 
     write(10,*) 'lun_pp%wttopounit' 
     write(10,*) lun_pp%wttopounit
     write(10,*) 'lun_pp%wtgcell' 
     write(10,*) lun_pp%wtgcell
     write(10,*) 'lun_pp%active' 
     write(10,*) lun_pp%active
     write(10,*) 'lun_pp%active' 
     write(10,*) lun_pp%active
     write(10,*) 'lun_pp%wtgcell' 
     write(10,*) lun_pp%wtgcell
     write(10,*) 'col_pp%wttopounit' 
     write(10,*) col_pp%wttopounit
     write(10,*) 'col_pp%wtgcell' 
     write(10,*) col_pp%wtgcell
     write(10,*) 'col_pp%active' 
     write(10,*) col_pp%active
     write(10,*) 'col_pp%active' 
     write(10,*) col_pp%active
     write(10,*) 'col_pp%wtgcell' 
     write(10,*) col_pp%wtgcell
     write(10,*) 'veg_pp%wtlunit' 
     write(10,*) veg_pp%wtlunit
     write(10,*) 'veg_pp%wttopounit' 
     write(10,*) veg_pp%wttopounit
     write(10,*) 'veg_pp%wtgcell' 
     write(10,*) veg_pp%wtgcell
     write(10,*) 'veg_pp%active' 
     write(10,*) veg_pp%active
     write(10,*) 'veg_pp%active' 
     write(10,*) veg_pp%active
     write(10,*) 'veg_pp%wtlunit' 
     write(10,*) veg_pp%wtlunit
     write(10,*) 'veg_pp%wtgcell' 
     write(10,*) veg_pp%wtgcell
     close(10)
end subroutine 
subroutine update_vars_check_weights(gpu)
     use ColumnType, only : col_pp 
     use LandunitType, only : lun_pp 
     use VegetationType, only : veg_pp 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_check_weights.txt"
     else
          file='cpu_check_weights.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_pp%active, & 
     !$acc col_pp%wtgcell )
     !$acc update self(& 
     !$acc lun_pp%active, & 
     !$acc lun_pp%wtgcell )
     !$acc update self(& 
     !$acc veg_pp%active, & 
     !$acc veg_pp%wtlunit, & 
     !$acc veg_pp%wtgcell )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_pp%active' 
     write(10,*) col_pp%active
     write(10,*) 'col_pp%wtgcell' 
     write(10,*) col_pp%wtgcell
     write(10,*) 'lun_pp%active' 
     write(10,*) lun_pp%active
     write(10,*) 'lun_pp%wtgcell' 
     write(10,*) lun_pp%wtgcell
     write(10,*) 'veg_pp%active' 
     write(10,*) veg_pp%active
     write(10,*) 'veg_pp%wtlunit' 
     write(10,*) veg_pp%wtlunit
     write(10,*) 'veg_pp%wtgcell' 
     write(10,*) veg_pp%wtgcell
     close(10)
end subroutine 
subroutine update_vars_setFilters(gpu)
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_setFilters.txt"
     else
          file='cpu_setFilters.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     end if 
     !! CPU print statements !! 
     close(10)
end subroutine 
subroutine update_vars_initialize_new_columns(gpu)
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : soilhydrology_vars 
     use ColumnDataType, only : col_es 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_initialize_new_columns.txt"
     else
          file='cpu_initialize_new_columns.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osoi_ice, & 
     !$acc col_ws%h2osoi_vol )
     !$acc update self(& 
     !$acc soilhydrology_vars%wa_col )
     !$acc update self(& 
     !$acc col_es%t_soisno )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osoi_ice' 
     write(10,*) col_ws%h2osoi_ice
     write(10,*) 'col_ws%h2osoi_vol' 
     write(10,*) col_ws%h2osoi_vol
     write(10,*) 'soilhydrology_vars%wa_col' 
     write(10,*) soilhydrology_vars%wa_col
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     close(10)
end subroutine 
subroutine update_vars_dyn_hwcontent_final(gpu)
     use GridcellDataType, only : grc_wf 
     use GridcellDataType, only : grc_ef 
     use GridcellDataType, only : grc_ws 
     use GridcellDataType, only : grc_es 
     use VegetationDataType, only : veg_ws 
     use ColumnDataType, only : col_es 
     use ColumnDataType, only : col_ws 
     use clm_instMod, only : soilhydrology_vars 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_dyn_hwcontent_final.txt"
     else
          file='cpu_dyn_hwcontent_final.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc grc_wf%qflx_liq_dynbal, & 
     !$acc grc_wf%qflx_ice_dynbal )
     !$acc update self(& 
     !$acc grc_ef%eflx_dynbal )
     !$acc update self(& 
     !$acc grc_ws%liq1, & 
     !$acc grc_ws%ice1, & 
     !$acc grc_ws%liq2, & 
     !$acc grc_ws%ice2 )
     !$acc update self(& 
     !$acc grc_es%heat2, & 
     !$acc grc_es%heat1 )
     !$acc update self(& 
     !$acc veg_ws%h2ocan, & 
     !$acc veg_ws%h2ocan )
     !$acc update self(& 
     !$acc col_es%t_soisno, & 
     !$acc col_es%t_h2osfc, & 
     !$acc col_es%t_soisno )
     !$acc update self(& 
     !$acc col_ws%h2osoi_liq, & 
     !$acc col_ws%h2osfc, & 
     !$acc col_ws%h2osoi_liq )
     !$acc update self(& 
     !$acc soilhydrology_vars%wa_col )
     end if 
     !! CPU print statements !! 
     write(10,*) 'grc_wf%qflx_liq_dynbal' 
     write(10,*) grc_wf%qflx_liq_dynbal
     write(10,*) 'grc_wf%qflx_ice_dynbal' 
     write(10,*) grc_wf%qflx_ice_dynbal
     write(10,*) 'grc_ef%eflx_dynbal' 
     write(10,*) grc_ef%eflx_dynbal
     write(10,*) 'grc_ws%liq1' 
     write(10,*) grc_ws%liq1
     write(10,*) 'grc_ws%ice1' 
     write(10,*) grc_ws%ice1
     write(10,*) 'grc_ws%liq2' 
     write(10,*) grc_ws%liq2
     write(10,*) 'grc_ws%ice2' 
     write(10,*) grc_ws%ice2
     write(10,*) 'grc_es%heat2' 
     write(10,*) grc_es%heat2
     write(10,*) 'grc_es%heat1' 
     write(10,*) grc_es%heat1
     write(10,*) 'veg_ws%h2ocan' 
     write(10,*) veg_ws%h2ocan
     write(10,*) 'veg_ws%h2ocan' 
     write(10,*) veg_ws%h2ocan
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_es%t_h2osfc' 
     write(10,*) col_es%t_h2osfc
     write(10,*) 'col_es%t_soisno' 
     write(10,*) col_es%t_soisno
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'col_ws%h2osfc' 
     write(10,*) col_ws%h2osfc
     write(10,*) 'col_ws%h2osoi_liq' 
     write(10,*) col_ws%h2osoi_liq
     write(10,*) 'soilhydrology_vars%wa_col' 
     write(10,*) soilhydrology_vars%wa_col
     close(10)
end subroutine 
subroutine update_vars_dyn_cnbal_patch(gpu)
     use clm_varctl      , only : use_c13, use_c14
     use GridcellDataType, only : c14_grc_cf 
     use ColumnDataType, only : col_cf 
     use ColumnDataType, only : c14_col_cf 
     use VegetationDataType, only : c14_veg_cf 
     use ColumnDataType, only : c13_col_cf 
     use VegetationDataType, only : veg_nf 
     use GridcellDataType, only : c13_grc_cf 
     use VegetationDataType, only : veg_pf 
     use clm_instMod, only : photosyns_vars 
     use ColumnDataType, only : col_pf 
     use VegetationDataType, only : c13_veg_cf 
     use ColumnDataType, only : col_nf 
     use GridcellDataType, only : grc_cf 
     use GridcellDataType, only : grc_pf 
     use GridcellDataType, only : grc_nf 
     use VegetationDataType, only : veg_cf 
     use clm_instMod, only : cnstate_vars 
     use VegetationDataType, only : veg_ns 
     use VegetationDataType, only : veg_ps 
     use clm_instMod, only : canopystate_vars 
     use VegetationDataType, only : veg_cs 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_dyn_cnbal_patch.txt"
     else
          file='cpu_dyn_cnbal_patch.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc c14_grc_cf%dwt_prod10c_gain, & 
     !$acc c14_grc_cf%dwt_prod100c_gain, & 
     !$acc c14_grc_cf%dwt_conv_cflux, & 
     !$acc c14_grc_cf%dwt_seedc_to_leaf, & 
     !$acc c14_grc_cf%dwt_seedc_to_deadstem )
     !$acc update self(& 
     !$acc col_cf%dwt_prod10c_gain, & 
     !$acc col_cf%dwt_livecrootc_to_cwdc, & 
     !$acc col_cf%dwt_prod100c_gain, & 
     !$acc col_cf%dwt_conv_cflux, & 
     !$acc col_cf%dwt_deadcrootc_to_cwdc, & 
     !$acc col_cf%dwt_frootc_to_litr_cel_c, & 
     !$acc col_cf%dwt_frootc_to_litr_lig_c, & 
     !$acc col_cf%dwt_frootc_to_litr_met_c, & 
     !$acc col_cf%dwt_slash_cflux )
     !$acc update self(& 
     !$acc c14_col_cf%dwt_prod10c_gain, & 
     !$acc c14_col_cf%dwt_livecrootc_to_cwdc, & 
     !$acc c14_col_cf%dwt_prod100c_gain, & 
     !$acc c14_col_cf%dwt_conv_cflux, & 
     !$acc c14_col_cf%dwt_deadcrootc_to_cwdc, & 
     !$acc c14_col_cf%dwt_frootc_to_litr_cel_c, & 
     !$acc c14_col_cf%dwt_frootc_to_litr_lig_c, & 
     !$acc c14_col_cf%dwt_frootc_to_litr_met_c, & 
     !$acc c14_col_cf%dwt_slash_cflux )
     !$acc update self(& 
     !$acc c14_veg_cf%dwt_prod10c_gain, & 
     !$acc c14_veg_cf%dwt_prod100c_gain, & 
     !$acc c14_veg_cf%dwt_conv_cflux, & 
     !$acc c14_veg_cf%dwt_seedc_to_leaf, & 
     !$acc c14_veg_cf%dwt_seedc_to_deadstem )
     !$acc update self(& 
     !$acc c13_col_cf%dwt_prod10c_gain, & 
     !$acc c13_col_cf%dwt_livecrootc_to_cwdc, & 
     !$acc c13_col_cf%dwt_prod100c_gain, & 
     !$acc c13_col_cf%dwt_conv_cflux, & 
     !$acc c13_col_cf%dwt_deadcrootc_to_cwdc, & 
     !$acc c13_col_cf%dwt_frootc_to_litr_cel_c, & 
     !$acc c13_col_cf%dwt_frootc_to_litr_lig_c, & 
     !$acc c13_col_cf%dwt_frootc_to_litr_met_c, & 
     !$acc c13_col_cf%dwt_slash_cflux )
     !$acc update self(& 
     !$acc veg_nf%dwt_prod100n_gain, & 
     !$acc veg_nf%dwt_seedn_to_npool, & 
     !$acc veg_nf%dwt_prod10n_gain, & 
     !$acc veg_nf%dwt_seedn_to_deadstem, & 
     !$acc veg_nf%dwt_conv_nflux, & 
     !$acc veg_nf%dwt_seedn_to_leaf, & 
     !$acc veg_nf%plant_ndemand, & 
     !$acc veg_nf%plant_nalloc, & 
     !$acc veg_nf%avail_retransn )
     !$acc update self(& 
     !$acc c13_grc_cf%dwt_prod10c_gain, & 
     !$acc c13_grc_cf%dwt_prod100c_gain, & 
     !$acc c13_grc_cf%dwt_conv_cflux, & 
     !$acc c13_grc_cf%dwt_seedc_to_leaf, & 
     !$acc c13_grc_cf%dwt_seedc_to_deadstem )
     !$acc update self(& 
     !$acc veg_pf%dwt_conv_pflux, & 
     !$acc veg_pf%dwt_seedp_to_deadstem, & 
     !$acc veg_pf%dwt_seedp_to_ppool, & 
     !$acc veg_pf%dwt_prod10p_gain, & 
     !$acc veg_pf%dwt_prod100p_gain, & 
     !$acc veg_pf%dwt_seedp_to_leaf, & 
     !$acc veg_pf%plant_palloc, & 
     !$acc veg_pf%avail_retransp, & 
     !$acc veg_pf%plant_pdemand )
     !$acc update self(& 
     !$acc photosyns_vars%c13_psnsun_patch, & 
     !$acc photosyns_vars%psnsun_patch, & 
     !$acc photosyns_vars%rc13_psnsun_patch, & 
     !$acc photosyns_vars%alphapsnsha_patch, & 
     !$acc photosyns_vars%alphapsnsun_patch, & 
     !$acc photosyns_vars%c13_psnsha_patch, & 
     !$acc photosyns_vars%c14_psnsha_patch, & 
     !$acc photosyns_vars%rc13_canair_patch, & 
     !$acc photosyns_vars%rc13_psnsha_patch, & 
     !$acc photosyns_vars%psnsha_patch, & 
     !$acc photosyns_vars%c14_psnsun_patch )
     !$acc update self(& 
     !$acc col_pf%dwt_conv_pflux, & 
     !$acc col_pf%dwt_frootp_to_litr_met_p, & 
     !$acc col_pf%dwt_frootp_to_litr_lig_p, & 
     !$acc col_pf%dwt_livecrootp_to_cwdp, & 
     !$acc col_pf%dwt_slash_pflux, & 
     !$acc col_pf%dwt_prod10p_gain, & 
     !$acc col_pf%dwt_prod100p_gain, & 
     !$acc col_pf%dwt_frootp_to_litr_cel_p, & 
     !$acc col_pf%dwt_deadcrootp_to_cwdp )
     !$acc update self(& 
     !$acc c13_veg_cf%dwt_prod10c_gain, & 
     !$acc c13_veg_cf%dwt_prod100c_gain, & 
     !$acc c13_veg_cf%dwt_conv_cflux, & 
     !$acc c13_veg_cf%dwt_seedc_to_leaf, & 
     !$acc c13_veg_cf%dwt_seedc_to_deadstem )
     !$acc update self(& 
     !$acc col_nf%dwt_prod100n_gain, & 
     !$acc col_nf%dwt_livecrootn_to_cwdn, & 
     !$acc col_nf%dwt_prod10n_gain, & 
     !$acc col_nf%dwt_deadcrootn_to_cwdn, & 
     !$acc col_nf%dwt_frootn_to_litr_met_n, & 
     !$acc col_nf%dwt_conv_nflux, & 
     !$acc col_nf%dwt_slash_nflux, & 
     !$acc col_nf%dwt_frootn_to_litr_lig_n, & 
     !$acc col_nf%dwt_frootn_to_litr_cel_n )
     !$acc update self(& 
     !$acc grc_cf%dwt_prod10c_gain, & 
     !$acc grc_cf%dwt_prod100c_gain, & 
     !$acc grc_cf%dwt_conv_cflux, & 
     !$acc grc_cf%dwt_seedc_to_leaf, & 
     !$acc grc_cf%dwt_seedc_to_deadstem )
     !$acc update self(& 
     !$acc grc_pf%dwt_conv_pflux, & 
     !$acc grc_pf%dwt_seedp_to_deadstem, & 
     !$acc grc_pf%dwt_seedp_to_ppool, & 
     !$acc grc_pf%dwt_prod10p_gain, & 
     !$acc grc_pf%dwt_prod100p_gain, & 
     !$acc grc_pf%dwt_seedp_to_leaf )
     !$acc update self(& 
     !$acc grc_nf%dwt_prod100n_gain, & 
     !$acc grc_nf%dwt_seedn_to_npool, & 
     !$acc grc_nf%dwt_prod10n_gain, & 
     !$acc grc_nf%dwt_seedn_to_deadstem, & 
     !$acc grc_nf%dwt_conv_nflux, & 
     !$acc grc_nf%dwt_seedn_to_leaf )
     !$acc update self(& 
     !$acc veg_cf%dwt_prod10c_gain, & 
     !$acc veg_cf%dwt_prod100c_gain, & 
     !$acc veg_cf%dwt_conv_cflux, & 
     !$acc veg_cf%dwt_seedc_to_leaf, & 
     !$acc veg_cf%dwt_seedc_to_deadstem )
     !$acc update self(& 
     !$acc cnstate_vars%lfpftd_patch, & 
     !$acc cnstate_vars%p_allometry_patch, & 
     !$acc cnstate_vars%bgtr_patch, & 
     !$acc cnstate_vars%alloc_pnow_patch, & 
     !$acc cnstate_vars%bglfr_patch, & 
     !$acc cnstate_vars%downreg_patch, & 
     !$acc cnstate_vars%onset_swi_patch, & 
     !$acc cnstate_vars%annsum_potential_gpp_patch, & 
     !$acc cnstate_vars%dormant_flag_patch, & 
     !$acc cnstate_vars%bglfr_leaf_patch, & 
     !$acc cnstate_vars%tempavg_t2m_patch, & 
     !$acc cnstate_vars%tempmax_retransp_patch, & 
     !$acc cnstate_vars%onset_counter_patch, & 
     !$acc cnstate_vars%tempsum_potential_gpp_patch, & 
     !$acc cnstate_vars%rc14_atm_patch, & 
     !$acc cnstate_vars%offset_fdd_patch, & 
     !$acc cnstate_vars%days_active_patch, & 
     !$acc cnstate_vars%onset_gdd_patch, & 
     !$acc cnstate_vars%offset_counter_patch, & 
     !$acc cnstate_vars%bglfr_froot_patch, & 
     !$acc cnstate_vars%onset_gddflag_patch, & 
     !$acc cnstate_vars%n_allometry_patch, & 
     !$acc cnstate_vars%onset_fdd_patch, & 
     !$acc cnstate_vars%offset_flag_patch, & 
     !$acc cnstate_vars%annavg_t2m_patch, & 
     !$acc cnstate_vars%lgsf_patch, & 
     !$acc cnstate_vars%onset_flag_patch, & 
     !$acc cnstate_vars%tempmax_retransn_patch, & 
     !$acc cnstate_vars%c_allometry_patch, & 
     !$acc cnstate_vars%offset_swi_patch, & 
     !$acc cnstate_vars%annmax_retransp_patch, & 
     !$acc cnstate_vars%annmax_retransn_patch )
     !$acc update self(& 
     !$acc veg_ns%livestemn, & 
     !$acc veg_ns%deadstemn, & 
     !$acc veg_ns%totvegn, & 
     !$acc veg_ns%dispvegn, & 
     !$acc veg_ns%deadstemn_xfer, & 
     !$acc veg_ns%retransn, & 
     !$acc veg_ns%npool, & 
     !$acc veg_ns%frootn_xfer, & 
     !$acc veg_ns%leafn_xfer, & 
     !$acc veg_ns%livestemn_storage, & 
     !$acc veg_ns%frootn, & 
     !$acc veg_ns%leafn, & 
     !$acc veg_ns%ntrunc, & 
     !$acc veg_ns%livecrootn_storage, & 
     !$acc veg_ns%frootn_storage, & 
     !$acc veg_ns%leafn_storage, & 
     !$acc veg_ns%livestemn_xfer, & 
     !$acc veg_ns%totpftn, & 
     !$acc veg_ns%deadstemn_storage, & 
     !$acc veg_ns%deadcrootn_storage, & 
     !$acc veg_ns%deadcrootn, & 
     !$acc veg_ns%livecrootn_xfer, & 
     !$acc veg_ns%deadcrootn_xfer, & 
     !$acc veg_ns%storvegn, & 
     !$acc veg_ns%livecrootn, & 
     !$acc veg_ns%deadstemn )
     !$acc update self(& 
     !$acc veg_ps%livecrootp_xfer, & 
     !$acc veg_ps%deadcrootp_storage, & 
     !$acc veg_ps%deadcrootp_xfer, & 
     !$acc veg_ps%livecrootp, & 
     !$acc veg_ps%retransp, & 
     !$acc veg_ps%leafp_storage, & 
     !$acc veg_ps%deadstemp_storage, & 
     !$acc veg_ps%deadstemp_xfer, & 
     !$acc veg_ps%ppool, & 
     !$acc veg_ps%livecrootp_storage, & 
     !$acc veg_ps%leafp, & 
     !$acc veg_ps%livestemp_storage, & 
     !$acc veg_ps%livestemp_xfer, & 
     !$acc veg_ps%leafp_xfer, & 
     !$acc veg_ps%livestemp, & 
     !$acc veg_ps%totpftp, & 
     !$acc veg_ps%dispvegp, & 
     !$acc veg_ps%ptrunc, & 
     !$acc veg_ps%totvegp, & 
     !$acc veg_ps%deadcrootp, & 
     !$acc veg_ps%frootp_xfer, & 
     !$acc veg_ps%storvegp, & 
     !$acc veg_ps%frootp_storage, & 
     !$acc veg_ps%frootp, & 
     !$acc veg_ps%deadstemp, & 
     !$acc veg_ps%deadstemp )
     !$acc update self(& 
     !$acc canopystate_vars%laisun_patch, & 
     !$acc canopystate_vars%laisha_patch )
     !$acc update self(& 
     !$acc veg_cs%deadstemc )
     end if 
     !! CPU print statements !! 
     if ( use_c14 ) then
     write(10,*) 'c14_grc_cf%dwt_prod10c_gain' 
     write(10,*) c14_grc_cf%dwt_prod10c_gain
     write(10,*) 'c14_grc_cf%dwt_prod100c_gain' 
     write(10,*) c14_grc_cf%dwt_prod100c_gain
     write(10,*) 'c14_grc_cf%dwt_conv_cflux' 
     write(10,*) c14_grc_cf%dwt_conv_cflux
     write(10,*) 'c14_grc_cf%dwt_seedc_to_leaf' 
     write(10,*) c14_grc_cf%dwt_seedc_to_leaf
     write(10,*) 'c14_grc_cf%dwt_seedc_to_deadstem' 
     write(10,*) c14_grc_cf%dwt_seedc_to_deadstem
     endif
     write(10,*) 'col_cf%dwt_prod10c_gain' 
     write(10,*) col_cf%dwt_prod10c_gain
     write(10,*) 'col_cf%dwt_livecrootc_to_cwdc' 
     write(10,*) col_cf%dwt_livecrootc_to_cwdc
     write(10,*) 'col_cf%dwt_prod100c_gain' 
     write(10,*) col_cf%dwt_prod100c_gain
     write(10,*) 'col_cf%dwt_conv_cflux' 
     write(10,*) col_cf%dwt_conv_cflux
     write(10,*) 'col_cf%dwt_deadcrootc_to_cwdc' 
     write(10,*) col_cf%dwt_deadcrootc_to_cwdc
     write(10,*) 'col_cf%dwt_frootc_to_litr_cel_c' 
     write(10,*) col_cf%dwt_frootc_to_litr_cel_c
     write(10,*) 'col_cf%dwt_frootc_to_litr_lig_c' 
     write(10,*) col_cf%dwt_frootc_to_litr_lig_c
     write(10,*) 'col_cf%dwt_frootc_to_litr_met_c' 
     write(10,*) col_cf%dwt_frootc_to_litr_met_c
     write(10,*) 'col_cf%dwt_slash_cflux' 
     write(10,*) col_cf%dwt_slash_cflux
     if ( use_c14 ) then
     write(10,*) 'c14_col_cf%dwt_prod10c_gain' 
     write(10,*) c14_col_cf%dwt_prod10c_gain
     write(10,*) 'c14_col_cf%dwt_livecrootc_to_cwdc' 
     write(10,*) c14_col_cf%dwt_livecrootc_to_cwdc
     write(10,*) 'c14_col_cf%dwt_prod100c_gain' 
     write(10,*) c14_col_cf%dwt_prod100c_gain
     write(10,*) 'c14_col_cf%dwt_conv_cflux' 
     write(10,*) c14_col_cf%dwt_conv_cflux
     write(10,*) 'c14_col_cf%dwt_deadcrootc_to_cwdc' 
     write(10,*) c14_col_cf%dwt_deadcrootc_to_cwdc
     write(10,*) 'c14_col_cf%dwt_frootc_to_litr_cel_c' 
     write(10,*) c14_col_cf%dwt_frootc_to_litr_cel_c
     write(10,*) 'c14_col_cf%dwt_frootc_to_litr_lig_c' 
     write(10,*) c14_col_cf%dwt_frootc_to_litr_lig_c
     write(10,*) 'c14_col_cf%dwt_frootc_to_litr_met_c' 
     write(10,*) c14_col_cf%dwt_frootc_to_litr_met_c
     write(10,*) 'c14_col_cf%dwt_slash_cflux' 
     write(10,*) c14_col_cf%dwt_slash_cflux
     write(10,*) 'c14_veg_cf%dwt_prod10c_gain' 
     write(10,*) c14_veg_cf%dwt_prod10c_gain
     write(10,*) 'c14_veg_cf%dwt_prod100c_gain' 
     write(10,*) c14_veg_cf%dwt_prod100c_gain
     write(10,*) 'c14_veg_cf%dwt_conv_cflux' 
     write(10,*) c14_veg_cf%dwt_conv_cflux
     write(10,*) 'c14_veg_cf%dwt_seedc_to_leaf' 
     write(10,*) c14_veg_cf%dwt_seedc_to_leaf
     write(10,*) 'c14_veg_cf%dwt_seedc_to_deadstem' 
     write(10,*) c14_veg_cf%dwt_seedc_to_deadstem
     endif
     if ( use_c14 ) then
     write(10,*) 'c13_col_cf%dwt_prod10c_gain' 
     write(10,*) c13_col_cf%dwt_prod10c_gain
     write(10,*) 'c13_col_cf%dwt_livecrootc_to_cwdc' 
     write(10,*) c13_col_cf%dwt_livecrootc_to_cwdc
     write(10,*) 'c13_col_cf%dwt_prod100c_gain' 
     write(10,*) c13_col_cf%dwt_prod100c_gain
     write(10,*) 'c13_col_cf%dwt_conv_cflux' 
     write(10,*) c13_col_cf%dwt_conv_cflux
     write(10,*) 'c13_col_cf%dwt_deadcrootc_to_cwdc' 
     write(10,*) c13_col_cf%dwt_deadcrootc_to_cwdc
     write(10,*) 'c13_col_cf%dwt_frootc_to_litr_cel_c' 
     write(10,*) c13_col_cf%dwt_frootc_to_litr_cel_c
     write(10,*) 'c13_col_cf%dwt_frootc_to_litr_lig_c' 
     write(10,*) c13_col_cf%dwt_frootc_to_litr_lig_c
     write(10,*) 'c13_col_cf%dwt_frootc_to_litr_met_c' 
     write(10,*) c13_col_cf%dwt_frootc_to_litr_met_c
     write(10,*) 'c13_col_cf%dwt_slash_cflux' 
     write(10,*) c13_col_cf%dwt_slash_cflux
     endif
     write(10,*) 'veg_nf%dwt_prod100n_gain' 
     write(10,*) veg_nf%dwt_prod100n_gain
     write(10,*) 'veg_nf%dwt_seedn_to_npool' 
     write(10,*) veg_nf%dwt_seedn_to_npool
     write(10,*) 'veg_nf%dwt_prod10n_gain' 
     write(10,*) veg_nf%dwt_prod10n_gain
     write(10,*) 'veg_nf%dwt_seedn_to_deadstem' 
     write(10,*) veg_nf%dwt_seedn_to_deadstem
     write(10,*) 'veg_nf%dwt_conv_nflux' 
     write(10,*) veg_nf%dwt_conv_nflux
     write(10,*) 'veg_nf%dwt_seedn_to_leaf' 
     write(10,*) veg_nf%dwt_seedn_to_leaf
     write(10,*) 'veg_nf%plant_ndemand' 
     write(10,*) veg_nf%plant_ndemand
     write(10,*) 'veg_nf%plant_nalloc' 
     write(10,*) veg_nf%plant_nalloc
     write(10,*) 'veg_nf%avail_retransn' 
     write(10,*) veg_nf%avail_retransn
     if ( use_c13 ) then
     write(10,*) 'c13_grc_cf%dwt_prod10c_gain' 
     write(10,*) c13_grc_cf%dwt_prod10c_gain
     write(10,*) 'c13_grc_cf%dwt_prod100c_gain' 
     write(10,*) c13_grc_cf%dwt_prod100c_gain
     write(10,*) 'c13_grc_cf%dwt_conv_cflux' 
     write(10,*) c13_grc_cf%dwt_conv_cflux
     write(10,*) 'c13_grc_cf%dwt_seedc_to_leaf' 
     write(10,*) c13_grc_cf%dwt_seedc_to_leaf
     write(10,*) 'c13_grc_cf%dwt_seedc_to_deadstem' 
     write(10,*) c13_grc_cf%dwt_seedc_to_deadstem
     endif
     write(10,*) 'veg_pf%dwt_conv_pflux' 
     write(10,*) veg_pf%dwt_conv_pflux
     write(10,*) 'veg_pf%dwt_seedp_to_deadstem' 
     write(10,*) veg_pf%dwt_seedp_to_deadstem
     write(10,*) 'veg_pf%dwt_seedp_to_ppool' 
     write(10,*) veg_pf%dwt_seedp_to_ppool
     write(10,*) 'veg_pf%dwt_prod10p_gain' 
     write(10,*) veg_pf%dwt_prod10p_gain
     write(10,*) 'veg_pf%dwt_prod100p_gain' 
     write(10,*) veg_pf%dwt_prod100p_gain
     write(10,*) 'veg_pf%dwt_seedp_to_leaf' 
     write(10,*) veg_pf%dwt_seedp_to_leaf
     write(10,*) 'veg_pf%plant_palloc' 
     write(10,*) veg_pf%plant_palloc
     write(10,*) 'veg_pf%avail_retransp' 
     write(10,*) veg_pf%avail_retransp
     write(10,*) 'veg_pf%plant_pdemand' 
     write(10,*) veg_pf%plant_pdemand
     write(10,*) 'photosyns_vars%c13_psnsun_patch' 
     write(10,*) photosyns_vars%c13_psnsun_patch
     write(10,*) 'photosyns_vars%psnsun_patch' 
     write(10,*) photosyns_vars%psnsun_patch
     write(10,*) 'photosyns_vars%rc13_psnsun_patch' 
     write(10,*) photosyns_vars%rc13_psnsun_patch
     write(10,*) 'photosyns_vars%alphapsnsha_patch' 
     write(10,*) photosyns_vars%alphapsnsha_patch
     write(10,*) 'photosyns_vars%alphapsnsun_patch' 
     write(10,*) photosyns_vars%alphapsnsun_patch
     write(10,*) 'photosyns_vars%c13_psnsha_patch' 
     write(10,*) photosyns_vars%c13_psnsha_patch
     write(10,*) 'photosyns_vars%c14_psnsha_patch' 
     write(10,*) photosyns_vars%c14_psnsha_patch
     write(10,*) 'photosyns_vars%rc13_canair_patch' 
     write(10,*) photosyns_vars%rc13_canair_patch
     write(10,*) 'photosyns_vars%rc13_psnsha_patch' 
     write(10,*) photosyns_vars%rc13_psnsha_patch
     write(10,*) 'photosyns_vars%psnsha_patch' 
     write(10,*) photosyns_vars%psnsha_patch
     write(10,*) 'photosyns_vars%c14_psnsun_patch' 
     write(10,*) photosyns_vars%c14_psnsun_patch
     write(10,*) 'col_pf%dwt_conv_pflux' 
     write(10,*) col_pf%dwt_conv_pflux
     write(10,*) 'col_pf%dwt_frootp_to_litr_met_p' 
     write(10,*) col_pf%dwt_frootp_to_litr_met_p
     write(10,*) 'col_pf%dwt_frootp_to_litr_lig_p' 
     write(10,*) col_pf%dwt_frootp_to_litr_lig_p
     write(10,*) 'col_pf%dwt_livecrootp_to_cwdp' 
     write(10,*) col_pf%dwt_livecrootp_to_cwdp
     write(10,*) 'col_pf%dwt_slash_pflux' 
     write(10,*) col_pf%dwt_slash_pflux
     write(10,*) 'col_pf%dwt_prod10p_gain' 
     write(10,*) col_pf%dwt_prod10p_gain
     write(10,*) 'col_pf%dwt_prod100p_gain' 
     write(10,*) col_pf%dwt_prod100p_gain
     write(10,*) 'col_pf%dwt_frootp_to_litr_cel_p' 
     write(10,*) col_pf%dwt_frootp_to_litr_cel_p
     write(10,*) 'col_pf%dwt_deadcrootp_to_cwdp' 
     write(10,*) col_pf%dwt_deadcrootp_to_cwdp
     if ( use_c13 ) then
     write(10,*) 'c13_veg_cf%dwt_prod10c_gain' 
     write(10,*) c13_veg_cf%dwt_prod10c_gain
     write(10,*) 'c13_veg_cf%dwt_prod100c_gain' 
     write(10,*) c13_veg_cf%dwt_prod100c_gain
     write(10,*) 'c13_veg_cf%dwt_conv_cflux' 
     write(10,*) c13_veg_cf%dwt_conv_cflux
     write(10,*) 'c13_veg_cf%dwt_seedc_to_leaf' 
     write(10,*) c13_veg_cf%dwt_seedc_to_leaf
     write(10,*) 'c13_veg_cf%dwt_seedc_to_deadstem' 
     write(10,*) c13_veg_cf%dwt_seedc_to_deadstem
     endif
     write(10,*) 'col_nf%dwt_prod100n_gain' 
     write(10,*) col_nf%dwt_prod100n_gain
     write(10,*) 'col_nf%dwt_livecrootn_to_cwdn' 
     write(10,*) col_nf%dwt_livecrootn_to_cwdn
     write(10,*) 'col_nf%dwt_prod10n_gain' 
     write(10,*) col_nf%dwt_prod10n_gain
     write(10,*) 'col_nf%dwt_deadcrootn_to_cwdn' 
     write(10,*) col_nf%dwt_deadcrootn_to_cwdn
     write(10,*) 'col_nf%dwt_frootn_to_litr_met_n' 
     write(10,*) col_nf%dwt_frootn_to_litr_met_n
     write(10,*) 'col_nf%dwt_conv_nflux' 
     write(10,*) col_nf%dwt_conv_nflux
     write(10,*) 'col_nf%dwt_slash_nflux' 
     write(10,*) col_nf%dwt_slash_nflux
     write(10,*) 'col_nf%dwt_frootn_to_litr_lig_n' 
     write(10,*) col_nf%dwt_frootn_to_litr_lig_n
     write(10,*) 'col_nf%dwt_frootn_to_litr_cel_n' 
     write(10,*) col_nf%dwt_frootn_to_litr_cel_n
     write(10,*) 'grc_cf%dwt_prod10c_gain' 
     write(10,*) grc_cf%dwt_prod10c_gain
     write(10,*) 'grc_cf%dwt_prod100c_gain' 
     write(10,*) grc_cf%dwt_prod100c_gain
     write(10,*) 'grc_cf%dwt_conv_cflux' 
     write(10,*) grc_cf%dwt_conv_cflux
     write(10,*) 'grc_cf%dwt_seedc_to_leaf' 
     write(10,*) grc_cf%dwt_seedc_to_leaf
     write(10,*) 'grc_cf%dwt_seedc_to_deadstem' 
     write(10,*) grc_cf%dwt_seedc_to_deadstem
     write(10,*) 'grc_pf%dwt_conv_pflux' 
     write(10,*) grc_pf%dwt_conv_pflux
     write(10,*) 'grc_pf%dwt_seedp_to_deadstem' 
     write(10,*) grc_pf%dwt_seedp_to_deadstem
     write(10,*) 'grc_pf%dwt_seedp_to_ppool' 
     write(10,*) grc_pf%dwt_seedp_to_ppool
     write(10,*) 'grc_pf%dwt_prod10p_gain' 
     write(10,*) grc_pf%dwt_prod10p_gain
     write(10,*) 'grc_pf%dwt_prod100p_gain' 
     write(10,*) grc_pf%dwt_prod100p_gain
     write(10,*) 'grc_pf%dwt_seedp_to_leaf' 
     write(10,*) grc_pf%dwt_seedp_to_leaf
     write(10,*) 'grc_nf%dwt_prod100n_gain' 
     write(10,*) grc_nf%dwt_prod100n_gain
     write(10,*) 'grc_nf%dwt_seedn_to_npool' 
     write(10,*) grc_nf%dwt_seedn_to_npool
     write(10,*) 'grc_nf%dwt_prod10n_gain' 
     write(10,*) grc_nf%dwt_prod10n_gain
     write(10,*) 'grc_nf%dwt_seedn_to_deadstem' 
     write(10,*) grc_nf%dwt_seedn_to_deadstem
     write(10,*) 'grc_nf%dwt_conv_nflux' 
     write(10,*) grc_nf%dwt_conv_nflux
     write(10,*) 'grc_nf%dwt_seedn_to_leaf' 
     write(10,*) grc_nf%dwt_seedn_to_leaf
     write(10,*) 'veg_cf%dwt_prod10c_gain' 
     write(10,*) veg_cf%dwt_prod10c_gain
     write(10,*) 'veg_cf%dwt_prod100c_gain' 
     write(10,*) veg_cf%dwt_prod100c_gain
     write(10,*) 'veg_cf%dwt_conv_cflux' 
     write(10,*) veg_cf%dwt_conv_cflux
     write(10,*) 'veg_cf%dwt_seedc_to_leaf' 
     write(10,*) veg_cf%dwt_seedc_to_leaf
     write(10,*) 'veg_cf%dwt_seedc_to_deadstem' 
     write(10,*) veg_cf%dwt_seedc_to_deadstem
     write(10,*) 'cnstate_vars%lfpftd_patch' 
     write(10,*) cnstate_vars%lfpftd_patch
     write(10,*) 'cnstate_vars%p_allometry_patch' 
     write(10,*) cnstate_vars%p_allometry_patch
     write(10,*) 'cnstate_vars%bgtr_patch' 
     write(10,*) cnstate_vars%bgtr_patch
     write(10,*) 'cnstate_vars%alloc_pnow_patch' 
     write(10,*) cnstate_vars%alloc_pnow_patch
     write(10,*) 'cnstate_vars%bglfr_patch' 
     write(10,*) cnstate_vars%bglfr_patch
     write(10,*) 'cnstate_vars%downreg_patch' 
     write(10,*) cnstate_vars%downreg_patch
     write(10,*) 'cnstate_vars%onset_swi_patch' 
     write(10,*) cnstate_vars%onset_swi_patch
     write(10,*) 'cnstate_vars%annsum_potential_gpp_patch' 
     write(10,*) cnstate_vars%annsum_potential_gpp_patch
     write(10,*) 'cnstate_vars%dormant_flag_patch' 
     write(10,*) cnstate_vars%dormant_flag_patch
     write(10,*) 'cnstate_vars%bglfr_leaf_patch' 
     write(10,*) cnstate_vars%bglfr_leaf_patch
     write(10,*) 'cnstate_vars%tempavg_t2m_patch' 
     write(10,*) cnstate_vars%tempavg_t2m_patch
     write(10,*) 'cnstate_vars%tempmax_retransp_patch' 
     write(10,*) cnstate_vars%tempmax_retransp_patch
     write(10,*) 'cnstate_vars%onset_counter_patch' 
     write(10,*) cnstate_vars%onset_counter_patch
     write(10,*) 'cnstate_vars%tempsum_potential_gpp_patch' 
     write(10,*) cnstate_vars%tempsum_potential_gpp_patch
     write(10,*) 'cnstate_vars%rc14_atm_patch' 
     write(10,*) cnstate_vars%rc14_atm_patch
     write(10,*) 'cnstate_vars%offset_fdd_patch' 
     write(10,*) cnstate_vars%offset_fdd_patch
     write(10,*) 'cnstate_vars%days_active_patch' 
     write(10,*) cnstate_vars%days_active_patch
     write(10,*) 'cnstate_vars%onset_gdd_patch' 
     write(10,*) cnstate_vars%onset_gdd_patch
     write(10,*) 'cnstate_vars%offset_counter_patch' 
     write(10,*) cnstate_vars%offset_counter_patch
     write(10,*) 'cnstate_vars%bglfr_froot_patch' 
     write(10,*) cnstate_vars%bglfr_froot_patch
     write(10,*) 'cnstate_vars%onset_gddflag_patch' 
     write(10,*) cnstate_vars%onset_gddflag_patch
     write(10,*) 'cnstate_vars%n_allometry_patch' 
     write(10,*) cnstate_vars%n_allometry_patch
     write(10,*) 'cnstate_vars%onset_fdd_patch' 
     write(10,*) cnstate_vars%onset_fdd_patch
     write(10,*) 'cnstate_vars%offset_flag_patch' 
     write(10,*) cnstate_vars%offset_flag_patch
     write(10,*) 'cnstate_vars%annavg_t2m_patch' 
     write(10,*) cnstate_vars%annavg_t2m_patch
     write(10,*) 'cnstate_vars%lgsf_patch' 
     write(10,*) cnstate_vars%lgsf_patch
     write(10,*) 'cnstate_vars%onset_flag_patch' 
     write(10,*) cnstate_vars%onset_flag_patch
     write(10,*) 'cnstate_vars%tempmax_retransn_patch' 
     write(10,*) cnstate_vars%tempmax_retransn_patch
     write(10,*) 'cnstate_vars%c_allometry_patch' 
     write(10,*) cnstate_vars%c_allometry_patch
     write(10,*) 'cnstate_vars%offset_swi_patch' 
     write(10,*) cnstate_vars%offset_swi_patch
     write(10,*) 'cnstate_vars%annmax_retransp_patch' 
     write(10,*) cnstate_vars%annmax_retransp_patch
     write(10,*) 'cnstate_vars%annmax_retransn_patch' 
     write(10,*) cnstate_vars%annmax_retransn_patch
     write(10,*) 'veg_ns%livestemn' 
     write(10,*) veg_ns%livestemn
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ns%totvegn' 
     write(10,*) veg_ns%totvegn
     write(10,*) 'veg_ns%dispvegn' 
     write(10,*) veg_ns%dispvegn
     write(10,*) 'veg_ns%deadstemn_xfer' 
     write(10,*) veg_ns%deadstemn_xfer
     write(10,*) 'veg_ns%retransn' 
     write(10,*) veg_ns%retransn
     write(10,*) 'veg_ns%npool' 
     write(10,*) veg_ns%npool
     write(10,*) 'veg_ns%frootn_xfer' 
     write(10,*) veg_ns%frootn_xfer
     write(10,*) 'veg_ns%leafn_xfer' 
     write(10,*) veg_ns%leafn_xfer
     write(10,*) 'veg_ns%livestemn_storage' 
     write(10,*) veg_ns%livestemn_storage
     write(10,*) 'veg_ns%frootn' 
     write(10,*) veg_ns%frootn
     write(10,*) 'veg_ns%leafn' 
     write(10,*) veg_ns%leafn
     write(10,*) 'veg_ns%ntrunc' 
     write(10,*) veg_ns%ntrunc
     write(10,*) 'veg_ns%livecrootn_storage' 
     write(10,*) veg_ns%livecrootn_storage
     write(10,*) 'veg_ns%frootn_storage' 
     write(10,*) veg_ns%frootn_storage
     write(10,*) 'veg_ns%leafn_storage' 
     write(10,*) veg_ns%leafn_storage
     write(10,*) 'veg_ns%livestemn_xfer' 
     write(10,*) veg_ns%livestemn_xfer
     write(10,*) 'veg_ns%totpftn' 
     write(10,*) veg_ns%totpftn
     write(10,*) 'veg_ns%deadstemn_storage' 
     write(10,*) veg_ns%deadstemn_storage
     write(10,*) 'veg_ns%deadcrootn_storage' 
     write(10,*) veg_ns%deadcrootn_storage
     write(10,*) 'veg_ns%deadcrootn' 
     write(10,*) veg_ns%deadcrootn
     write(10,*) 'veg_ns%livecrootn_xfer' 
     write(10,*) veg_ns%livecrootn_xfer
     write(10,*) 'veg_ns%deadcrootn_xfer' 
     write(10,*) veg_ns%deadcrootn_xfer
     write(10,*) 'veg_ns%storvegn' 
     write(10,*) veg_ns%storvegn
     write(10,*) 'veg_ns%livecrootn' 
     write(10,*) veg_ns%livecrootn
     write(10,*) 'veg_ns%deadstemn' 
     write(10,*) veg_ns%deadstemn
     write(10,*) 'veg_ps%livecrootp_xfer' 
     write(10,*) veg_ps%livecrootp_xfer
     write(10,*) 'veg_ps%deadcrootp_storage' 
     write(10,*) veg_ps%deadcrootp_storage
     write(10,*) 'veg_ps%deadcrootp_xfer' 
     write(10,*) veg_ps%deadcrootp_xfer
     write(10,*) 'veg_ps%livecrootp' 
     write(10,*) veg_ps%livecrootp
     write(10,*) 'veg_ps%retransp' 
     write(10,*) veg_ps%retransp
     write(10,*) 'veg_ps%leafp_storage' 
     write(10,*) veg_ps%leafp_storage
     write(10,*) 'veg_ps%deadstemp_storage' 
     write(10,*) veg_ps%deadstemp_storage
     write(10,*) 'veg_ps%deadstemp_xfer' 
     write(10,*) veg_ps%deadstemp_xfer
     write(10,*) 'veg_ps%ppool' 
     write(10,*) veg_ps%ppool
     write(10,*) 'veg_ps%livecrootp_storage' 
     write(10,*) veg_ps%livecrootp_storage
     write(10,*) 'veg_ps%leafp' 
     write(10,*) veg_ps%leafp
     write(10,*) 'veg_ps%livestemp_storage' 
     write(10,*) veg_ps%livestemp_storage
     write(10,*) 'veg_ps%livestemp_xfer' 
     write(10,*) veg_ps%livestemp_xfer
     write(10,*) 'veg_ps%leafp_xfer' 
     write(10,*) veg_ps%leafp_xfer
     write(10,*) 'veg_ps%livestemp' 
     write(10,*) veg_ps%livestemp
     write(10,*) 'veg_ps%totpftp' 
     write(10,*) veg_ps%totpftp
     write(10,*) 'veg_ps%dispvegp' 
     write(10,*) veg_ps%dispvegp
     write(10,*) 'veg_ps%ptrunc' 
     write(10,*) veg_ps%ptrunc
     write(10,*) 'veg_ps%totvegp' 
     write(10,*) veg_ps%totvegp
     write(10,*) 'veg_ps%deadcrootp' 
     write(10,*) veg_ps%deadcrootp
     write(10,*) 'veg_ps%frootp_xfer' 
     write(10,*) veg_ps%frootp_xfer
     write(10,*) 'veg_ps%storvegp' 
     write(10,*) veg_ps%storvegp
     write(10,*) 'veg_ps%frootp_storage' 
     write(10,*) veg_ps%frootp_storage
     write(10,*) 'veg_ps%frootp' 
     write(10,*) veg_ps%frootp
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'veg_ps%deadstemp' 
     write(10,*) veg_ps%deadstemp
     write(10,*) 'canopystate_vars%laisun_patch' 
     write(10,*) canopystate_vars%laisun_patch
     write(10,*) 'canopystate_vars%laisha_patch' 
     write(10,*) canopystate_vars%laisha_patch
     write(10,*) 'veg_cs%deadstemc' 
     write(10,*) veg_cs%deadstemc
     close(10)
end subroutine 
subroutine update_vars_dyn_cnbal_column(gpu)
     use ColumnDataType, only : col_cs 
     use ColumnDataType, only : col_ns 
     use ColumnDataType, only : col_ps 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_dyn_cnbal_column.txt"
     else
          file='cpu_dyn_cnbal_column.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_cs%ctrunc_vr, & 
     !$acc col_cs%decomp_cpools_vr, & 
     !$acc col_cs%dyn_cbal_adjustments )
     !$acc update self(& 
     !$acc col_ns%ntrunc_vr, & 
     !$acc col_ns%dyn_nbal_adjustments, & 
     !$acc col_ns%decomp_npools_vr, & 
     !$acc col_ns%sminn_vr )
     !$acc update self(& 
     !$acc col_ps%labilep_vr, & 
     !$acc col_ps%decomp_ppools_vr, & 
     !$acc col_ps%ptrunc_vr, & 
     !$acc col_ps%secondp_vr, & 
     !$acc col_ps%dyn_pbal_adjustments, & 
     !$acc col_ps%solutionp_vr )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_cs%ctrunc_vr' 
     write(10,*) col_cs%ctrunc_vr
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     write(10,*) 'col_cs%dyn_cbal_adjustments' 
     write(10,*) col_cs%dyn_cbal_adjustments
     write(10,*) 'col_ns%ntrunc_vr' 
     write(10,*) col_ns%ntrunc_vr
     write(10,*) 'col_ns%dyn_nbal_adjustments' 
     write(10,*) col_ns%dyn_nbal_adjustments
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'col_ns%sminn_vr' 
     write(10,*) col_ns%sminn_vr
     write(10,*) 'col_ps%labilep_vr' 
     write(10,*) col_ps%labilep_vr
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     write(10,*) 'col_ps%ptrunc_vr' 
     write(10,*) col_ps%ptrunc_vr
     write(10,*) 'col_ps%secondp_vr' 
     write(10,*) col_ps%secondp_vr
     write(10,*) 'col_ps%dyn_pbal_adjustments' 
     write(10,*) col_ps%dyn_pbal_adjustments
     write(10,*) 'col_ps%solutionp_vr' 
     write(10,*) col_ps%solutionp_vr
     close(10)
end subroutine 
subroutine update_vars_CarbonStateUpdateDynPatch(gpu)
     use GridcellDataType, only : grc_cs 
     use ColumnDataType, only : col_cs 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_CarbonStateUpdateDynPatch.txt"
     else
          file='cpu_CarbonStateUpdateDynPatch.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc grc_cs%seedc )
     !$acc update self(& 
     !$acc col_cs%decomp_cpools_vr )
     end if 
     !! CPU print statements !! 
     write(10,*) 'grc_cs%seedc' 
     write(10,*) grc_cs%seedc
     write(10,*) 'col_cs%decomp_cpools_vr' 
     write(10,*) col_cs%decomp_cpools_vr
     close(10)
end subroutine 
subroutine update_vars_NitrogenStateUpdateDynPatch(gpu)
     use ColumnDataType, only : col_ns 
     use GridcellDataType, only : grc_ns 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_NitrogenStateUpdateDynPatch.txt"
     else
          file='cpu_NitrogenStateUpdateDynPatch.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc col_ns%decomp_npools_vr )
     !$acc update self(& 
     !$acc grc_ns%seedn )
     end if 
     !! CPU print statements !! 
     write(10,*) 'col_ns%decomp_npools_vr' 
     write(10,*) col_ns%decomp_npools_vr
     write(10,*) 'grc_ns%seedn' 
     write(10,*) grc_ns%seedn
     close(10)
end subroutine 
subroutine update_vars_PhosphorusStateUpdateDynPatch(gpu)
     use GridcellDataType, only : grc_ps 
     use ColumnDataType, only : col_ps 
     implicit none 
     integer, intent(in) :: gpu
     character(len=256) :: file
     if(gpu>0) then
          file="gpu_PhosphorusStateUpdateDynPatch.txt"
     else
          file='cpu_PhosphorusStateUpdateDynPatch.txt'
     end if
     open(UNIT=10, STATUS='REPLACE', FILE=file)
     if(gpu>0) then
     !$acc update self(& 
     !$acc grc_ps%seedp )
     !$acc update self(& 
     !$acc col_ps%decomp_ppools_vr )
     end if 
     !! CPU print statements !! 
     write(10,*) 'grc_ps%seedp' 
     write(10,*) grc_ps%seedp
     write(10,*) 'col_ps%decomp_ppools_vr' 
     write(10,*) col_ps%decomp_ppools_vr
     close(10)
end subroutine 
end module verificationMod
