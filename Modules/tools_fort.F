!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROMS ROUTINES
!!
!! copied from actual ROMS scripts
!!
!! compile with:
!! "cpp R_tools_fort.F R_tools_fort.f"
!! "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m R_tools_fort R_tools_fort.f" for python use
!!
!! print R_tools_fort.rho_eos.__doc__
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!# include "tools_fort_routines/zlevs.F"
# include "tools_fort_routines/sigma_to_z_intr.F"
# include "tools_fort_routines/sigma_to_z_intr_bot.F"
# include "tools_fort_routines/sigma_to_z_intr_bot_2d.F"
# include "tools_fort_routines/sigma_to_z_intr_bounded.F"
# include "tools_fort_routines/sigma_to_z_intr_sfc.F"
# include "tools_fort_routines/single_connect.F"
!!! diverse !!!
# include "tools_fort_routines/roms_read_write.F"
# include "tools_fort_routines/srtopo.F"
!!! compute hraw !!!
# include "tools_fort_routines/compute_hraw.F"

!!! smooth !!!
# include "tools_fort_routines/topo_smooth_subs.F"
# include "tools_fort_routines/tools_topo.F"
# include "tools_fort_routines/smooth.F"

!!! lsmooth !!!
# define LOG_SMOOTHING
# include "tools_fort_routines/smooth.F"
# undef LOG_SMOOTHING

!!! legacy_lsmooth !!!
# define LEGACY_LOG_SMOOTH
# include "tools_fort_routines/smooth.F"
# undef LEGACY_LOG_SMOOTH

!!! lsmooth2 !!!
# define LOG_SMOOTH_2
# include "tools_fort_routines/smooth.F"
# undef LOG_SMOOTH_2

!!! lsmooth1 !!!
# define LOG_SMOOTH_1
# include "tools_fort_routines/smooth.F"
# undef LOG_SMOOTH_1

!!! non_decreasing_rx_cond !!!
# define NON_DECREASING_RX_COND
# include "tools_fort_routines/cond_rx0_topo.F"
# include "tools_fort_routines/smooth.F"
# undef NON_DECREASING_RX_COND

# include "tools_fort_routines/r2r_match_topo.F"

!!! mask !!!
# include "tools_fort_routines/gshhs_to_roms_mask.F"

# include "tools_fort_routines/r2r_init.F"
# include "tools_fort_routines/r2r_bry.F"
# include "tools_fort_routines/etch_into_land.F"
# include "tools_fort_routines/set_depth.F"
# include "tools_fort_routines/def_roms_file.F"
# include "tools_fort_routines/def_bry_file.F"
# include "tools_fort_routines/r2r_subs.F"
# include "tools_fort_routines/r2r_vert_interp.F"
# include "tools_fort_routines/r2r_rotate.F"
# include "tools_fort_routines/r2r_bry_rotate.F"

!!! tools used for smoothing !!!
# include "tools_fort_routines/smooth_utils.F"
# include "tools_fort_routines/tiling.F"
# include "tools_fort_routines/read_string.F"
# include "tools_fort_routines/lenstr.F"
# include "tools_fort_routines/mrg_zone_subs.F"
# include "tools_fort_routines/r2r_interp_init.F"
# include "tools_fort_routines/r2r_bry_interp.F"
!!!!!!!!!!

