module timeinfoMod

  use shr_kind_mod       , only : r8 => shr_kind_r8
  use clm_time_manager

  !variables needed for time steps
  implicit none
  real(r8) :: dtime_mod
  real(r8) :: dayspyr_mod
  integer  :: year_mod
  integer  :: mon_mod
  integer  :: day_mod
  integer  :: secs_mod
  integer  :: nstep_mod
  integer  :: jday_mod


  !$acc declare create( dtime_mod   )
  !$acc declare create(dayspyr_mod )
  !$acc declare create(year_mod     )
  !$acc declare create(mon_mod      )
  !$acc declare create(day_mod      )
  !$acc declare create(secs_mod     )
  !$acc declare create(nstep_mod    )
  !$acc declare create(jday_mod     )
contains

  subroutine set_time_vars()
    implicit none

    nstep_mod = get_nstep()
    dtime_mod = real( get_step_size(), r8 )
    dayspyr_mod =  get_days_per_year()

    call get_curr_date(year_mod, mon_mod, day_mod, secs_mod)
    jday_mod    = get_curr_calday()

    

  end subroutine



end MODULE
