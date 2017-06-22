module time_data

implicit none

integer(kind=4) :: iclock,iclock_init
integer(kind=4) :: iclock_gather1,iclock_gather2,iclock_gather3
integer(kind=4) :: iclock_solve1,iclock_solve2,iclock_solve3
integer(kind=4) :: iclock_scatter1,iclock_scatter2,iclock_scatter3
integer(kind=4) :: iclock_i1,iclock_i2,iclock_i3,iclock_i4
real   (kind=8) :: dtime,dtime_init
real   (kind=8) :: dtime_gather1,dtime_gather2,dtime_gather3
real   (kind=8) :: dtime_scatter1,dtime_scatter2,dtime_scatter3
real   (kind=8) :: dtime_solve1,dtime_solve2,dtime_solve3
real   (kind=8) :: dtime_i1,dtime_i2,dtime_i3,dtime_i4

end module 
