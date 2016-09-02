
!###############################################################################

program AlphaRelate

  use AlphaRelateMod

  implicit none

  include "mkl_vml.f90"

  real(real32) :: start, finish

  call cpu_time(start)
  call AlphaRelateTitles

  call ReadParam
  call ReadData

  if (AMake) then
    call MakeAMatrix
  end if

  if (AInvMake) then
    call MakeInvAMatrix
  end if

  if (GMake .or. GInvMake) then
    call MakeG
  end if

  if (HInvMake .or. HMake) then
    call MakeHAndHInv
  end if

  call cpu_time(finish)
  print *," "
  print '("  Time duration of AlphaRelate = ",f20.4," seconds.")',finish-start
  print *," "

end program

!###############################################################################
