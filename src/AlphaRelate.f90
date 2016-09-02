
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

  if (MakeA) then
    call MakeAMatrix
  end if

  if (MakeInvA) then
    call MakeInvAMatrix
  end if

  if (MakeG .or. MakeInvG) then
    call MakeGAndInvGMatrix
  end if

  if (MakeH .or. MakeInvH) then
    call MakeHAndInvHMatrix
  end if

  call cpu_time(finish)
  print *," "
  print '("  Time duration of AlphaRelate = ",f20.4," seconds.")',finish-start
  print *," "

end program

!###############################################################################
