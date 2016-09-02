
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

  if (MakeAInv) then
    call MakeInvAMatrix
  end if

  if (MakeG .or. MakeGInv) then
    call MakeGAndGInvMatrix
  end if

  if (MakeHInv .or. MakeH) then
    call MakeHAndHInvMatrix
  end if

  call cpu_time(finish)
  print *," "
  print '("  Time duration of AlphaRelate = ",f20.4," seconds.")',finish-start
  print *," "

end program

!###############################################################################
