
!###############################################################################

program AlphaRelate
  use ISO_Fortran_env
  use AlphaRelateModule

  implicit none

  real(real32) :: start, finish
  type(AlphaRelateSpec) :: Spec
  type(AlphaRelateData) :: Data

  call cpu_time(start)
  call AlphaRelateTitle

  Spec = AlphaRelateSpec(SpecFile="AlphaRelateSpec.txt")
  Data = AlphaRelateData(Spec=Spec)

  if (Spec%PedInbreeding) then
    call Data%CalcPedInbreeding()
    call Data%WritePedInbreeding(File="PedigreeInbreeding.txt")
  end if

  !if (Spec%MakeA) then
  !  call MakeAMatrix(Spec=Spec, Data=Data, ???)
  !end if

  !if (Spec%MakeInvA) then
  !  call MakeInvAMatrix
  !end if

  !if (Spec%MakeG .or. Spec%MakeInvG) then
  !  call MakeGAndInvGMatrix(Spec=Spec, x=Data, ???)
  !end if

  !if (MakeH .or. MakeInvH) then
  !  call MakeHAndInvHMatrix
  !end if

  call cpu_time(finish)
  print *," "
  print '("  Time duration of AlphaRelate = ",f20.4," seconds.")',finish-start
  print *," "

end program

!###############################################################################
