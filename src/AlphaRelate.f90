
!###############################################################################

program AlphaRelate

  use AlphaRelateModule

  implicit none

  real(real32) :: start, finish
  type(AlphaRelateSpec) :: Spec
  type(AlphaRelateData) :: Data

  integer(int32) :: i
  integer(int32) :: PedInbreedingUnit

  call cpu_time(start)
  call AlphaRelateTitles

  Spec = AlphaRelateSpec(SpecFile="AlphaRelateSpec.txt")
  Data = AlphaRelateData(Spec=Spec)

  if (Spec%MakePedInb) then
    allocate(Data%PedInbreeding(0:Data%nIndPed))
    Data%PedInbreeding = PedInbreeding(RecPed=Data%RecPed%id, n=Data%nIndPed)
    open(newunit=PedInbreedingUnit, file="PedigreeInbreeding.txt", status="unknown")
    do i = 1, Data%nIndPed
      write(PedInbreedingUnit, "(a20,f20.16)") Data%RecPed%OriginalId(i), Data%PedInbreeding(i)
    end do
    close(PedInbreedingUnit)
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
