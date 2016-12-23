
!###############################################################################

program AlphaRelate
  use ISO_Fortran_env
  use AlphaHouseMod, only : PrintElapsedTime
  use AlphaRelateModule

  implicit none

  real(real32) :: StartTime, EndTime
  type(AlphaRelateSpec) :: Spec
  type(AlphaRelateData) :: Data

  call cpu_time(StartTime)
  call AlphaRelateTitle

  ! TODO: get spec file from command line, if not use default
  Spec = AlphaRelateSpec(SpecFile="AlphaRelateSpec.txt")
  Data = AlphaRelateData(Spec=Spec)

  if (Spec%PedInbreeding) then
    call Data%CalcPedInbreeding()
    call Data%WritePedInbreeding(File="PedigreeInbreeding.txt", Spec=Spec)
  end if

  if (Spec%PedNrm) then
    call Data%CalcPedNrm(Spec=Spec)
    call Data%WritePedNrm(File="PedigreeNrm.txt", Spec=Spec)
  end if

  !if (Spec%PedNrmInv) then
  !  call Data%CalcPedNrmInv()
  !  call Data%WritePedNrmInv(File="PedigreeNrmInv.txt", Mat=Spec%PedNrmInvMat, Ija=Spec%PedNrmInvIja)
  !end if

  call cpu_time(EndTime)
  call AlphaRelateTitle
  call PrintElapsedTime(StartTime, EndTime)

end program

!###############################################################################
