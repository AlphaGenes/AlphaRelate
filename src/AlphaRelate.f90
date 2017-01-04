
!###############################################################################

program AlphaRelate
  use ISO_Fortran_env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use ConstantModule, only : FILELENGTH
  use AlphaHouseMod, only : PrintElapsedTime
  use AlphaRelateModule

  implicit none

  integer(int32) :: nArg
  real(real32) :: StartTime, EndTime
  character(len=FILELENGTH) :: SpecFile
  type(AlphaRelateSpec) :: Spec
  type(AlphaRelateData) :: Data

  call cpu_time(StartTime)
  call AlphaRelateTitle

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " Processing specifications ..."
  nArg = command_argument_count()
  if (nArg > 0) then
    call get_command_argument(1, SpecFile)
  else
    SpecFile = "AlphaRelateSpec.txt"
  end if
  write(STDOUT, "(2a)") " Using specification file: ", trim(SpecFile)
  Spec = AlphaRelateSpec(SpecFile=trim(SpecFile))

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " Processing data ..."
  Data = AlphaRelateData(Spec=Spec)

  if (Spec%PedInbreeding) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree inbreeding ..."
    call Data%CalcPedInbreeding
    call Data%WritePedInbreeding(File="PedigreeInbreeding.txt", Spec=Spec)
  end if

  if (Spec%PedNrm) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree NRM ..."
    call Data%CalcPedNrm(Spec=Spec)
    call Data%WritePedNrm(File="PedigreeNrm.txt", Spec=Spec)
  end if

  if (Spec%PedNrmInv) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree NRM inverse ..."
    call Data%CalcPedNrmInv
    call Data%WritePedNrmInv(File="PedigreeNrmInv.txt", Spec=Spec)
  end if

  call cpu_time(EndTime)
  call AlphaRelateTitle
  call PrintElapsedTime(StartTime, EndTime)

end program

!###############################################################################
