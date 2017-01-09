
!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaRelate.f90
!
! DESCRIPTION:
!> @brief    Relationships
!
!> @details  Calculate relationships among individuals from different sources, i.e.,
!!           pedigree, marker genotypes, or haplotypes.
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date     December 19, 2016
!
!> @version  0.0.1 (alpha)
!
!-------------------------------------------------------------------------------

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
  call Spec%Read(SpecFile=trim(SpecFile))

  write(STDOUT, "(a)") ""
  write(STDOUT, "(a)") " Processing data ..."
  call Data%Read(Spec=Spec)

  if (Spec%PedInbreeding) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree inbreeding ..."
    call Data%CalcPedInbreeding
    call Data%PedInbreeding%Write(File="PedigreeInbreeding.txt", OutputFormat=Spec%OutputFormat)
  end if

  if (Spec%PedNrm) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree NRM ..."
    call Data%CalcPedNrm(Spec=Spec)
    call Data%PedNrm%Write(File="PedigreeNrm.txt", OutputFormat=Spec%OutputFormat, Ija=Spec%PedNrmIja)
  end if

  if (Spec%PedNrmInv) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree NRM inverse ..."
    call Data%CalcPedNrmInv
    call Data%PedNrmInv%Write(File="PedigreeNrmInv.txt", OutputFormat=Spec%OutputFormat, Ija=Spec%PedNrmInvIja)
  end if

  ! if (Spec%GenInbreeding) then
  !   write(STDOUT, "(a)") ""
  !   write(STDOUT, "(a)") " Calculating genotype inbreeding ..."
  !   call Data%CalcGenInbreeding
  !   call Data%WriteGenInbreeding(File="GenotypeInbreeding.txt", Spec=Spec)
  ! end if

  ! if (Spec%GenNrm) then
  !   write(STDOUT, "(a)") ""
  !   write(STDOUT, "(a)") " Calculating genotype NRM ..."
  !   call Data%CalcGenNrm(Spec=Spec)
  !   call Data%WriteGenNrm(File="GenotypeNrm.txt", Spec=Spec)
  ! end if

  ! if (Spec%GenNrmInv) then
  !   write(STDOUT, "(a)") ""
  !   write(STDOUT, "(a)") " Calculating genotype NRM inverse ..."
  !   call Data%CalcGenNrmInv
  !   call Data%WriteGenNrmInv(File="GenotypeNrmInv.txt", Spec=Spec)
  ! end if

  call cpu_time(EndTime)
  call AlphaRelateTitle
  call PrintElapsedTime(StartTime, EndTime)

end program

!###############################################################################
