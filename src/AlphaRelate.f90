#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#define DASH "\"
#define COPY "copy"
#define MD "md"
#define RMDIR "RMDIR /S /Q"
#define RM "del"
#define RENAME "MOVE /Y"
#define SH "BAT"
#define EXE ".exe"
#define NULL " >NUL"

#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#define DASH "/"
#define COPY "cp"
#define MD "mkdir"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#define SH "sh"
#define EXE ""
#define NULL ""

#endif

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
!> @date     January 14, 2016
!
!> @version  0.1.0 (alpha)
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
  logical :: InversionSucceded
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
  if (Spec%PedigreeGiven) then
    call Data%RecPed%Write(File=trim(Spec%OutputBasename)//trim(Spec%PedigreeFile)//"_Recoded.txt")
  end if

  if (Spec%PedInbreeding) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree inbreeding ..."
    call Data%CalcPedInbreeding
    call Data%PedInbreeding%Write(File=trim(Spec%OutputBasename)//"PedigreeInbreeding.txt", OutputFormat=Spec%OutputFormat)
  end if

  if (Spec%PedNrm) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree NRM ..."
    call Data%CalcPedNrm(Spec=Spec)
    call Data%PedNrm%Write(File=trim(Spec%OutputBasename)//"PedigreeNrm.txt", OutputFormat=Spec%OutputFormat, Ija=Spec%PedNrmIja)
  end if

  if (Spec%PedNrmInv) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating pedigree NRM inverse ..."
    call Data%CalcPedNrmInv
    call Data%PedNrmInv%Write(File=trim(Spec%OutputBasename)//"PedigreeNrmInv.txt", OutputFormat=Spec%OutputFormat, Ija=Spec%PedNrmInvIja)
  end if

  if (Spec%GenNrm) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating genotype NRM ..."
    call Data%CalcGenNrm(Spec=Spec)
    call Data%GenNrm%Write(File=trim(Spec%OutputBasename)//"GenotypeNrm.txt", OutputFormat=Spec%OutputFormat, Ija=Spec%GenNrmIja)
    if ((trim(Spec%GenNrmType) == "vanraden1" .or. &
         trim(Spec%GenNrmType) == "vanraden2" .or. &
         trim(Spec%GenNrmType) == "yang")    .and. &
         .not. Spec%AlleleFreqGiven) then
      call Data%AlleleFreq%Write(File=trim(Spec%OutputBasename)//"AlleleFreq.txt")
    end if
  end if

  if (Spec%GenInbreeding) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating genotype inbreeding ..."
    call Data%CalcGenInbreeding(Spec=Spec)
    call Data%GenInbreeding%Write(File=trim(Spec%OutputBasename)//"GenotypeInbreeding.txt", OutputFormat=Spec%OutputFormat)
  end if

  if (Spec%GenNrmInv) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating genotype NRM inverse ..."
    call Data%CalcGenNrmInv(Spec=Spec, Info=InversionSucceded)
    if (InversionSucceded) then
      call Data%GenNrmInv%Write(File=trim(Spec%OutputBasename)//"GenotypeNrmInv.txt", OutputFormat=Spec%OutputFormat)
    else
      write(STDERR, "(a)") " ERROR: Inversion of genotype NRM failed"
      write(STDERR, "(a)") " "
      stop 1
    end if
  end if

  if (Spec%HapIbdNrm) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating haplotype IBD NRM ..."
    call Data%CalcHapIbdNrm(Spec=Spec)
    call Data%HapIbdNrm%Write(File=trim(Spec%OutputBasename)//"HaplotypeIbdNrm.txt", OutputFormat=Spec%OutputFormat, Ija=Spec%HapIbdNrmIja)
  end if

  if (Spec%HapIbdInbreeding) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating haplotype IBD inbreeding ..."
    call Data%CalcHapIbdInbreeding(Spec=Spec)
    call Data%HapIbdInbreeding%Write(File=trim(Spec%OutputBasename)//"HaplotypeIbdInbreeding.txt", OutputFormat=Spec%OutputFormat)
  end if

  if (Spec%HapIbdNrmInv) then
    write(STDOUT, "(a)") ""
    write(STDOUT, "(a)") " Calculating haplotype IBD NRM inverse ..."
    call Data%CalcHapIbdNrmInv(Spec=Spec, Info=InversionSucceded)
    if (InversionSucceded) then
      call Data%HapIbdNrmInv%Write(File=trim(Spec%OutputBasename)//"HaplotypeIbdNrmInv.txt", OutputFormat=Spec%OutputFormat)
    else
      write(STDERR, "(a)") " ERROR: Inversion of haplotype IBD NRM failed"
      write(STDERR, "(a)") " "
      stop 1
    end if
  end if

  call cpu_time(EndTime)
  call AlphaRelateTitle
  call PrintElapsedTime(StartTime, EndTime)

end program

!###############################################################################
