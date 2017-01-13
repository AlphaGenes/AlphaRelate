
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaRelateModule.f90
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
! REVISION HISTORY:
! 2016-12-19 GGorjanc - Initial setup as portable/documented module
!
!-------------------------------------------------------------------------------
module AlphaRelateModule
  use ISO_Fortran_env, STDIN => input_unit, STDOUT => output_unit, STDERR => error_unit
  use ConstantModule, only : FILELENGTH, SPECOPTIONLENGTH, IDLENGTH, IDINTLENGTH,&
                             EMPTYID, MISSINGGENOTYPECODE
  use AlphaHouseMod, only : CountLines, Char2Int, Char2Double, Int2Char, Real2Char,&
                            ParseToFirstWhitespace, SplitLineIntoTwoParts, ToLower, FindLoc
  use PedigreeModule, only : PedigreeHolder, RecodedPedigreeArray, MakeRecodedPedigreeArray
  use GenotypeModule, only : Genotype
  use HaplotypeModule, only : Haplotype
  use Blas95
  use Lapack95
  use F95_precision

  implicit none

  private
  ! Types
  public :: AlphaRelateTitle, AlphaRelateSpec, AlphaRelateData, IndSet, Inbreeding
  public :: Nrm, AlleleFreq, LocusWeight, GenotypeArray, HaplotypeArray
  ! Functions
  public :: MatchId, PedInbreedingRecursive, PedInbreedingMeuwissenLuo, PedNrm, PedNrmTimesVector, PedNrmInv

  !> @brief Set of individuals holder
  type IndSet
    integer(int32)                                     :: nInd
    character(len=IDLENGTH), allocatable, dimension(:) :: OriginalId
    integer(int32), allocatable, dimension(:)          :: Id
    contains
      procedure :: Init    => InitIndSet
      procedure :: Destroy => DestroyIndSet
      procedure :: MatchId => MatchIdIndSet
      procedure :: Write   => WriteIndSet
      procedure :: Read    => ReadIndSet
  end type

  !> @brief Inbreeding holder
  type Inbreeding
    ! @todo howto to extend the IndSet class and inherit some of the methods (Init and MatchId are pretty much the same)
    integer(int32)                                     :: nInd
    character(len=IDLENGTH), allocatable, dimension(:) :: OriginalId
    integer(int32), allocatable, dimension(:)          :: Id
    real(real64), allocatable, dimension(:)            :: Inb
    contains
      procedure :: Init    => InitInbreeding
      procedure :: Destroy => DestroyInbreeding
      procedure :: MatchId => MatchIdInbreeding
      procedure :: Write   => WriteInbreeding
      procedure :: Read    => ReadInbreeding
  end type

  !> @brief Numerator relationship (or its inverse) holder
  type Nrm
    ! @todo howto to extend the IndSet class and inherit some of the methods (Init and MatchId are pretty much the same)
    integer(int32)                                     :: nInd
    character(len=IDLENGTH), allocatable, dimension(:) :: OriginalId
    integer(int32), allocatable, dimension(:)          :: Id
    real(real64), allocatable, dimension(:, :)         :: Nrm
    ! @todo use packed storage?
    ! - see https://software.intel.com/en-us/node/468672
    ! - see https://software.intel.com/en-us/node/471382#C62A5095-C2EE-4AAD-AAA6-589230521A55
    ! @todo: create dense and sparse version?
    contains
      procedure :: Init    => InitNrm
      procedure :: Destroy => DestroyNrm
      procedure :: MatchId => MatchIdNrm
      procedure :: Write   => WriteNrm
      procedure :: Read    => ReadNrm
  end type

  !> @brief Genotype data set holder
  type GenotypeArray
    ! @todo howto to extend the IndSet class and inherit some of the methods (Init and MatchId are pretty much the same)
    integer(int32)                                     :: nInd
    character(len=IDLENGTH), allocatable, dimension(:) :: OriginalId
    integer(int32), allocatable, dimension(:)          :: Id
    integer(int32)                                     :: nLoc
    type(Genotype), allocatable, dimension(:)          :: Genotype
    real(real64), allocatable, dimension(:, :)         :: GenotypeReal
    contains
      procedure :: Init                       => InitGenotypeArray
      procedure :: Destroy                    => DestroyGenotypeArray
      procedure :: MatchId                    => MatchIdGenotypeArray
      procedure :: Write                      => WriteGenotypeArray
      procedure :: Read                       => ReadGenotypeArray
      procedure :: WriteReal                  => WriteGenotypeArrayReal
      procedure :: ReadReal                   => ReadGenotypeArrayReal
      procedure :: MakeGenotypeReal
      procedure :: CenterGenotypeReal
      procedure :: CenterAndScaleGenotypeReal
      procedure :: WeightGenotypeReal
  end type

  !> @brief Haplotype data set holder
  type HaplotypeArray
    ! @todo howto to extend the IndSet class and inherit some of the methods (Init and MatchId are pretty much the same)
    integer(int32)                                     :: nInd
    character(len=IDLENGTH), allocatable, dimension(:) :: OriginalId
    integer(int32), allocatable, dimension(:)          :: Id
    integer(int32)                                     :: nLoc
    type(Haplotype), allocatable, dimension(:)         :: Haplotype
    real(real64), allocatable, dimension(:, :)         :: HaplotypeReal
    contains
      ! procedure :: Init                       => InitHaplotypeArray
      ! procedure :: Destroy                    => DestroyHaplotypeArray
      ! procedure :: MatchId                    => MatchIdHaplotypeArray
      ! procedure :: Write                      => WriteHaplotypeArray
      ! procedure :: Read                       => ReadHaplotypeArray
      ! procedure :: WriteReal                  => WriteHaplotypeArrayReal
      ! procedure :: ReadReal                   => ReadHaplotypeArrayReal
      ! procedure :: MakeHaplotypeReal
      ! procedure :: CenterHaplotypeReal
      ! procedure :: CenterAndScaleHaplotypeReal
  end type

  ! @brief Allele frequencies
  type AlleleFreq
    integer(int32)                          :: nLoc
    real(real64), allocatable, dimension(:) :: Value
    contains
      procedure :: Init    => InitAlleleFreq
      procedure :: Destroy => DestroyAlleleFreq
      procedure :: Write   => WriteAlleleFreq
      procedure :: Read    => ReadAlleleFreq
  end type

  ! @brief Locus weights
  type LocusWeight
    integer(int32)                          :: nLoc
    real(real64), allocatable, dimension(:) :: Value
    contains
      procedure :: Init    => InitLocusWeight
      procedure :: Destroy => DestroyLocusWeight
      procedure :: Write   => WriteLocusWeight
      procedure :: Read    => ReadLocusWeight
  end type

  !> @brief AlphaRelate specifications
  type AlphaRelateSpec
    character(len=FILELENGTH) :: SpecFile, PedigreeFile, GenotypeFile, HaplotypeFile
    character(len=FILELENGTH) :: PedNrmSubsetFile, AlleleFreqFile, LocusWeightFile!, PedNrmOldFile
    character(len=FILELENGTH) :: OutputBasename
    character(len=SPECOPTIONLENGTH) :: OutputFormat, GenNrmType

    logical :: PedigreeGiven, GenotypeGiven, HaplotypeGiven
    logical :: PedNrmSubsetGiven, PedNrmOldGiven, AlleleFreqGiven, AlleleFreqFixed, LocusWeightGiven

    logical :: PedInbreeding, PedNrm, PedNrmIja, PedNrmInv, PedNrmInvIja
    logical :: GenInbreeding, GenNrm, GenNrmIja, GenNrmInv, GenNrmInvIja
    ! logical :: HapInbreeding, HapNrm, HapNrmIja, HapNrmInv, HapNrmInvIja
    logical :: FudgeGenNrmDiag, BlendGenNrmWithPedNrm!, FudgeHapNrmDiag, BlendHapNrmWithPedNrm

    integer(int32) :: nLoc

    real(real64) :: AlleleFreqFixedValue
    real(real64) :: FudgeGenNrmDiagValue, BlendGenNrmWithPedNrmFactor(2)!, FudgeHapNrmDiagValue, BlendHapNrmWithPedNrmFactor(2)
    contains
      procedure :: Init => InitAlphaRelateSpec
      procedure :: Read => ReadAlphaRelateSpec
  end type

  !> @brief AlphaRelate data
  type AlphaRelateData
    ! Pedigree stuff
    type(RecodedPedigreeArray) :: RecPed
    type(IndSet)               :: PedNrmSubset
    type(Inbreeding)           :: PedInbreeding
    type(Nrm)                  :: PedNrm
    type(Nrm)                  :: PedNrmInv
    ! Genome stuff
    type(AlleleFreq)           :: AlleleFreq
    type(LocusWeight)          :: LocusWeight
    ! Genotype stuff
    type(GenotypeArray)        :: Gen
    type(Inbreeding)           :: GenInbreeding
    type(Nrm)                  :: GenNrm
    type(Nrm)                  :: GenNrmInv
    ! Haplotype stuff
    type(HaplotypeArray)       :: Hap
    ! type(Inbreeding)           :: HapInbreeding
    ! type(Nrm)                  :: HapNrm
    ! type(Nrm)                  :: HapNrmInv

    contains
      procedure :: Read              => ReadAlphaRelateData
      procedure :: Destroy           => DestroyAlphaRelateData
      procedure :: Write             => WriteAlphaRelateData
      procedure :: CalcPedInbreeding => CalcPedInbreedingAlphaRelateData
      procedure :: CalcPedNrm        => CalcPedNrmAlphaRelateData
      procedure :: CalcPedNrmInv     => CalcPedNrmInvAlphaRelateData
      procedure :: CalcAlleleFreq    => CalcAlleleFreqAlphaRelateData
      procedure :: CalcGenInbreeding => CalcGenInbreedingAlphaRelateData
      procedure :: CalcGenNrm        => CalcGenNrmAlphaRelateData
      procedure :: CalcGenNrmInv     => CalcGenNrmInvAlphaRelateData
  end type

  contains

    !###########################################################################

    subroutine AlphaRelateTitle ! not pure due to IO
      implicit none
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            *     AlphaRelate     *                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                                                                              "
      write(STDOUT, "(a)") "           Software for calculating relationships among individuals           "
      write(STDOUT, "(a)") "                       http://AlphaGenes.Roslin.ed.ac.uk                      "
      write(STDOUT, "(a)") "                                 No liability                                 "
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "                       Commit:   "//TOSTRING(COMMIT),"                        "
      write(STDOUT, "(a)") "                       Compiled: "//__DATE__//", "//__TIME__
      write(STDOUT, "(a)") ""
    end subroutine

    !###########################################################################

    ! IndSet type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  IndSet constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine InitIndSet(This, nInd, OriginalId, OriginalIdSuperset, Skip)
        implicit none

        ! Arguments
        class(IndSet), intent(out) :: This                                     !< @return IndSet holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(nInd)      !< Original Id of individuals in the      set (note that this should not have 0th margin)
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip


        ! Other
        integer(int32) :: Ind, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        This%nInd = nInd

        if (allocated(This%OriginalId)) then
          deallocate(This%OriginalId)
        end if
        allocate(This%OriginalId(0:nInd))
        This%OriginalId(0) = EMPTYID
        if (present(OriginalId)) then
          This%OriginalId(1:nInd) = OriginalId
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        This%Id(0) = 0
        if (present(OriginalIdSuperset)) then
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id(1:nInd) = [(Ind, Ind = 1, nInd)]
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  IndSet destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine DestroyIndSet(This)
        implicit none
        class(IndSet), intent(inout) :: This !< @return IndSet holder
        This%nInd = 0
        if (allocated(This%OriginalId)) then
         deallocate(This%OriginalId)
        end if
        if (allocated(This%Id)) then
         deallocate(This%Id)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Find location of a set of original Id in the superset of original Id
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine MatchIdIndSet(This, OriginalIdSuperset, Skip)
        implicit none

        ! Arguments
        class(IndSet), intent(inout) :: This                         !< @return IndSet holder
        character(len=IDLENGTH), intent(in) :: OriginalIdSuperset(:) !< The superset of individual ids
        integer(int32), intent(in), optional :: Skip                 !< How many elements of OriginalIdSuperset to skip

        ! Other
        integer(int32) :: n, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        n = This%nInd
        This%Id(0) = 0
        This%Id(1:n) = MatchId(IdSet=This%OriginalId(1:n),& ! to handle "0th margin"
                               IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Find location of a set of original Id in the superset of original Id
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure function MatchId(IdSet, IdSuperset) result(Result)
        implicit none

        ! Arguments
        character(len=IDLENGTH), intent(in) :: IdSet(:)      !< A        set of individual ids
        character(len=IDLENGTH), intent(in) :: IdSuperset(:) !< The superset of individual ids
        integer(int32), allocatable, dimension(:) :: Result  !< @return Locations

        ! Other
        integer(int32) :: n, i

        n = size(IdSet)
        allocate(Result(n))
        do i = 1, n
          Result(i) = FindLoc(Val=IdSet(i), Vec=IdSuperset)
        end do
      end function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief   Write IndSet to a file or stdout
      !> @author  Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date    January 4, 2017
      !-------------------------------------------------------------------------
      subroutine WriteIndSet(This, File) ! not pure due to IO
        implicit none
        class(IndSet), intent(in) :: This              !< IndSet holder
        character(len=*), intent(in), optional :: File !< File that will hold a set of Original Id and internal integer sequence

        integer(int32) :: Unit, Ind
        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if
        do Ind = 1, This%nInd
          write(Unit, "(a, i)") This%OriginalId(Ind), This%Id(Ind)
        end do
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read IndSet from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      subroutine ReadIndSet(This, File) ! not pure due to IO
        implicit none
        class(IndSet), intent(out) :: This   !< @return IndSet holder
        character(len=*), intent(in) :: File !< File that holds a set of Original Id (internal integer sequence is not read)

        integer(int32) :: nInd, Ind, Unit

        nInd = CountLines(File)
        call This%Init(nInd=nInd)
        open(newunit=Unit, file=trim(File), action="read", status="old")
        do Ind = 1, nInd
          read(Unit, *) This%OriginalId(Ind)
        end do
        close(Unit)
      end subroutine

      !#########################################################################

    !###########################################################################

    ! Inbreeding type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Inbreeding constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine InitInbreeding(This, nInd, OriginalId, OriginalIdSuperset, Skip, InbInput)
        implicit none

        ! Arguments
        class(Inbreeding), intent(out) :: This                                 !< @return Inbreeding holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(:)         !< Original Id of individuals in the      set (note that this should not have 0th margin)
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip
        real(real64), intent(in), optional :: InbInput(nInd)                   !< Inbreeding coefficients (note that this should not have 0th margin)

        ! Other
        integer(int32) :: Ind, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        This%nInd = nInd

        if (allocated(This%OriginalId)) then
          deallocate(This%OriginalId)
        end if
        allocate(This%OriginalId(0:nInd))
        This%OriginalId(0) = EMPTYID
        if (present(OriginalId)) then
          This%OriginalId(1:nInd) = OriginalId
        else
          This%OriginalId(1:nInd) = EMPTYID
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        This%Id(0) = 0
        if (present(OriginalIdSuperset)) then
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id(1:nInd) = [(Ind, Ind = 1, nInd)]
        end if

        if (allocated(This%Inb)) then
          deallocate(This%Inb)
        end if
        allocate(This%Inb(0:nInd))
        This%Inb(0) = -1.0d0
        if (present(InbInput)) then
          This%Inb(1:nInd) = InbInput
        else
          This%Inb(1:nInd) = 0.0d0
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Inbreeding destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine DestroyInbreeding(This)
        implicit none
        class(Inbreeding), intent(inout) :: This !< @return Inbreeding holder
        This%nInd = 0
        if (allocated(This%OriginalId)) then
         deallocate(This%OriginalId)
        end if
        if (allocated(This%Id)) then
         deallocate(This%Id)
        end if
        if (allocated(This%Inb)) then
         deallocate(This%Inb)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Find location of a set of original Id in the superset of original Id
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine MatchIdInbreeding(This, OriginalIdSuperset, Skip)
        implicit none

        ! Arguments
        class(Inbreeding), intent(inout) :: This                     !< @return Inbreeding holder
        character(len=IDLENGTH), intent(in) :: OriginalIdSuperset(:) !< The superset of individual ids
        integer(int32), intent(in), optional :: Skip                 !< How many elements of OriginalIdSuperset to skip

        ! Other
        integer(int32) :: n, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        n = This%nInd
        This%Id(0) = 0
        This%Id(1:n) = MatchId(IdSet=This%OriginalId(1:n),& ! to handle "0th margin"
                               IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write inbreeding to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine WriteInbreeding(This, File, OutputFormat) ! not pure due to IO
        implicit none
        class(Inbreeding), intent(in) :: This                  !< Inbreeding holder
        character(len=*), intent(in), optional :: File         !< File that will hold Original Id and inbreeding coefficients
        character(len=*), intent(in), optional :: OutputFormat !< Format for inbreeding coefficients, default is "f"

        integer(int32) :: Unit, Ind
        character(len=:), allocatable :: OutputFormatInternal, Fmt

        if (present(OutputFormat)) then
          OutputFormatInternal = OutputFormat
        else
          OutputFormatInternal = "f"
        endif

        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if
        Fmt = "(a"//Int2Char(IDLENGTH)//", "//OutputFormatInternal//")"
        do Ind = 1, This%nInd
          write(Unit, Fmt) This%OriginalId(Ind), This%Inb(Ind)
        end do
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read inbreeding from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine ReadInbreeding(This, File) ! not pure due to IO
        implicit none
        class(Inbreeding), intent(out) :: This !< @return Inbreeding holder
        character(len=*), intent(in) :: File   !< File that holds Original Id and pedigree inbreeding

        integer(int32) :: nInd, Ind, Unit

        nInd = CountLines(File)
        call This%Init(nInd=nInd)
        open(newunit=Unit, file=trim(File), action="read", status="old")
        do Ind = 1, nInd
          read(Unit, *) This%OriginalId(Ind), This%Inb(Ind)
        end do
        close(Unit)
      end subroutine

      !#########################################################################

    !###########################################################################

    ! Nrm type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Nrm constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine InitNrm(This, nInd, OriginalId, OriginalIdSuperset, Skip, NrmInput)
        implicit none

        ! Arguments
        class(Nrm), intent(out) :: This                                        !< @return Nrm holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(nInd)      !< Original Id of individuals in the      set (note that this should not have 0th margin)
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip
        real(real64), intent(in), optional :: NrmInput(nInd, nInd)             !< Relationship coefficients (note that this should not have 0th margin)

        ! Other
        integer(int32) :: Ind, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        This%nInd = nInd

        if (allocated(This%OriginalId)) then
          deallocate(This%OriginalId)
        end if
        allocate(This%OriginalId(0:nInd))
        This%OriginalId(0) = EMPTYID
        if (present(OriginalId)) then
          This%OriginalId(1:nInd) = OriginalId
        else
          This%OriginalId(1:nInd) = EMPTYID
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        This%Id(0) = 0
        if (present(OriginalIdSuperset)) then
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id(1:nInd) = [(Ind, Ind = 1, nInd)]
        end if

        if (allocated(This%Nrm)) then
          deallocate(This%Nrm)
        end if
        allocate(This%Nrm(0:nInd, 0:nInd))
        This%Nrm(0:nInd, 0) = 0.0d0
        This%Nrm(0, 0:nInd) = 0.0d0
        if (present(NrmInput)) then
          This%Nrm(1:nInd, 1:nInd) = NrmInput
        else
          This%Nrm(1:nInd, 1:nInd) = 0.0d0
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Nrm destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine DestroyNrm(This)
        implicit none
        class(Nrm), intent(inout) :: This !< @return Nrm holder
        This%nInd = 0
        if (allocated(This%OriginalId)) then
         deallocate(This%OriginalId)
        end if
        if (allocated(This%Id)) then
         deallocate(This%Id)
        end if
        if (allocated(This%Nrm)) then
         deallocate(This%Nrm)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Find location of a set of original Id in the superset of original Id
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine MatchIdNrm(This, OriginalIdSuperset, Skip)
        implicit none

        ! Arguments
        class(Nrm), intent(inout) :: This                            !< @return Nrm holder
        character(len=IDLENGTH), intent(in) :: OriginalIdSuperset(:) !< The superset of individual ids
        integer(int32), intent(in), optional :: Skip                 !< How many elements of OriginalIdSuperset to skip

        ! Other
        integer(int32) :: n, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        n = This%nInd
        This%Id(0) = 0
        This%Id(1:n) = MatchId(IdSet=This%OriginalId(1:n),& ! to handle "0th margin"
                               IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write Nrm to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine WriteNrm(This, File, Ija, OutputFormat) ! not pure due to IO
        implicit none
        class(Nrm), intent(in) :: This                         !< Nrm Holder
        character(len=*), intent(in), optional :: File         !< File that will hold Original Id and NRM
        logical, intent(in), optional :: Ija                   !< Write in sparse ija format
        character(len=*), intent(in), optional :: OutputFormat !< Format for relationship values, default is "f"

        integer(int32) :: Unit, Unit2, Ind1, Ind2
        character(len=:), allocatable :: Fmt, OutputFormatInternal
        logical :: IjaInternal

        if (present(Ija)) then
          IjaInternal = Ija
        else
          IjaInternal = .false.
        end if

        if (present(OutputFormat)) then
          OutputFormatInternal = OutputFormat
        else
          OutputFormatInternal = "f"
        end if

        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if
        if (IjaInternal) then
          ! Original Ids
          if (present(File)) then
            open(newunit=Unit2, file=trim(File)//"_IdMap", action="write", status="unknown")
          else
            Unit2 = STDOUT
          end if
          Fmt = "(i"//Int2Char(IDINTLENGTH)//", a1, a"//Int2Char(IDLENGTH)//")"
          do Ind1 = 1, This%nInd
            write(Unit2, Fmt) Ind1, "", This%OriginalId(Ind1)
          end do
          if (present(File)) then
            close(Unit2)
          end if
          ! No. of individuals
          write(Unit, "(i)") This%nInd
          ! Triplets
          Fmt = "(2i"//Int2Char(IDINTLENGTH)//", "//OutputFormatInternal//")"
          do Ind1 = 1, This%nInd
            do Ind2 = Ind1, This%nInd
              if (abs(This%Nrm(Ind2, Ind1)) .gt. 0.d0) then
                write(Unit, Fmt) Ind2, Ind1, This%Nrm(Ind2, Ind1)
              end if
            end do
          end do
        else
          Fmt = "(a"//Int2Char(IDLENGTH)//", "//Int2Char(This%nInd)//OutputFormatInternal//")"
          do Ind1 = 1, This%nInd
            write(Unit, Fmt) This%OriginalId(Ind1), This%Nrm(1:This%nInd, Ind1)
          end do
        end if
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read Nrm from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine ReadNrm(This, File, Ija) ! not pure due to IO
        implicit none
        class(Nrm), intent(out) :: This      !< @return Nrm holder
        character(len=*), intent(in) :: File !< File that holds Original Id and NRM
        logical, intent(in) :: Ija           !< Read from a sparse ija format?

        integer(int32) :: n, Line, nLine, Unit, Unit2, Ind1, Ind2
        character(len=:), allocatable :: Fmt

        nLine = CountLines(File)

        open(newunit=Unit, file=trim(File), action="read", status="old")
        if (Ija) then
          ! No. of individuals
          read(Unit, *) n
          call This%Init(nInd=n)
          open(newunit=Unit2, file=trim(File)//"_IdMap", action="read", status="old")
          Fmt = "(i"//Int2Char(IDINTLENGTH)//", a"//Int2Char(IDLENGTH)//")"
          do Ind1 = 1, This%nInd
            read(Unit2, *) Ind2, This%OriginalId(Ind1) ! Ind2 just placeholder here
          end do
          close(Unit2)
          ! Triplets
          do Line = 1, (nLine - 1)
            read(Unit, *) Ind2, Ind1, This%Nrm(Ind2, Ind1)
            This%Nrm(Ind1,Ind2) = This%Nrm(Ind2, Ind1)
          end do
        else
          n = nLine
          call This%Init(nInd=n)
          do Ind1 = 1, n
            read(Unit, *) This%OriginalId(Ind1), This%Nrm(1:n, Ind1)
          end do
        end if
        close(Unit)
      end subroutine

      !#########################################################################

    !###########################################################################

    ! AlleleFreq type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  AlleleFreq constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine InitAlleleFreq(This, nLoc, AlleleFreqInput)
        implicit none

        ! Arguments
        class(AlleleFreq), intent(out) :: This                      !< @return AlleleFreq holder
        integer(int32), intent(in) :: nLoc                          !< Number of loci
        real(real64), intent(in), optional :: AlleleFreqInput(nLoc) !< Allele frequencies

        ! Init numbers
        This%nLoc = nLoc

        ! Init Value
        if (allocated(This%Value)) then
          deallocate(This%Value)
        end if
        allocate(This%Value(This%nLoc))
        if (present(AlleleFreqInput)) then
          This%Value = AlleleFreqInput
        else
          This%Value = 0.0d0
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  AlleleFreq destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine DestroyAlleleFreq(This)
        implicit none
        class(AlleleFreq), intent(inout) :: This !< @return AlleleFreq holder
        This%nLoc = 0
        if (allocated(This%Value)) then
         deallocate(This%Value)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write allele freqs to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 10, 2016
      !-------------------------------------------------------------------------
      subroutine WriteAlleleFreq(This, File) ! not pure due to IO
        implicit none
        class(AlleleFreq), intent(in) :: This           !< AlleleFre Holder
        character(len=*), intent(in), optional :: File !< File that will hold allele freqs

        integer(int32) :: Unit, Loc

        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if
        do Loc = 1, This%nLoc
          write(Unit, "(i, f)") Loc, This%Value(Loc)
        end do
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read allele freqs from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 10, 2016
      !-------------------------------------------------------------------------
      subroutine ReadAlleleFreq(This, File, nLoc) ! not pure due to IO
        implicit none
        class(AlleleFreq), intent(out) :: This       !< AlleleFreq Holder
        character(len=*), intent(in) :: File         !< File that holds allele freqs
        integer(int32), intent(in), optional :: nLoc !< Number of loci

        integer(int32) :: Unit, LocLoop, Loc

        if (present(nLoc)) then
          call This%Init(nLoc=nLoc)
        else
          call This%Init(nLoc=CountLines(File))
        end if
        open(newunit=Unit, file=trim(File), action="read", status="old")
        do LocLoop = 1, This%nLoc
          read(Unit, *) Loc, This%Value(Loc)
        end do
        close(Unit)
      end subroutine

      !#########################################################################

    !###########################################################################

    ! Locus weights type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  LocusWeight constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine InitLocusWeight(This, nLoc, LocusWeightInput)
        implicit none

        ! Arguments
        class(LocusWeight), intent(out) :: This                      !< @return LocusWeight holder
        integer(int32), intent(in) :: nLoc                           !< Number of loci
        real(real64), intent(in), optional :: LocusWeightInput(nLoc) !< Locus weights

        ! Init numbers
        This%nLoc = nLoc

        ! Init Value
        if (allocated(This%Value)) then
          deallocate(This%Value)
        end if
        allocate(This%Value(This%nLoc))
        if (present(LocusWeightInput)) then
          This%Value = LocusWeightInput
        else
          This%Value = 0.0d0
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  LocusWeight destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine DestroyLocusWeight(This)
        implicit none
        class(LocusWeight), intent(inout) :: This !< @return LocusWeight holder
        This%nLoc = 0
        if (allocated(This%Value)) then
         deallocate(This%Value)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write locus weights to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 10, 2016
      !-------------------------------------------------------------------------
      subroutine WriteLocusWeight(This, File) ! not pure due to IO
        implicit none
        class(LocusWeight), intent(in) :: This         !< LocusWeight Holder
        character(len=*), intent(in), optional :: File !< File that will hold loucs weights

        integer(int32) :: Unit, Loc

        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if
        do Loc = 1, This%nLoc
          write(Unit, "(i, f)") Loc, This%Value(Loc)
        end do
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read locus weights from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 10, 2016
      !-------------------------------------------------------------------------
      subroutine ReadLocusWeight(This, File, nLoc) ! not pure due to IO
        implicit none
        class(LocusWeight), intent(out) :: This      !< LocusWeight Holder
        character(len=*), intent(in) :: File         !< File that holds locus weights
        integer(int32), intent(in), optional :: nLoc !< Number of loci

        integer(int32) :: Unit, LocLoop, Loc

        if (present(nLoc)) then
          call This%Init(nLoc=nLoc)
        else
          call This%Init(nLoc=CountLines(File))
        end if

        open(newunit=Unit, file=trim(File), action="read", status="old")
        do LocLoop = 1, This%nLoc
          read(Unit, *) Loc, This%Value(Loc)
        end do
        close(Unit)
      end subroutine

      !#########################################################################

    !###########################################################################

    ! GenotypeArray type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  GenotypeArray constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine InitGenotypeArray(This, nInd, nLoc, OriginalId, OriginalIdSuperset, Skip, IntegerInput, GenotypeInput)
        implicit none

        ! Arguments
        class(GenotypeArray), intent(out) :: This                              !< @return GenotypeArray holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        integer(int32), intent(in) :: nLoc                                     !< Number of loci in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(nInd)      !< Original Id of individuals in the      set (note that this should not have 0th margin)
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip
        integer(int8), intent(in), optional :: IntegerInput(nLoc, nInd)        !< Genotypes as an array of integers (note that this should not have 0th margin)
        type(Genotype), intent(in), optional :: GenotypeInput(nInd)            !< Genotypes as genotype type (note that this should not have 0th margin)

        ! Other
        integer(int32) :: Ind, Start
        integer(int8) :: TempInt(nLoc)
        type(Genotype) :: TempGeno

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        ! Init numbers
        This%nInd = nInd
        This%nLoc = nLoc

        ! Init OriginalId
        if (allocated(This%OriginalId)) then
          deallocate(This%OriginalId)
        end if
        allocate(This%OriginalId(0:nInd))
        This%OriginalId(0) = EMPTYID
        if (present(OriginalId)) then
          This%OriginalId(1:nInd) = OriginalId
        else
          This%OriginalId(1:nInd) = EMPTYID
        end if

        ! Init Id
        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        This%Id(0) = 0
        if (present(OriginalIdSuperset)) then
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id(1:nInd) = [(Ind, Ind = 1, nInd)]
        end if

        ! Init Genotype
        ! NOTE: this one is allocated here as it is the way we primarily store genotypes
        if (allocated(This%Genotype)) then
          deallocate(This%Genotype)
        end if
        allocate(This%Genotype(0:nInd))
        TempInt = MISSINGGENOTYPECODE
        TempGeno = Genotype(Geno=TempInt)
        This%Genotype(0) = TempGeno
        if (present(IntegerInput) .or. present(GenotypeInput)) then
          if (present(IntegerInput)) then
            do Ind = 1, nInd
              This%Genotype(Ind) = Genotype(Geno=IntegerInput(:, Ind))
            end do
          end if
          if (present(GenotypeInput)) then
            do Ind = 1, nInd
              This%Genotype(Ind) = GenotypeInput(Ind)
            end do
          end if
        else
          do Ind = 1, nInd
            This%Genotype(Ind) = TempGeno
          end do
        end if

        ! Init GenotypeReal
        ! NOTE: this one is not allocated here as it is not the way we primarily store genotypes and can be potentially quite memory consuming
        ! Use call This%MakeGenotypeReal() if needed
        if (allocated(This%GenotypeReal)) then
          deallocate(This%GenotypeReal)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  GenotypeArray destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine DestroyGenotypeArray(This)
        implicit none
        class(GenotypeArray), intent(inout) :: This !< @return GenotypeArray holder
        This%nInd = 0
        This%nLoc = 0
        if (allocated(This%OriginalId)) then
         deallocate(This%OriginalId)
        end if
        if (allocated(This%Id)) then
         deallocate(This%Id)
        end if
        if (allocated(This%Genotype)) then
         deallocate(This%Genotype)
        end if
        if (allocated(This%GenotypeReal)) then
         deallocate(This%GenotypeReal)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Find location of a set of original Id in the superset of original Id
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 4, 2017
      !-------------------------------------------------------------------------
      pure subroutine MatchIdGenotypeArray(This, OriginalIdSuperset, Skip)
        implicit none

        ! Arguments
        class(GenotypeArray), intent(inout) :: This                  !< @return GenotypeArray holder
        character(len=IDLENGTH), intent(in) :: OriginalIdSuperset(:) !< The superset of individual ids
        integer(int32), intent(in), optional :: Skip                 !< How many elements of OriginalIdSuperset to skip

        ! Other
        integer(int32) :: n, Start

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        n = This%nInd
        This%Id(0) = 0
        This%Id(1:n) = MatchId(IdSet=This%OriginalId(1:n),& ! to handle "0th margin"
                               IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write GenotypeArray genotypes to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2017
      !-------------------------------------------------------------------------
      subroutine WriteGenotypeArray(This, File) ! not pure due to IO
        implicit none
        class(GenotypeArray), intent(in) :: This       !< GenotypeArray Holder
        character(len=*), intent(in), optional :: File !< File that will hold Genotypes

        integer(int32) :: Unit, Ind
        character(len=:), allocatable :: Fmt

        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if

        Fmt = "(a"//Int2Char(IDLENGTH)//", "//Int2Char(This%nLoc)//"i2)"
        do Ind = 1, This%nInd
          write(Unit, Fmt) This%OriginalId(Ind), This%Genotype(Ind)%ToIntegerArray()
        end do
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write GenotypeArray genotypes (the GenotypeReal component) to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2017
      !-------------------------------------------------------------------------
      subroutine WriteGenotypeArrayReal(This, File) ! not pure due to IO
        implicit none
        class(GenotypeArray), intent(in) :: This       !< GenotypeArray Holder
        character(len=*), intent(in), optional :: File !< File that will hold Genotypes

        integer(int32) :: Unit, Ind
        character(len=:), allocatable :: Fmt

        if (present(File)) then
          open(newunit=Unit, file=trim(File), action="write", status="unknown")
        else
          Unit = STDOUT
        end if

        Fmt = "(a"//Int2Char(IDLENGTH)//", "//Int2Char(This%nLoc)//"f)"
        do Ind = 1, This%nInd
          write(Unit, Fmt) This%OriginalId(Ind), This%GenotypeReal(:, Ind)
        end do
        if (present(File)) then
          close(Unit)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read GenotypeArray from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2017
      !-------------------------------------------------------------------------
      subroutine ReadGenotypeArray(This, File, nLoc) ! not pure due to IO
        implicit none
        class(GenotypeArray), intent(out) :: This !< @return GenotypeArray holder
        character(len=*), intent(in) :: File      !< File that holds genotypes
        integer(int32), intent(in) :: nLoc        !< Number of loci to read in

        integer(int32) :: nInd, Ind, Unit
        integer(int8) :: Temp(nLoc)

        nInd = CountLines(File)
        call This%Init(nInd=nInd, nLoc=nLoc)
        open(newunit=Unit, file=trim(File), action="read", status="old")
        do Ind = 1, This%nInd
          read(Unit, *) This%OriginalId(Ind), Temp
          This%Genotype(Ind) = Genotype(Geno=Temp)
        end do
        close(Unit)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read GenotypeArray (the GenotypeReal component) from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2017
      !-------------------------------------------------------------------------
      subroutine ReadGenotypeArrayReal(This, File, nLoc) ! not pure due to IO
        implicit none
        class(GenotypeArray), intent(out) :: This !< @return GenotypeArray holder
        character(len=*), intent(in) :: File      !< File that holds genotypes
        integer(int32), intent(in) :: nLoc        !< Number of loci to read in

        integer(int32) :: nInd, Ind, Unit

        nInd = CountLines(File)
        call This%Init(nInd=nInd, nLoc=nLoc)
        allocate(This%GenotypeReal(nLoc, 0:nInd))
        This%GenotypeReal(:, 0) = 0.0d0
        open(newunit=Unit, file=trim(File), action="read", status="old")
        do Ind = 1, This%nInd
          read(Unit, *) This%OriginalId(Ind), This%GenotypeReal(:, Ind)
        end do
        close(Unit)
      end subroutine

      !#########################################################################
      !-------------------------------------------------------------------------
      !> @brief  Make genotypes real (build real64 array of genotypes)
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine MakeGenotypeReal(This)
        implicit none
        class(GenotypeArray), intent(inout) :: This !< @return GenotypeArray holder

        integer(int32) :: Ind

        if (allocated(This%GenotypeReal)) then
          deallocate(This%GenotypeReal)
        end if
        allocate(This%GenotypeReal(This%nLoc, 0:This%nInd))

        This%GenotypeReal(:, 0) = dble(MISSINGGENOTYPECODE)
        ! @todo could we handle missing genotypes here somehow?
        do Ind = 1, This%nInd
          This%GenotypeReal(:, Ind) = dble(This%Genotype(Ind)%ToIntegerArray())
        end do
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Center genotypes real
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine CenterGenotypeReal(This, AlleleFreq)
        implicit none
        class(GenotypeArray), intent(inout) :: This                  !< @return GenotypeArray holder
        real(real64), intent(in), dimension(This%nLoc) :: AlleleFreq !< Allele frequencies for centering

        integer(int32) :: Ind, Loc
        real(real64), allocatable, dimension(:) :: Mean

        if (.not. allocated(This%GenotypeReal)) then
          call This%MakeGenotypeReal
        end if

        allocate(Mean(This%nLoc))
        Mean = 2.0d0 * AlleleFreq
        This%GenotypeReal(:, 0) = 0.0d0
        ! @todo could we not assume that GenotypeReal has no missing values
        !       (it should have been cleaned prior to this program) and then
        !       we can avoid the ifs and simplify computatation a lot
        ! do Ind = 1, This%nInd
        !   This%GenotypeReal(:, Ind) = This%GenotypeReal(:, Ind) - AlleleFreq2
        ! end do
        do Ind = 1, This%nInd
          do Loc = 1, This%nLoc
            if ((This%GenotypeReal(Loc, Ind) .ge. 0.0d0) .and. (This%GenotypeReal(Loc, Ind) .le. 2.0d0)) then
              This%GenotypeReal(Loc, Ind) = This%GenotypeReal(Loc, Ind) - Mean(Loc)
            else
              This%GenotypeReal(Loc, Ind) = 0.0d0
            end if
          end do
        end do
        deallocate(Mean)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Center and scale genotypes real
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine CenterAndScaleGenotypeReal(This, AlleleFreq)
        implicit none
        class(GenotypeArray), intent(inout) :: This                  !< @return GenotypeArray holder
        real(real64), intent(in), dimension(This%nLoc) :: AlleleFreq !< Allele frequencies for centering and scaling

        integer(int32) :: Ind, Loc
        real(real64), allocatable, dimension(:) :: Mean, StDev

        if (.not. allocated(This%GenotypeReal)) then
          call This%MakeGenotypeReal
        end if

        allocate(Mean(This%nLoc))
        allocate(StDev(This%nLoc))
        Mean = 2.0d0 * AlleleFreq
        StDev = sqrt(Mean * (1.0d0 - AlleleFreq))

        This%GenotypeReal(:, 0) = 0.0d0
        ! @todo could we not assume that GenotypeReal has no missing values
        !       (it should have been cleaned prior to this program) and then
        !       we can avoid the ifs and simplify computatation a lot
        ! do Ind = 1, This%nInd
        !   This%GenotypeReal(:, Ind) = (This%GenotypeReal(:, Ind) - AlleleFreq2) / Scale
        ! end do
        do Ind = 1, This%nInd
          do Loc = 1, This%nLoc
            if (StDev(Loc) > tiny(StDev(Loc))) then ! to avoid dividing by zero
              if ((This%GenotypeReal(Loc, Ind) .ge. 0.0d0) .and. (This%GenotypeReal(Loc, Ind) .le. 2.0d0)) then
                This%GenotypeReal(Loc, Ind) = (This%GenotypeReal(Loc, Ind) - Mean(Loc)) / StDev(Loc)
              else
                This%GenotypeReal(Loc, Ind) = 0.0d0
              end if
            end if
          end do
        end do
        deallocate(Mean)
        deallocate(StDev)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Weight genotypes real
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine WeightGenotypeReal(This, Weight)
        implicit none
        class(GenotypeArray), intent(inout) :: This              !< @return GenotypeArray holder
        real(real64), intent(in), dimension(This%nLoc) :: Weight !< Locus weights

        integer(int32) :: Ind

        if (.not. allocated(This%GenotypeReal)) then
          call This%MakeGenotypeReal
        end if

        do Ind = 1, This%nInd
          This%GenotypeReal(:, Ind) = This%GenotypeReal(:, Ind) * Weight
        end do
      end subroutine

      !#########################################################################

    !###########################################################################

    ! AlphaRelateSpec type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  AlphaRelateSpec constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine InitAlphaRelateSpec(This)
        implicit none

        ! Arguments
        class(AlphaRelateSpec), intent(out) :: This !< @return AlphaRelateSpec holder

        ! Defaults
        This%SpecFile         = "None"
        This%PedigreeFile     = "None"
        This%PedNrmSubsetFile = "None"
        ! This%PedNrmOldFile    = "None"
        This%GenotypeFile     = "None"
        This%HaplotypeFile    = "None"
        This%AlleleFreqFile   = "None"
        This%LocusWeightFile  = "None"
        This%OutputBasename   = ""

        This%GenNrmType       = "None"
        This%OutputFormat     = "f"

        This%PedigreeGiven     = .false.
        This%PedNrmSubsetGiven = .false.
        This%PedNrmOldGiven    = .false.
        This%GenotypeGiven     = .false.
        This%HaplotypeGiven    = .false.
        This%AlleleFreqGiven   = .false.
        This%AlleleFreqFixed   = .false.
        This%LocusWeightGiven  = .false.

        This%PedInbreeding = .false.
        This%PedNrm        = .false.
        This%PedNrmIja     = .false.
        This%PedNrmInv     = .false.
        This%PedNrmInvIja  = .false.

        This%GenInbreeding         = .false.
        This%GenNrm                = .false.
        This%GenNrmIja             = .false.
        This%GenNrmInv             = .false.
        This%GenNrmInvIja          = .false.
        This%FudgeGenNrmDiag       = .false.
        This%BlendGenNrmWithPedNrm = .false.

        ! This%HapInbreeding      = .false.
        ! This%HapNrm             = .false.
        ! This%HapNrmIja          = .false.
        ! This%HapNrmInv          = .false.
        ! This%HapNrmInvIja       = .false.
        ! This%FudgeHapNrmDiag    = .false.
        ! This%BlendHapWithPedNrm = .false.

        This%nLoc             = 0

        This%AlleleFreqFixedValue        = 0.5d0
        This%FudgeGenNrmDiagValue        = 0.0d0
        This%BlendGenNrmWithPedNrmFactor = [1.0d0, 0.0d0]
        ! This%FudgeHapNrmDiagValue        = 0.0d0
        ! This%BlendHapNrmWithPedNrmFactor = [1.0d0, 0.0d0]
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read AlphaRelateSpec from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine ReadAlphaRelateSpec(This, SpecFile) ! not pure due to IO
        implicit none

        ! Arguments
        class(AlphaRelateSpec), intent(out) :: This !< @return AlphaRelateSpec holder
        character(len=*), intent(in) :: SpecFile    !< Spec file; when missing, a stub with defaults is created

        ! Other
        character(len=:), allocatable :: DumString
        character(len=SPECOPTIONLENGTH) :: Line
        character(len=SPECOPTIONLENGTH) :: First
        character(len=SPECOPTIONLENGTH), allocatable, dimension(:) :: Second

        integer(int32) :: SpecUnit, Stat

        ! Defaults
        call This%Init

        This%SpecFile = SpecFile
        open(newunit=SpecUnit, file=trim(This%SpecFile), action="read", status="old")

        Stat = 0
        ReadSpec: do while (Stat == 0)
          read(SpecUnit, "(a)", iostat=Stat) Line
          if (len_trim(Line) == 0) then
            cycle
          end if
          call SplitLineIntoTwoParts(trim(adjustl(Line)), First, Second)
          DumString = ParseToFirstWhitespace(First)
          ! @todo why (len_trim(Line) == 0)? if we use (len_trim(Line) == 0) above
          if (First(1:1) == "=" .or. len_trim(Line) == 0) then
            cycle
          else
            select case (ToLower(trim(DumString)))

              case ("outputbasename")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "none") then
                    write(STDOUT, "(a)") " Not using output basename"
                  else
                    write(This%OutputBasename, *) trim(adjustl(Second(1)))
                    This%OutputBasename = adjustl(This%OutputBasename)
                    write(STDOUT, "(2a)") " Using output basename: ", trim(This%OutputBasename)
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a string for OutputBasename, i.e., OutputBasename, AnalysisX"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("pedigreefile")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "none") then
                    write(STDOUT, "(a)") " Not using pedigree file"
                  else
                    This%PedigreeGiven = .true.
                    write(This%PedigreeFile, *) trim(adjustl(Second(1)))
                    This%PedigreeFile = adjustl(This%PedigreeFile)
                    write(STDOUT, "(2a)") " Using pedigree file: ", trim(This%PedigreeFile)
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a file for PedigreeFile, i.e., PedigreeFile, Pedigree.txt"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("genotypefile")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "none") then
                    write(STDOUT, "(a)") " Not using genotype file"
                  else
                    This%GenotypeGiven = .true.
                    write(This%GenotypeFile, *) trim(adjustl(Second(1)))
                    This%GenotypeFile = adjustl(This%GenotypeFile)
                    write(STDOUT, "(2a)") " Using genotype file: ", trim(This%GenotypeFile)
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a file for GenotypeFile, i.e., GenotypeFile, Genotype.txt"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              ! case ("haplotypefile")
              !   if (allocated(Second)) then
              !     if (ToLower(trim(adjustl(Second(1)))) == "none") then
              !       write(STDOUT, "(a)") " Not using haplotype file"
              !     else
              !       This%HaplotypeGiven = .true.
              !       write(This%HaplotypeFile, *) trim(adjustl(Second(1)))
              !       This%HaplotypeFile = adjustl(This%HaplotypeFile)
              !       write(STDOUT, "(2a)") " Using haplotype file: ", trim(This%HaplotypeFile)
              !     end if
              !   else
              !     write(STDERR, "(a)") " ERROR: Must specify a file for HaplotypeFile, i.e., HaplotypeFile, Haplotype.txt"
              !     write(STDERR, "(a)") ""
              !     stop 1
              !   end if

              case ("locusweightfile")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "none") then
                    write(STDOUT, "(a)") " Not using locus weights file"
                  else
                    This%LocusWeightGiven = .true.
                    write(This%LocusWeightFile, *) trim(adjustl(Second(1)))
                    This%LocusWeightFile = adjustl(This%LocusWeightFile)
                    write(STDOUT, "(2a)") " Using locus weight file: ", trim(This%LocusWeightFile)
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a file for LocusWeightFile, i.e., LocusWeightFile, LocusWeight.txt"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("allelefreqfile")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "none") then
                    write(STDOUT, "(a)") " Not using precalculated/fixed allele frequencies file"
                  else
                    This%AlleleFreqGiven = .true.
                    if (ToLower(trim(adjustl(Second(1)))) == "fixed") then
                      This%AlleleFreqFixed = .true.
                      if (size(Second) > 1) then
                        This%AlleleFreqFixedValue = Char2Double(trim(adjustl(Second(2))), "(f20.16)")
                      else
                        This%AlleleFreqFixedValue = 0.5d0
                      end if
                      write(STDOUT, "(2a)") " Using fixed allele frequency: ", Real2Char(This%AlleleFreqFixedValue, "(f6.4)")
                    else
                      write(This%AlleleFreqFile, *) trim(adjustl(Second(1)))
                      This%AlleleFreqFile = adjustl(This%AlleleFreqFile)
                      write(STDOUT, "(2a)") " Using allele frequencies file: ", trim(This%AlleleFreqFile)
                    end if
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a file for AlleleFreqFile, i.e., AlleleFreqFile, AlleleFreq.txt"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("numberofloci")
                if (allocated(Second)) then
                  This%nLoc = Char2Int(trim(adjustl(Second(1))))
                  write(STDOUT, "(a, i0)") " Number of loci: ", This%nLoc
                else
                  write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfLoci, i.e., NumberOfLoci, 10"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              ! case ("numberoftraits")
              !   if (allocated(Second)) then
              !     This%nTrait = Char2Int(trim(adjustl(Second(1))))
              !     write(STDOUT, "(a, i)") " Number of traits: ", This%nTrait
              !   else
              !     write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfTraits, i.e., NumberOfTraits, 1"
              !     write(STDERR, "(a)") ""
              !     stop 1
              !   end if

              case ("outputformat")
                if (allocated(Second)) then
                  write(This%OutputFormat, *) trim(adjustl(Second(1)))
                  This%OutputFormat = adjustl(This%OutputFormat)
                  write(STDOUT, "(2a)") " Output precision: ", trim(This%OutputFormat)
                else
                  write(STDERR, "(a)") " ERROR: Must specify Fortran format for OutputFormat, i.e., OutputFormat, f16.8"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("pedinbreeding")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "yes") then
                    This%PedInbreeding = .true.
                    write(STDOUT, "(a)") " Calculate pedigree inbreeding: Yes"
                  else
                    write(STDOUT, "(a)") " Calculate pedigree inbreeding: No"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify Yes/No for PedInbreeding, i.e., PedInbreeding, Yes"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("pednrm")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "yes") then
                    This%PedNrm = .true.
                    write(STDOUT, "(a)") " Calculate pedigree NRM: Yes"
                    if (size(Second) > 1) then
                      if (ToLower(trim(adjustl(Second(2)))) == "ija") then
                        This%PedNrmIja = .true.
                        write(STDOUT, "(a)") " Write pedigree NRM format: ija"
                      end if
                    else
                      write(STDOUT, "(a)") " Write pedigree NRM format: matrix"
                    end if
                  else
                    write(STDOUT, "(a)") " Calculate pedigree NRM: No"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify Yes/No[,Ija] for PedNrm, i.e., PedNrm, Yes or PedNrm, Yes, Ija"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("pednrmsubsetfile")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "none") then
                    write(STDOUT, "(a)") " Not using pedigree NRM set file"
                  else
                    This%PedNrmSubsetGiven = .true.
                    write(This%PedNrmSubsetFile, *) trim(adjustl(Second(1)))
                    This%PedNrmSubsetFile = adjustl(This%OutputFormat)
                    write(STDOUT, "(2a)") " Using pedigree NRM set file: ", trim(This%PedNrmSubsetFile)
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify a file for PedNrmSubsetFile, i.e., PedNrmSubsetFile, PedNrmSubset.txt"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              ! @todo
              ! n = CountLines(This%SpecFile)
              ! if (n > 25) then
              !   write(STDOUT, "(a)") " BEWARE: Using an old A matrix is an experimental feature"
              !   write(STDOUT, "(a)") " BEWARE: - It requires id of individuals to be numeric and sequential and no unknown parents"
              !   write(STDOUT, "(a)") " BEWARE: - It requires the old A matrix between the parents of individuals whose A matrix will be built"
              !   write(STDOUT, "(a)") " BEWARE: - It switches off creation of other matrices (exit after AMat is done)"
              !   write(STDOUT, "(a)") " "
              !   read(SpecUnit, *) DumC, This%OldAMatFile, This%OldAMatNInd
              !   This%OldPedNrmFile = .true.
              ! end if

              case ("pednrminv")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "yes") then
                    This%PedNrmInv = .true.
                    write(STDOUT, "(a)") " Calculate pedigree NRM inverse: Yes"
                    if (size(Second) > 1) then
                      if (ToLower(trim(adjustl(Second(2)))) == "ija") then
                        This%PedNrmInvIja = .true.
                        write(STDOUT, "(a)") " Write pedigree NRM inverse format: ija"
                      end if
                    else
                      write(STDOUT, "(a)") " Write pedigree NRM inverse format: matrix"
                    end if
                  else
                    write(STDOUT, "(a)") " Calculate pedigree NRM inverse: No"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify Yes/No[,Format] for PedNrmInv, i.e., PedNrmInv, Yes or PedNrmInv, Yes, Ija"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("gennrm")
                if (allocated(Second)) then
                  if (ToLower(trim(adjustl(Second(1)))) == "yes") then
                    This%GenNrm = .true.
                    write(STDOUT, "(a)") " Calculate genotype NRM: Yes"
                    if (size(Second) > 1) then
                      if (ToLower(trim(adjustl(Second(2)))) == "ija") then
                        This%GenNrmIja = .true.
                        write(STDOUT, "(a)") " Write genotype NRM format: ija"
                      end if
                    else
                      write(STDOUT, "(a)") " Write genotype NRM format: matrix"
                    end if
                  else
                    write(STDOUT, "(a)") " Calculate genotype NRM: No"
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify Yes/No[,Ija] for GenNrm, i.e., GenNrm, Yes or GenNrm, Yes, Ija"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("gennrmtype")
                if (allocated(Second)) then
                  write(This%GenNrmType, *) ToLower(trim(adjustl(Second(1))))
                  This%GenNrmType = adjustl(This%GenNrmType)
                  if (trim(This%GenNrmType) /= "vanraden"        .and. &
                      trim(This%GenNrmType) /= "vanraden1"       .and. &
                      trim(This%GenNrmType) /= "vanraden2"       .and. &
                      trim(This%GenNrmType) /= "yang"            .and. &
                      trim(This%GenNrmType) /= "nejati-javaremi") then
                      ! trim(This%GenNrmType) /= "day-williams") then
                    write(STDERR, "(a)")  " ERROR: GenNrmType must be either VanRaden=VanRaden1, VanRaden2, Yang, or Nejati-Javaremi"
                    write(STDERR, "(a)") " ERROR: |"//trim(This%GenNrmType)//"|"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if
                  if (trim(This%GenNrmType) == "vanraden") then
                    This%GenNrmType = "vanraden1"
                  end if
                  write(STDOUT, "(2a)") " Genotype NRM type: ", trim(This%GenNrmType)
                else
                  write(STDERR, "(a)") " ERROR: Must specify a method for GenNrmType, i.e., GenNrmType, VanRaden"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("fudgegennrmdiag")
                if (allocated(Second)) then
                  This%FudgeGenNrmDiag = .true.
                  This%FudgeGenNrmDiagValue = Char2Double(trim(adjustl(Second(1))), "(f20.16)")
                  write(STDOUT, "(2a)") " Fudge genotype NRM diagonal: ", Real2Char(This%FudgeGenNrmDiagValue, "(f6.4)")
                else
                  write(STDERR, "(a)") " ERROR: Must specify a value for FudgeGenNrmDiag, i.e., FudgeGenNrmDiag, 0.001"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              case ("blendgennrmwithpednrm")
                if (allocated(Second)) then
                  if (size(Second) == 2) then
                    This%PedNrm = .true.
                    This%BlendGenNrmWithPedNrm = .true.
                    This%BlendGenNrmWithPedNrmFactor(1) = Char2Double(trim(adjustl(Second(1))), "(f20.16)")
                    This%BlendGenNrmWithPedNrmFactor(2) = Char2Double(trim(adjustl(Second(2))), "(f20.16)")
                    write(STDOUT, "(a, 2f)") " Blend genotype NRM with pedigree NRM: ", This%BlendGenNrmWithPedNrmFactor
                  else
                    write(STDERR, "(a)") " ERROR: Must specify two values for BlendGenNrmWithPedNrm, i.e., BlendGenNrmWithPedNrm, 0.95, 0.05"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if
                else
                  write(STDERR, "(a)") " ERROR: Must specify two values for BlendGenNrmWithPedNrm, i.e., BlendGenNrmWithPedNrm, 0.95, 0.05"
                  write(STDERR, "(a)") ""
                  stop 1
                end if

              ! case ("fudgehapnrmdiag")
              !   if (allocated(Second)) then
              !     This%FudgeHapNrmDiag = .true.
              !     This%FudgeHapNrmDiagValue = Char2Double(trim(adjustl(Second(1))), "(f20.16)")
              !     write(STDOUT, "(2a)") " Fudge haplotype NRM diagonal: ", Real2Char(This%FudgeHapNrmDiagValue, "(f6.4)")
              !   else
              !     write(STDERR, "(a)") " ERROR: Must specify a value for FudgeHapNrmDiag, i.e., FudgeHapNrmDiag, 0.001"
              !     write(STDERR, "(a)") ""
              !     stop 1
              !   end if

              ! case ("blendhapnrmwithpednrm")
              !   if (allocated(Second)) then
              !     if (size(Second) == 2) then
              !       This%PedNrm = .true.
              !       This%BlendHapNrmWithPedNrm = .true.
              !       This%BlendHapNrmWithPedNrmFactor(1) = Char2Double(trim(adjustl(Second(1))), "(f20.16)")
              !       This%BlendHapNrmWithPedNrmFactor(2) = Char2Double(trim(adjustl(Second(2))), "(f20.16)")
              !       write(STDOUT, "(a, 2f)") " Blend haplotype NRM with pedigree NRM: ", This%BlendHapNrmWithPedNrmFactor
              !     else
              !       write(STDERR, "(a)") " ERROR: Must specify two values for BlendHapNrmWithPedNrm, i.e., BlendHapNrmWithPedNrm, 0.95, 0.05"
              !       write(STDERR, "(a)") ""
              !       stop 1
              !     end if
              !   else
              !     write(STDERR, "(a)") " ERROR: Must specify two values for BlendHapNrmWithPedNrm, i.e., BlendHapNrmWithPedNrm, 0.95, 0.05"
              !     write(STDERR, "(a)") ""
              !     stop 1
              !   end if

              ! if (This%HFullMat .or. This%HIJA) then
              !   This%MakeH = .true.
              !   This%MakeG = .true.
              !   This%MakeA = .true.
              ! end if

              ! read(SpecUnit,*) DumC, Option
              ! This%InvHFullMat = trim(adjustl(Option)) == "Yes"

              ! read(SpecUnit,*) DumC, Option
              ! This%InvHIJA = trim(adjustl(Option)) == "Yes"

              ! if (This%InvHFullMat .or. This%InvHIJA) then
              !   This%MakeInvH = .true.
              !   This%MakeG    = .true.
              !   This%MakeA    = .true.
              !   This%MakeInvA = .true.
              ! end if

              case ("stop")
                write(STDOUT, "(3a)") " NOTE: Encountered Stop specification - the rest of specifications will be ignored"
                write(STDOUT, "(a)") " "
                exit

              case default
                write(STDOUT, "(3a)") " NOTE: Specification '", trim(Line), "' was ignored"
                write(STDOUT, "(a)") " "
            end select
          end if
        end do ReadSpec
        close(SpecUnit)

        if ((This%PedInbreeding .or. This%PedNrm .or. This%PedNrmInv .or. This%PedNrmSubsetGiven)&
            .and. .not. This%PedigreeGiven) then
          write(STDERR, "(a)") " ERROR: Must provide pedigree file to calculate pedigree inbreeding, NRM, or NRM inverse"
          write(STDERR, "(a)") ""
          stop 1
        end if

        if ((This%GenInbreeding .or. This%GenNrm .or. This%GenNrmInv)&
            .and. .not. This%GenotypeGiven) then
          write(STDERR, "(a)") " ERROR: Must provide genotype file to calculate genotype inbreeding, NRM, or NRM inverse"
          write(STDERR, "(a)") ""
          stop 1
        end if

        ! if ((This%HapInbreeding .or. This%HapNrm .or. This%HapNrmInv)&
        !     .and. .not. This%HaplotypeGiven) then
        !   write(STDERR, "(a)") " ERROR: Must provide haplotype file to calculate haplotype inbreeding, NRM, or NRM inverse"
        !   write(STDERR, "(a)") ""
        !   stop 1
        ! end if

        if (This%BlendGenNrmWithPedNrm .and. .not. This%PedigreeGiven) then
          write(STDERR, "(a)") " ERROR: Must provide pedigree file to blend genotype NRM with pedigree NRM"
          write(STDERR, "(a)") ""
          stop 1
        end if

        ! if (This%BlendHapNrm .and. .not. This%PedigreeGiven) then
        !   write(STDERR, "(a)") " ERROR: Must provide pedigree file to blend haplotype NRM with pedigree NRM"
        !   write(STDERR, "(a)") ""
        !   stop 1
        ! end if

        ! if ((This%MakeG .or. This%MakeInvG .or. This%MakeH .or. This%MakeInvH) .and. .not. This%GenotypeGiven) then
        !   write(STDOUT, "(a)") " NOTE: To create G or H matrix, a genotype file must be given --> ommited G or H."
        !   write(STDOUT, "(a)") " "
        !   This%MakeG    = .false.
        !   This%MakeInvG = .false.
        !   This%MakeH    = .false.
        !   This%MakeInvH = .false.
        ! end if

        ! if ((This%MakeA .or. This%MakeInvA .or. This%MakeH .or. This%MakeInvH) .and. .not. This%PedigreeGiven) then
        !   write(STDOUT, "(a)") " NOTE: To create A or H matrix, a pedigree file must be given --> ommited A or H."
        !   write(STDOUT, "(a)") " "
        !   This%MakeA    = .false.
        !   This%MakeInvA = .false.
        !   This%MakeH    = .false.
        !   This%MakeInvH = .false.
        ! end if
      end subroutine

    !###########################################################################

    ! AlphaRelateData type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Read AlphaRelateData from files
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine ReadAlphaRelateData(This, Spec) ! not pure due to IO
        implicit none

        ! Arguments
        class(AlphaRelateData), intent(out) :: This !< @return AlphaRelateData holder
        type(AlphaRelateSpec), intent(in) :: Spec   !< Specifications

        ! Other
        type(PedigreeHolder) :: PedObj

        if (Spec%PedigreeGiven) then
          ! Read in the pedigree
          PedObj = PedigreeHolder(Spec%PedigreeFile)

          ! Sort and recode pedigree
          call PedObj%MakeRecodedPedigreeArray(RecPed=This%RecPed)
          write(STDOUT, "(a1, i8, a)") " ", This%RecPed%nInd," individuals in pedigree"
          call This%RecPed%Write(File=trim(Spec%PedigreeFile)//"_Recoded.txt")

          ! Free some memory
          call PedObj%DestroyPedigree

          ! Handle subset
          if (Spec%PedNrmSubsetGiven) then
            call This%PedNrmSubset%Read(File=Spec%PedNrmSubsetFile)
            write(STDOUT, "(a1, i8, a)") " ", This%PedNrmSubset%nInd," individuals in the pedigree NRM subset"
            call This%PedNrmSubset%MatchId(OriginalIdSuperset=This%RecPed%OriginalId, Skip=1) ! skip=1 because of the "0th margin" in This%RecPed%OriginalId
            block ! @todo make this block a subroutine - it is 99% copied bellow - perhaps Pedigree&Individual types handle this much better?
              integer(int32) :: Ind
              logical :: IdMatchNotFound
              IdMatchNotFound = .false.
              do Ind = 1, This%PedNrmSubset%nInd
                if (This%PedNrmSubset%Id(Ind) == 0) then
                  write(STDERR, "(2a)") " ERROR: No match found in pedigree for a pedigree NRM subset identification: ", trim(This%PedNrmSubset%OriginalId(Ind))
                  IdMatchNotFound = .true.
                end if
              end do
              if (IdMatchNotFound) then
                write(STDERR, "(a)")  ""
                stop 1
              end if
            end block
          end if
        end if

        if (Spec%GenotypeGiven) then

          call This%Gen%Read(File=Spec%GenotypeFile, nLoc=Spec%nLoc)
          write(STDOUT, "(a1, i8, a)") " ", This%Gen%nInd," individuals with genotypes"

          if (Spec%PedigreeGiven) then
            call This%Gen%MatchId(OriginalIdSuperset=This%RecPed%OriginalId, Skip=1) ! skip=1 because of the "0th margin" in This%RecPed%OriginalId
            block ! @todo make this block a subroutine - it is 99% copied bellow - perhaps Pedigree&Individual types handle this much better?
              integer(int32) :: Ind
              logical :: IdMatchNotFound
              IdMatchNotFound = .false.
              do Ind = 1, This%Gen%nInd
                if (This%Gen%Id(Ind) == 0) then
                  write(STDERR, "(2a)") " ERROR: No match found in pedigree for a genotype identification: ", trim(This%Gen%OriginalId(Ind))
                  IdMatchNotFound = .true.
                end if
              end do
              if (IdMatchNotFound) then
                write(STDERR, "(a)")  ""
                stop 1
              end if
            end block
          end if

        end if

        if (Spec%GenotypeGiven .or. Spec%HaplotypeGiven) then
          if (Spec%AlleleFreqGiven) then
            if (Spec%AlleleFreqFixed) then
              call This%AlleleFreq%Init(nLoc=This%Gen%nLoc)
              This%AlleleFreq%Value = Spec%AlleleFreqFixedValue
            else
              call This%AlleleFreq%Read(File=trim(Spec%AlleleFreqFile), nLoc=This%Gen%nLoc)
            end if
          end if
          if (Spec%LocusWeightGiven) then
            call This%LocusWeight%Read(File=trim(Spec%LocusWeightFile), nLoc=This%Gen%nLoc)
          end if
        end if

        ! if (Spec%PedigreeGiven .and. Spec%GenotypeGiven) then
        !   ! These three vectors use the Pedigree animals as base,
        !   ! i.e. after reordering, the index for the nth pedigree animal is n.
        !   allocate(This%MapAnimal(1:(This%nAnisP + This%nAnisG)))
        !   allocate(This%MapToG(1:(This%nAnisP + This%nAnisG)))
        !   allocate(This%AnimalsInBoth(1:This%nAnisP + This%nAnisG)) !@todo, should this be 1:(This%nAnisP + This%nAnisG)?
        !   This%MapAnimal = 0
        !   This%MapToG = .false.
        !   This%AnimalsInBoth = .false.
        !   This%nAnisH = This%nAnisP
        !   do i = 1, This%nAnisP
        !     This%MapAnimal(i) = i
        !   end do
        !
        !   This%AnimalsInBoth = .false.
        !   ! Match genotyped individuals to pedigree
        !   do i = 1, This%nAnisG
        !     GenoInPed = 0
        !     do j = 1, This%nAnisP
        !       if (trim(This%IdGeno(i)) == trim(Id(j))) then ! @todo: can I include Id() into the Data object?
        !         This%MapToG(j) = .true.
        !         This%MapAnimal(j) = i
        !         This%AnimalsInBoth(j) = .true.
        !         GenoInPed = 1
        !         exit
        !       end if
        !     end do
        !     if (GenoInPed == 0) then
        !       This%nAnisH = This%nAnisH + 1
        !       This%MapAnimal(This%nAnisH) = i
        !       This%MapToG(This%nAnisH) = .true.
        !       write(STDOUT, "(2a)") " Genotyped individual not in the pedigree file: ", trim(This%IdGeno(i))
        !       write(STDOUT, "(a)")  " "
        !       ! stop 1
        !     end if
        !   end do
        ! end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write AlphaRelateData to files or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine WriteAlphaRelateData(This, Basename) ! not pure due to IO
        implicit none
        class(AlphaRelateData), intent(in) :: This         !< AlphaRelateData holder
        character(len=*), intent(in), optional :: Basename !< Basename for produced files

        ! The best is to be in line with the AlphaRelateData type definition

        ! Pedigree stuff

        if (allocated(This%RecPed%OriginalId)) then
          if (present(Basename)) then
            call This%RecPed%Write(File=trim(Basename)//"RecodedPedigree.txt")
          else
            write(STDOUT, "(a)") "Recoded pedigree"
            call This%RecPed%Write
          end if
        end if

        if (allocated(This%PedNrmSubset%OriginalId)) then
          if (present(Basename)) then
            call This%PedNrmSubset%Write(File=trim(Basename)//"PedNrmSubset.txt")
          else
            write(STDOUT, "(a)") "Pedigree NRM subset"
            call This%PedNrmSubset%Write
          end if
        end if

        if (allocated(This%PedInbreeding%OriginalId)) then
          if (present(Basename)) then
            call This%PedInbreeding%Write(File=trim(Basename)//"PedInbreeding.txt")
          else
            write(STDOUT, "(a)") "Pedigree inbreeding"
            call This%PedInbreeding%Write
          end if
        end if

        if (allocated(This%PedNrm%OriginalId)) then
          if (present(Basename)) then
            call This%PedNrm%Write(File=trim(Basename)//"PedNrm.txt")
          else
            write(STDOUT, "(a)") "Pedigree NRM"
            call This%PedNrm%Write
          end if
        end if

        if (allocated(This%PedNrmInv%OriginalId)) then
          if (present(Basename)) then
            call This%PedNrmInv%Write(File=trim(Basename)//"PedNrmInv.txt")
          else
            write(STDOUT, "(a)") "Pedigree NRM inverse"
            call This%PedNrmInv%Write
          end if
        end if

        ! Genome stuff

        if (allocated(This%AlleleFreq%Value)) then
          if (present(Basename)) then
            call This%AlleleFreq%Write(File=trim(Basename)//"AlleleFreq.txt")
          else
            write(STDOUT, "(a)") "Allele frequencies"
            call This%AlleleFreq%Write
          end if
        end if

        if (allocated(This%LocusWeight%Value)) then
          if (present(Basename)) then
            call This%LocusWeight%Write(File=trim(Basename)//"LocusWeight.txt")
          else
            write(STDOUT, "(a)") "Locus weights"
            call This%LocusWeight%Write
          end if
        end if

        ! Genotype stuff

        if (allocated(This%Gen%OriginalId)) then
          if (present(Basename)) then
            call This%Gen%Write(File=trim(Basename)//"Genotype.txt")
          else
            write(STDOUT, "(a)") "Genotype"
            call This%Gen%Write
          end if
        end if

        if (allocated(This%Gen%GenotypeReal)) then
          if (present(Basename)) then
            call This%Gen%WriteReal(File=trim(Basename)//"GenotypeReal.txt")
          else
            write(STDOUT, "(a)") "Genotype (as real)"
            call This%Gen%WriteReal
          end if
        end if

        if (allocated(This%GenInbreeding%OriginalId)) then
          if (present(Basename)) then
            call This%GenInbreeding%Write(File=trim(Basename)//"GenInbreeding.txt")
          else
            write(STDOUT, "(a)") "Genotype inbreeding"
            call This%GenInbreeding%Write
          end if
        end if

        if (allocated(This%GenNrm%OriginalId)) then
          if (present(Basename)) then
            call This%PedNrm%Write(File=trim(Basename)//"GenNrm.txt")
          else
            write(STDOUT, "(a)") "Genotype NRM"
            call This%GenNrm%Write
          end if
        end if

        if (allocated(This%GenNrmInv%OriginalId)) then
          if (present(Basename)) then
            call This%GenNrmInv%Write(File=trim(Basename)//"GenNrmInv.txt")
          else
            write(STDOUT, "(a)") "Genotype NRM inverse"
            call This%GenNrmInv%Write
          end if
        end if

        ! Haplotype stuff

        ! @todo Write Haplotypes
        ! if (allocated(This%Hap%OriginalId)) then
        !   if (present(Basename)) then
        !     call This%Hap%Write(File=trim(Basename)//"Haplotype.txt")
        !   else
        !     write(STDOUT, "(a)") "Haplotype"
        !     call This%Hap%Write
        !   end if
        ! end if

        ! @todo Write Haplotypes results
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  AlphaRelateData destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine DestroyAlphaRelateData(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder

        ! The best is to be in line with the AlphaRelateData type definition

        ! Pedigree stuff

        if (allocated(This%RecPed%OriginalId)) then
          call This%RecPed%Destroy
        end if

        if (allocated(This%PedInbreeding%OriginalId)) then
          call This%PedInbreeding%Destroy
        end if

        if (allocated(This%PedNrm%OriginalId)) then
          call This%PedNrm%Destroy
        end if

        if (allocated(This%PedNrmSubset%OriginalId)) then
          call This%PedNrmSubset%Destroy
        end if

        if (allocated(This%PedNrmInv%OriginalId)) then
          call This%PedNrmInv%Destroy
        end if

        ! Genome stuff

        if (allocated(This%AlleleFreq%Value)) then
          call This%AlleleFreq%Destroy
        end if

        if (allocated(This%LocusWeight%Value)) then
          call This%LocusWeight%Destroy
        end if

        ! Genotype stuff

        if (allocated(This%Gen%OriginalId)) then
          call This%Gen%Destroy
        end if

        if (allocated(This%GenInbreeding%OriginalId)) then
          call This%GenInbreeding%Destroy
        end if

        if (allocated(This%GenNrm%OriginalId)) then
          call This%GenNrm%Destroy
        end if

        if (allocated(This%GenNrmInv%OriginalId)) then
          call This%GenNrmInv%Destroy
        end if

        ! Haplotype stuff

        ! if (allocated(This%Hap%OriginalId)) then
        !   call This%Hap%Destroy
        ! end if

        ! @todo
        ! type(Inbreeding)           :: HapInbreeding
        ! type(Nrm)                  :: HapNrm
        ! type(Nrm)                  :: HapNrmInv

      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree inbreeding on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcPedInbreedingAlphaRelateData(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder, note This%PedInbreeding(0) = -1.0!!!
        call This%PedInbreeding%Init(nInd=This%RecPed%nInd)
        This%PedInbreeding%OriginalId = This%RecPed%OriginalId
        This%PedInbreeding%Inb = PedInbreedingMeuwissenLuo(RecPed=This%RecPed%Id, n=This%PedInbreeding%nInd)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcPedNrmAlphaRelateData(This, Spec)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder
        type(AlphaRelateSpec), intent(in) :: Spec     !< Specifications

        if (Spec%PedNrmSubsetGiven) then
          ! Using the Colleau/Tier method to get Nrm for a subset of individuals
          call This%PedNrm%Init(nInd=This%PedNrmSubset%nInd)
          This%PedNrm%OriginalId = This%PedNrmSubset%OriginalId
          This%PedNrm%Id = This%PedNrmSubset%Id
          block
            integer(int32) :: Ind, xPos
            real(real64), allocatable, dimension(:) :: x, NrmCol
            allocate(     x(0:This%RecPed%nInd))
            allocate(NrmCol(0:This%RecPed%nInd))
            call This%CalcPedInbreeding
            x = 0.0d0
            ! @todo: this could be run in parallel (is it worth it?; x must be made private!!!)
            do Ind = 1, This%PedNrm%nInd
              xPos = This%PedNrmSubset%Id(Ind)
              x(xPos) = 1.0d0
              NrmCol = PedNrmTimesVector(RecPed=This%RecPed%Id, n=This%RecPed%nInd,&
                                         Inbreeding=This%PedInbreeding%Inb, Vector=x)
              This%PedNrm%Nrm(0:This%PedNrm%nInd, Ind) = NrmCol(This%PedNrmSubset%Id)
              x(xPos) = 0.0d0
            end do
            deallocate(NrmCol)
            deallocate(x)
          end block
        else if (Spec%PedNrmOldGiven) then
          ! @todo: this needs work
          ! @todo: Put this into a block
          ! type(Nrm) :: OldNrm
          ! integer(int32) :: Ind, MinOldId, MaxOldId
          ! logical :: OldIdUnknown
          ! @todo: read this already in the Data function!!!
          ! call ReadNrm(File=Spec%OldPedNrmFile, Ija=Spec%PedNrmIja,&
          !              OriginalId=OldNrm%OriginalId, Nrm=OldNrm%Nrm, n=OldNrm%nInd)
          ! MinOldId = 1
          ! MaxOldId = 1
          ! OldIdUnknown = .true.
          ! Ind = 0
          ! do while (OldIdUnknown)
          !   Ind = Ind + 1
          !   if (OldNrm%OriginalId(1)           == This%RecPed%OriginalId(Ind)) then
          !     MinOldId = Ind
          !     OldIdUnknown = .false.
          !   end if
          ! end do
          ! OldIdUnknown = .true.
          ! Ind = 0
          ! do while (OldIdUnknown)
          !   Ind = Ind + 1
          !   if (OldNrm%OriginalId(OldNrm%nInd) == This%RecPed%OriginalId(Ind)) then
          !     MaxOldId = Ind
          !     OldIdUnknown = .false.
          !   end if
          ! end do
          !
          ! This%PedNrm%nInd = This%RecPed%nInd - MaxOldId
          !
          ! if (allocated(This%PedNrm%OriginalId)) then
          !   deallocate(This%PedNrm%OriginalId)
          ! end if
          ! allocate(This%PedNrm%OriginalId(0:This%PedNrm%nInd))
          ! This%PedNrm%OriginalId(0) = EMPTYID
          ! This%PedNrm%OriginalId(1:This%PedNrm%nInd) = This%RecPed%OriginalId((MaxOldId + 1):This%RecPed%nInd)
          !
          ! if (allocated(This%PedNrm%Nrm)) then
          !   deallocate(This%PedNrm%Nrm)
          ! end if
          ! allocate(This%PedNrm%Nrm(0:This%PedNrm%nInd, 0:This%PedNrm%nInd))
          !
          ! This%PedNrm%Nrm = PedNrmWithOldNrm(RecPed=This%RecPed%Id, n=This%RecPed%nInd,&
          !                                      nNew=This%PedNrm%nInd,&
          !                                      OldNrm=OldNrm%Nrm, nOld=OldNrm%nInd,&
          !                                      MinOldId=MinOldId, MaxOldId=MaxOldId)
        else
          ! Standard method for all individuals
          call This%PedNrm%Init(nInd=This%RecPed%nInd)
          This%PedNrm%OriginalId = This%RecPed%OriginalId
          This%PedNrm%Nrm = PedNrm(RecPed=This%RecPed%Id, n=This%PedNrm%nInd)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM inverse on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcPedNrmInvAlphaRelateData(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder
        if (.not. allocated(This%PedInbreeding%Inb)) then
          call This%CalcPedInbreeding
        end if
        call This%PedNrmInv%Init(nInd=This%RecPed%nInd)
        This%PedNrmInv%OriginalId = This%RecPed%OriginalId
        This%PedNrmInv%Nrm = PedNrmInv(RecPed=This%RecPed%Id, n=This%PedNrmInv%nInd,&
                                       Inbreeding=This%PedInbreeding%Inb)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate allele frequencies on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcAlleleFreqAlphaRelateData(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder

        integer(int32) :: Ind, Loc
        integer(int32), allocatable, dimension(:) :: nObs

        call This%AlleleFreq%Init(nLoc=This%Gen%nLoc)
        ! @todo can we build a method that works on type(Genotype) so we would not need GenotypeReal for allele freq calculation?
        if (.not. allocated(This%Gen%GenotypeReal)) then
          call This%Gen%MakeGenotypeReal
        end if

        ! @todo could we not assume that GenotypeReal has no missing values
        !       (it should have been cleaned prior to this program) and then
        !       we can avoid the ifs and simplify computatation a lot
        ! do Ind = 1, This%Gen%nInd
        !   This%AlleleFreq%Value = This%AlleleFreq%Value + This%Gen%GenotypeReal(:, Ind)
        ! end do
        ! This%AlleleFreq%Value = This%AlleleFreq%Value / (2.0d0 * dble(This%Gen%nInd))
        allocate(nObs(This%Gen%nLoc))
        nObs = 0
        do Ind = 1, This%Gen%nInd
          do Loc = 1, This%Gen%nLoc
            if ((This%Gen%GenotypeReal(Loc, Ind) .ge. 0.0d0) .and. (This%Gen%GenotypeReal(Loc, Ind) .le. 2.0d0)) then
              This%AlleleFreq%Value(Loc) = This%AlleleFreq%Value(Loc) + This%Gen%GenotypeReal(Loc, Ind)
              nObs(Loc) = nObs(Loc) + 1
            end if
          end do
        end do
        do Loc = 1, This%Gen%nLoc
          if (nObs(Loc) .gt. 0) then
            This%AlleleFreq%Value(Loc) = This%AlleleFreq%Value(Loc) / (2.0d0 * dble(nObs(Loc)))
          else
            This%AlleleFreq%Value(Loc) = 0.0d0
          end if
        end do
        deallocate(nObs)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate genotype NRM on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcGenNrmAlphaRelateData(This, Spec)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder
        type(AlphaRelateSpec), intent(in) :: Spec     !< Specifications

        call This%GenNrm%Init(nInd=This%Gen%nInd)
        This%GenNrm%OriginalId = This%Gen%OriginalId
        This%GenNrm%Id = This%Gen%Id

        ! GEMM docs https://software.intel.com/en-us/node/468480
        ! NOTE: Cov(a|Z) = Cov(Zalpha|Z) = ZVar(alpha)Z', where Z(nInd, nLoc)
        ! NOTE: Z here is (nLoc, 0:nInd), hence need to compute Z'Z

        select case (trim(Spec%GenNrmType))
          case ("vanraden1")
            ! Setup
            if (.not. allocated(This%Gen%GenotypeReal)) then
              call This%Gen%MakeGenotypeReal
            end if
            if (.not. Spec%AlleleFreqGiven) then
              call This%CalcAlleleFreq
            end if
            ! Center
            call This%Gen%CenterGenotypeReal(AlleleFreq=This%AlleleFreq%Value)
            ! Weight
            if (Spec%LocusWeightGiven) then
              call This%Gen%WeightGenotypeReal(Weight=sqrt(This%LocusWeight%Value)) ! sqrt, because we do (Zsqrt(W))'(sqrt(W)Z) later
            end if
            ! Compute ~Z'Z
            call gemm(A=This%Gen%GenotypeReal, B=This%Gen%GenotypeReal, C=This%GenNrm%Nrm, TransA="T")
            ! Average over loci
            This%GenNrm%Nrm = This%GenNrm%Nrm / (2.0d0 * sum(This%AlleleFreq%Value * (1.0d0 - This%AlleleFreq%Value)))
            ! This%GenNrm%Nrm = GenNrmVanRaden1LoopOnGenotype(Genotype=This%Gen%Genotype,&
            !                                                 nInd=This%Gen%nInd,&
            !                                                 nLoc=This%Gen%nLoc,&
            !                                                 AlleleFreq=nLoc=This%AlleleFreq%Value)

          case ("vanraden2")
            ! Setup
            if (.not. allocated(This%Gen%GenotypeReal)) then
              call This%Gen%MakeGenotypeReal
            end if
            if (.not. Spec%AlleleFreqGiven) then
              call This%CalcAlleleFreq
            end if
            ! Center and scale
            call This%Gen%CenterAndScaleGenotypeReal(AlleleFreq=This%AlleleFreq%Value)
            ! Weight
            if (Spec%LocusWeightGiven) then
              call This%Gen%WeightGenotypeReal(Weight=sqrt(This%LocusWeight%Value)) ! sqrt, because we do (Zsqrt(W))'(sqrt(W)Z) later
            end if
            ! Compute ~Z'Z
            call gemm(A=This%Gen%GenotypeReal, B=This%Gen%GenotypeReal, C=This%GenNrm%Nrm, TransA="T")
            ! Average over loci
            This%GenNrm%Nrm = This%GenNrm%Nrm / dble(This%Gen%nLoc)

          case ("yang")
            block
              integer(int32) :: Ind, Loc
              real(real64) :: Diag(This%Gen%nInd), T2AFLoc(This%Gen%nLoc), Het(This%Gen%nLoc), Tmp
              ! Setup
              if (.not. allocated(This%Gen%GenotypeReal)) then
                call This%Gen%MakeGenotypeReal
              end if
              if (.not. Spec%AlleleFreqGiven) then
                call This%CalcAlleleFreq
              end if
              ! Prepare diagonal (need to do it before GenotypeReal gets centered and scaled!!!)
              Diag = 0.0d0
              T2AFLoc = 2.0d0 * This%AlleleFreq%Value
              Het = T2AFLoc * (1.0d0 - This%AlleleFreq%Value)
              do Ind = 1, This%GenNrm%nInd
                do Loc = 1, This%Gen%nLoc
                  if (Het(Loc) > tiny(Het(Loc))) then ! to avoid dividing by zero
                    if (Spec%LocusWeightGiven) then
                      Tmp = This%LocusWeight%Value(Loc)
                    else
                      Tmp = 1.0d0
                    end if
                    Diag(Ind) = Diag(Ind) + &
                      Tmp * (1.0d0 + &
                      ((This%Gen%GenotypeReal(Loc, Ind) * This%Gen%GenotypeReal(Loc, Ind)) - &
                       ((1.0d0 + T2AFLoc(Loc)) * This%Gen%GenotypeReal(Loc, Ind)) + &
                       (T2AFLoc(Loc) * This%AlleleFreq%Value(Loc))) / Het(Loc))
                  end if
                end do
              end do
              ! Center and scale
              call This%Gen%CenterAndScaleGenotypeReal(AlleleFreq=This%AlleleFreq%Value)
              ! Weight
              if (Spec%LocusWeightGiven) then
                call This%Gen%WeightGenotypeReal(Weight=sqrt(This%LocusWeight%Value)) ! sqrt, because we do (Zsqrt(W))'(sqrt(W)Z) later
              end if
              ! Compute ~Z'Z
              call gemm(A=This%Gen%GenotypeReal, B=This%Gen%GenotypeReal, C=This%GenNrm%Nrm, TransA="T")
              ! Modify diagonal
              do Ind = 1, This%GenNrm%nInd
                This%GenNrm%Nrm(Ind, Ind) = Diag(Ind)
              end do
              ! Average over loci
              This%GenNrm%Nrm = This%GenNrm%Nrm / dble(This%Gen%nLoc)
            end block

          case ("nejati-javaremi")
            block
              real(real64) :: AlleleFreqHalf(This%Gen%nLoc), Tmp
              ! Setup
              if (.not. allocated(This%Gen%GenotypeReal)) then
                call This%Gen%MakeGenotypeReal
              end if
              ! Center with allele freq of 0.5
              AlleleFreqHalf = 0.5d0
              call This%Gen%CenterGenotypeReal(AlleleFreq=AlleleFreqHalf)
              ! Weight
              if (Spec%LocusWeightGiven) then
                call This%Gen%WeightGenotypeReal(Weight=sqrt(This%LocusWeight%Value)) ! sqrt, because we do (Zsqrt(W))'(sqrt(W)Z) later
              end if
              ! Compute ~Z'Z
              call gemm(A=This%Gen%GenotypeReal, B=This%Gen%GenotypeReal, C=This%GenNrm%Nrm, TransA="T")
              ! Average over loci
              This%GenNrm%Nrm = This%GenNrm%Nrm / dble(This%Gen%nLoc)
              ! Modify scale from [-1,1] to [0,2]
              if (Spec%LocusWeightGiven) then
                Tmp = sum(This%LocusWeight%Value) / dble(This%Gen%nLoc)
              else
                Tmp = 1.0d0
              end if
              This%GenNrm%Nrm = This%GenNrm%Nrm + Tmp
              ! Make sure the "0th" margin is 0.0
              This%GenNrm%Nrm(0:This%GenNrm%nInd, 0) = 0.0d0
              This%GenNrm%Nrm(0, 0:This%GenNrm%nInd) = 0.0d0
            end block
        end select

        if (Spec%FudgeGenNrmDiag) then
          block
            integer(int32) :: Ind
            do Ind = 1, This%GenNrm%nInd
              This%GenNrm%Nrm(Ind, Ind) = This%GenNrm%Nrm(Ind, Ind) + Spec%FudgeGenNrmDiagValue
            end do
          end block
        end if

        if (Spec%BlendGenNrmWithPedNrm) then
          if (.not. allocated(This%PedNrm%Nrm)) then
            call This%CalcPedNrm(Spec=Spec)
          end if
          This%GenNrm%Nrm = Spec%BlendGenNrmWithPedNrmFactor(1) * This%GenNrm%Nrm + &
                            Spec%BlendGenNrmWithPedNrmFactor(2) * This%PedNrm%Nrm(This%GenNrm%Id, This%GenNrm%Id)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate genotype NRM - VanRaden1 loop on genotype type
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !> @details This function is not tested to work properly!
      !-------------------------------------------------------------------------
      pure function GenNrmVanRaden1LoopOnGenotype(GenotypeInput, nInd, nLoc, AlleleFreq) result(Nrm)
        implicit none

        ! Arguments
        type(Genotype), intent(in) :: GenotypeInput(0:nInd) !< Genotypes
        integer(int32), intent(in) :: nInd                  !< Number of individuals
        integer(int32), intent(in) :: nLoc                  !< Number of loci
        integer(int32), intent(in) :: AlleleFreq(nLoc)      !< Allele frequencies
        real(real64) :: Nrm(0:nInd, 0:nInd)                 !< @return Genotype NRM

        integer(int32) :: Ind1, Ind2
        real(real64) :: Scale, Genotype2(nLoc)

        Scale = 2.0d0 * sum(AlleleFreq * (1.0d0 - AlleleFreq))
        Nrm(0:nInd, 0) = 0.0d0
        Nrm(0, 0:nInd) = 0.0d0
        do Ind2 = 1, nInd
          Genotype2 = dble(GenotypeInput(Ind2)%ToIntegerArray()) - AlleleFreq
          Nrm(Ind2, Ind2) = dot(Genotype2, Genotype2) / Scale
          do Ind1 = Ind2 + 1, nInd
            Nrm(Ind1, Ind2) = dot(dble(GenotypeInput(Ind1)%ToIntegerArray()) - AlleleFreq, Genotype2) / Scale
            Nrm(Ind2, Ind1) = Nrm(Ind1, Ind2)
          end do
        end do
      end function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate genotype inbreeding on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcGenInbreedingAlphaRelateData(This, Spec)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder, note This%GenInbreeding(0) = -1.0!!!
        class(AlphaRelateSpec), intent(in) :: Spec    !< AlphaRelateSpecs

        integer(int32) :: Ind

        if (.not. allocated(This%GenNrm%Nrm)) then
          call This%CalcGenNrm(Spec=Spec)
        end if

        call This%GenInbreeding%Init(nInd=This%GenNrm%nInd)
        This%GenInbreeding%OriginalId = This%GenNrm%OriginalId
        This%GenInbreeding%Id = This%GenNrm%Id
        do Ind = 1, This%GenInbreeding%nInd
          This%GenInbreeding%Inb(Ind) = This%GenNrm%Nrm(Ind, Ind) - 1.0d0
        end do
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate genotype NRM inverse on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 9, 2016
      !-------------------------------------------------------------------------
!pure subroutine CalcGenNrmInvAlphaRelateData(This, Spec, Info)
subroutine CalcGenNrmInvAlphaRelateData(This, Spec, Info)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder
        class(AlphaRelateSpec), intent(in) :: Spec    !< AlphaRelateSpecs
        logical, intent(out) :: Info                  !< @return Success of inversion (.true.) or failure (.false.)

        integer(int32) :: Ind
        integer :: InfoInt

        if (.not. allocated(This%GenNrm%Nrm)) then
          call This%CalcGenNrm(Spec=Spec)
        end if
        call This%GenNrmInv%Init(nInd=This%GenNrm%nInd)
        This%GenNrmInv%OriginalId = This%GenNrm%OriginalId
        This%GenNrmInv%Id = This%GenNrm%Id

        This%GenNrmInv%Nrm = This%GenNrm%Nrm

        Info = .true.
        InfoInt = 0

        ! Cholesky factorization of a symmetric positive definite matrix
        ! https://software.intel.com/en-us/node/468690
call This%GenNrmInv%Write
        ! call potrf(A=This%GenNrmInv%Nrm, Info=InfoInt)
        call potrf(A=This%GenNrmInv%Nrm)
call This%GenNrmInv%Write
        if (InfoInt == 0) then
          ! Computes the inverse based on the Cholesky factor obtained with potrf()
          ! https://software.intel.com/en-us/node/468824
          ! call potri(A=This%GenNrmInv%Nrm, Info=InfoInt)
          call potri(A=This%GenNrmInv%Nrm)
call This%GenNrmInv%Write
          if (InfoInt == 0) then
            ! Fill the lower triangle
            ! @todo would not need these memory jump if the Nrm would be considered as symmetric
            do Ind = 1, This%GenNrmInv%nInd
              This%GenNrmInv%Nrm((Ind + 1):This%GenNrmInv%nInd, Ind) = This%GenNrmInv%Nrm(Ind, (Ind + 1):This%GenNrmInv%nInd)
            end do

call This%GenNrmInv%Write
          else
            Info = .false.
          end if
        else
          Info = .false.
        end if
      end subroutine

      !#########################################################################

    !###########################################################################

    ! @todo Old code

      !###########################################################################

      ! subroutine MakeHAndInvHMatrix
      !   ! Feature added by Stefan Hoj-Edwards, February 2016
      !   ! Making the Inverse H matrix ala Aguilar et al 201? and Christensen 2012

      !   ! Prerequisite and assumptions for this subroutine:
      !   ! There is given both a pedigree and genotype file, and there is an overlap
      !   ! of animals between two data sets.
      !   ! Diagonals of from both A have been collected during MakeA and MakeG,
      !   ! as well as average of A22.
      !   ! GMat is already calculated and loaded in memory.
      !   ! InvA is calculated an loaded into memory.
      !   !
      !   ! Further assumes that animals are ordered the same in both A and G.

      !   implicit none

      !   integer(int32) :: i,j,k,m,p,q,div,t1,t2,whichMat,nBoth
      !   integer(int32),allocatable :: MapToA11(:), MapToA22(:) !Gmap(:),

      !   real(real64) :: GMatavg, nom, denom, slope, intercept, Gmean, Amean, Hii
      !   real(real64),allocatable :: Gdiag(:), Hrow(:), A22(:,:), InvA22(:,:), G22(:,:), A11(:,:), A12(:,:), tmp(:,:), Gboth(:,:)

      !   character(len=1000) :: nChar,fmt1, fmt2,filout
      !   character(len=IDLENGTH),allocatable :: Ids(:)

      !   logical,allocatable :: AnimToWrite(:)

      !   nboth = count(AnimalsInBoth)
      !   ! Make H and/or InvH
      !   allocate(Ids(1:nAnisH))
      !   allocate(AnimToWrite(1:nAnisH))

      !   do i=1,nAnisH
      !     if (MapToG(i)) then
      !       Ids(i) = IdGeno(MapAnimal(i))
      !       AnimToWrite(i) = .true.
      !     else
      !       Ids(i) = Id(MapAnimal(i))
      !       AnimToWrite(i) = RecPed(MapAnimal(i),4)
      !     end if
      !   end do

      !   allocate(InvA22(nBoth,nBoth))
      !   allocate(MapToA22(nAnisH))
      !   if (MakeH) then
      !     allocate(A22(nBoth,nBoth))
      !   end if

      !   k = 0
      !   do i=1,nAnisP
      !     if (.not. AnimalsInBoth(i)) then
      !       cycle
      !     end if
      !     k = k + 1
      !     MapToA22(i) = k
      !     m = 0
      !     do j=1,nAnisP
      !       if (.not. AnimalsInBoth(j)) then
      !         cycle
      !       end if
      !       m = m + 1
      !       InvA22(m,k) = AMat(j,i)
      !     end do
      !   end do
      !   if (MakeH) then
      !     A22 = InvA22
      !   end if

      !   call invert(InvA22,size(InvA22,1),.true.,1)

      !   ! This is the G matrix in Legarra,
      !   ! Sadly, no genotypes where provided, instead the resulting G matrix was.
      !   if (.false.) then
      !     print *, "Overwriting G matrix with example in Legarra 2008!"
      !     do i=1,nAnisG
      !       do j=1,nAnisG
      !         if (i==j) then
      !           GMat(i,j,1) = 1
      !         else
      !           GMat(i,j,1) = 0.7
      !         end if
      !       end do
      !     end do
      !   end if

      !   whichMat = 0
      !   do t1=1,nTrait
      !     do t2=t1,nTrait
      !       whichMat = whichMat + 1

      !       write(*, '(" Starting on H matrix "i0" - "i0)') t1, t2

      !       ! Collect G22
      !       allocate(G22(nAnisG,nAnisG))

      !       G22 = 0.0d0
      !       do j=1,nAnisG
      !         do i=1,nAnisG
      !           nom = GMat(i,j,whichMat)
      !           if (i == j) then
      !             nom = nom - DiagFudge
      !           end if
      !           G22(i,j) = nom
      !         end do
      !       end do

      !       if (ScaleGByRegression) then
      !         allocate(Gdiag(0:nBoth))
      !         Gdiag=0.0d0
      !         GMatavg=0.0d0
      !         div=dble(nBoth**2)
      !         !allocate(Gmap(nBoth))

      !         k = 0
      !         do i=1,nAnisH
      !           if (.not. AnimalsInBoth(i)) then
      !             cycle
      !           end if
      !           k = k+1
      !           Gdiag(k) = G22(MapAnimal(i),MapAnimal(i))
      !           do j=1,nAnisH
      !             if (.not. AnimalsInBoth(j)) then
      !               cycle
      !             end if
      !             GMatavg=GMatavg + G22(MapAnimal(j),MapAnimal(i))/div
      !           end do
      !         end do
      !         Gdiag(0) = GMatavg

      !         ! Now do simple linear regression
      !         nom = 0.0d0
      !         denom = 0.0d0
      !         Gmean = sum(Gdiag) / dble(size(Gdiag, 1))
      !         Amean = sum(Adiag) / dble(size(Adiag, 1))
      !         do i=0,ubound(Adiag, 1)
      !           nom = nom + (Adiag(i) - Amean) * (Gdiag(i) - Gmean)
      !           denom = denom + (Adiag(i) - Amean)**2
      !         end do
      !         slope = nom / denom
      !         intercept = Amean - slope * Gmean

      !         ! Scale G
      !         G22 = slope * G22 + intercept
      !         !do i=1,nAnisG
      !         ! G22(i,i) = G22(i,i) + DiagFudge
      !         !end do
      !         print *, "Scaling of G:"
      !         write(*, "(a,f7.4,a,f7.4)"), " G* = G x ", slope, " + ", intercept
      !         deallocate(Gdiag)
      !       else
      !         do i=1,nAnisH
      !           if (.not. MapToG(i)) then
      !             cycle
      !           end if
      !           do j=1,nAnisH
      !             if (.not. MapToG(j)) then
      !               cycle
      !             end if
      !             if (AnimalsInBoth(i) .and. AnimalsInBoth(j)) then
      !               G22(MapAnimal(j),MapAnimal(i)) = ScaleGToA * G22(MapAnimal(j),MapAnimal(i)) + (1.0d0 - ScaleGToA) * AMat(j,i)
      !             end if
      !           end do
      !         end do
      !       end if

      !       do i=1,nAnisG
      !         G22(i,i) = G22(i,i) + DiagFudge
      !       end do

      !       allocate(Hrow(1:count(AnimToWrite)))

      !       if (MakeH) then

      !         allocate(A11(nAnisP-nBoth, nAnisP-nBoth))
      !         allocate(A12(nAnisP-nBoth, nBoth))
      !         allocate(MapToA11(nAnisP))
      !         allocate(tmp(nAnisP-nBoth, nBoth))
      !         allocate(Gboth(nBoth,nBoth))

      !         MapToA11 = 0
      !         k = 0
      !         p = 0
      !         do i=1,nAnisP
      !           if (AnimalsInBoth(i)) then
      !             p = p + 1
      !             q = 0
      !             do j=1,nAnisP
      !               if (.not. AnimalsInBoth(j)) then
      !                 cycle
      !               end if
      !               q = q + 1
      !               Gboth(q,p) = G22(MapAnimal(j),MapAnimal(i))
      !             end do
      !           else
      !             k = k+1
      !             m = 0
      !             MapToA11(i) = k
      !             do j=1,nAnisP
      !               if (AnimalsInBoth(j)) then
      !                 A12(k,MapAnimal(j)) = AMat(j,i)  !A12 is not symmetrical
      !               else
      !                 m = m+1
      !                 A11(m,k) = AMat(j,i)
      !               end if
      !             end do
      !           end if
      !         end do

      !         ! @todo: use DGEMM
      !         tmp = matmul(A12, InvA22)
      !         !tmp = matmul(matmul(tmp, (Gboth - A22)), transpose(tmp))
      !         tmp = matmul(tmp, (Gboth-A22))
      !         tmp = matmul(tmp, InvA22)
      !         !tmp = matmul(tmp, transpose(A12))

      !         A11 = A11 + matmul(tmp, transpose(A12))
      !         A12 = matmul(matmul(A12, InvA22), Gboth)

      !         deallocate(tmp)
      !         deallocate(Gboth)

      !         print *, "Start writing H matrices (full and/or ija)"

      !         if (HFullMat) then
      !           write(filout,'("HFullMatrix"i0,"-"i0".txt")') t1,t2
      !           write(nChar,*) nAnisH
      !           fmt1="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
      !           open(unit=202,file=trim(filout),status="unknown")
      !         end if

      !         if (HIJA) then
      !           write(filout,'("Hija"i0,"-"i0".txt")') t1,t2
      !           fmt2="(a20,a20,"//trim(adjustl(OutputFormat))//")"
      !           open(unit=204,file=trim(filout),status="unknown")
      !         end if

      !         do i=1,nAnisH
      !           if (AnimToWrite(i) .eq. .false.) then
      !             cycle
      !           end if
      !           Hrow = 0
      !           k = 0
      !           do j=1,nAnisH
      !             if (AnimToWrite(j) .eq. .false.) then
      !               cycle
      !             end if
      !             k = k + 1
      !             if (MapToG(i)) then
      !               if (MapToG(j)) then
      !                 Hii = G22(MapAnimal(i),MapAnimal(j))
      !               else
      !                 Hii = A12(MapToA11(j),MapAnimal(i)) ! Remember to transpose
      !               end if
      !             else
      !               if (MapToG(j)) then
      !                 Hii = A12(MapToA11(i),MapAnimal(j))
      !               else
      !                 Hii = A11(MapToA11(i),MapToA11(j))
      !               end if
      !             end if
      !             if (InvHIJA .and. i .le. j .and. Hii /= 0.0d0) then
      !               write(204,fmt2) Ids(i), Ids(j), Hii
      !             end if
      !             Hrow(k) = Hii
      !           end do
      !           if (HFullMat) then
      !             write(202,fmt1) Ids(i),Hrow(:)
      !           end if
      !         end do

      !       end if

      !       if (MakeInvH) then
      !         print *, "Start inverting scaled G matrix"
      !         call invert(G22, size(G22, 1), .true., 1)

      !         do i=1,nAnisH
      !           if (AnimToWrite(i) .eq. .false.) then
      !             cycle
      !           end if
      !           Hrow = 0
      !           k = 0
      !           do j=1,nAnisH
      !             if (AnimToWrite(j) .eq. .false.) then
      !               cycle
      !             end if
      !             k = k + 1
      !             if (MapToG(i) .and. MapToG(j)) then
      !               Hrow(k) = G22(MapAnimal(i),MapAnimal(j))
      !               if (i <= nAnisP .and. j <= nAnisP) then
      !                 Hrow(k) = Hrow(k) + InvAMat(i,j) - InvA22(MapToA22(i),MapToA22(j))
      !               end if
      !             else if (i <= nAnisP .and. j <= nAnisP) then !if (MapToG(i) .eq. .false. .and. MapToG(j) .eq. .false.  ) then
      !               Hrow(k) = InvAMat(i,j)
      !             end if
      !             if (InvHIJA .and. i .le. j .and. Hrow(k) /= 0.0d0) then
      !               write(204,fmt2) trim(Ids(i)), trim(Ids(j)), Hrow(k)
      !             end if
      !           end do
      !           if (InvHFullMat) then
      !             write(202,fmt1) Ids(i),Hrow(:)
      !           end if
      !         end do

      !#########################################################################

    !###########################################################################

    ! Functions

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree inbreeding using the recursive method
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   January 13, 2017
      !-------------------------------------------------------------------------
      pure function PedInbreedingRecursive(RecPed, n) result(f)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3, 0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n                !< Number of individuals in pedigree
        real(real64) :: f(0:n)                         !< @return Pedigree inbreeding, note PedInbreeding(0) = -1.0!!!

        ! Other
        integer(int32) :: Ind

        f(0) = -1.0d0

        do Ind = 1, n
          if (RecPed(2, Ind) .eq. 0 .or. RecPed(3, Ind) .eq. 0) then
            f(Ind) = 0.0d0
          else
            f(Ind) = 0.5d0 * PedNrmRecursive(Ind1=RecPed(2, Ind), Ind2=RecPed(3, Ind))
          end if
        end do

        contains

          !---------------------------------------------------------------------
          !> @brief  Calculate pedigree numerator relationship between two individuals
          !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
          !> @date   January 13, 2017
          !---------------------------------------------------------------------
          pure recursive function PedNrmRecursive(Ind1, Ind2) result(r)
            integer(int32), intent(in) :: Ind1 !< First  individual
            integer(int32), intent(in) :: Ind2 !< Second individual
            real(real64) :: r                  !< @return relationship between the individuals

            ! Needs access to f(0:n) and RecPed(3, 0:n)

            if (Ind1 .eq. 0 .or. Ind2 .eq. 0) then
              r = 0.0d0
            else if (Ind1 .eq. Ind2) then
              r = 1.0d0 + f(Ind1)
            else
              if (Ind1 .lt. Ind2) then
                r = 0.5d0 * (PedNrmRecursive(Ind1=Ind1, Ind2=RecPed(2, Ind2)) + &
                             PedNrmRecursive(Ind1=Ind1, Ind2=RecPed(3, Ind2)))
              else
                r = 0.5d0 * (PedNrmRecursive(Ind1=Ind2, Ind2=RecPed(2, Ind1)) + &
                             PedNrmRecursive(Ind1=Ind2, Ind2=RecPed(3, Ind1)))
              end if
            end if
          end function
      end function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree inbreeding using the Meuwissen and
      !!         Luo (1992, GSE 24: 305-313) method
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk & John Hickey, john.hickey@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure function PedInbreedingMeuwissenLuo(RecPed, n) result(f)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3, 0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n                !< Number of individuals in pedigree
        real(real64) :: f(0:n)                         !< @return Pedigree inbreeding, note PedInbreeding(0) = -1.0!!!

        ! Other
        integer(int32) :: i, is, id, j, k, ks, kd
        integer(int32) :: ped(3, 0:n), point(0:n)
        real(real64) :: l(n), d(n), fi, r

        point = 0
        l = 0.0d0
        d = 0.0d0

        f = 0.0d0
        ped(1, :) = RecPed(1, :)
        ped(2, :) = RecPed(2, :)
        ped(3, :) = RecPed(3, :)

        f(0) = -1.0d0
        do i = 1, n
          is = RecPed(2, i)
          id = RecPed(3, i)
          ped(2, i) = max(is, id)
          ped(3, i) = min(is, id)
          d(i) = 0.5d0 - 0.25d0 * (f(is) + f(id))
          if (is .eq. 0 .or. id .eq. 0) then
            f(i) = 0.0d0
          else if ((ped(2, i-1) .eq. ped(2, i)) .and. (ped(3, i-1) .eq. ped(3, i))) then
            f(i) = f(i-1)
          else
            fi = -1.0d0
            l(i) = 1.0d0
            j = i

            do while (j .ne. 0)
              k = j
              r = 0.5d0 * l(k)
              ks = ped(2, k)
              kd = ped(3, k)
              if (ks .gt. 0) then
                do while (point(k) .gt. ks)
                  k = point(k)
                end do
                l(ks) = l(ks) + r
                if (ks .ne. point(k)) then
                  point(ks) = point(k)
                  point(k) = ks
                end if
                if (kd .gt. 0) then
                  do while (point(k) .gt. kd)
                    k = point(k)
                  end do
                  l(kd) = l(kd) + r
                  if (kd .ne. point(k)) then
                    point(kd) = point(k)
                    point(k) = kd
                  end if
                end if
              end if
              fi = fi + l(j) * l(j) * d(j)
              l(j) = 0.0d0
              k = j
              j = point(j)
              point(k) = 0
            end do

            f(i) = fi
          end if
        end do
      end function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk & John Hickey, john.hickey@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !> @todo Metafounders and/or inbreeding for animals with unknown parents
      !-------------------------------------------------------------------------
      pure function PedNrm(RecPed, n) result(Nrm)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3, 0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n                !< Number of individuals in pedigree
        real(real64) :: Nrm(0:n, 0:n)                  !< @return Pedigree NRM

        ! Other
        integer(int32) :: Ind1, Ind2, Par1, Par2

        Nrm = 0.0d0
        do Ind2 = 1, n
            Par1 = min(RecPed(2, Ind2), RecPed(3, Ind2))
            Par2 = max(RecPed(3, Ind2), RecPed(2, Ind2))
            do Ind1 = 1, Ind2 - 1
                Nrm(Ind1, Ind2) = 0.50d0 * (Nrm(Ind1, Par1) + Nrm(Ind1, Par2))
                ! @todo would not need this memory jump if the Nrm would be considered as symmetric,
                !       but would also have to make sure algorithm is using only one triangle
                Nrm(Ind2, Ind1) = Nrm(Ind1, Ind2)
            end do
            Nrm(Ind2, Ind2) = 1.0d0 + 0.5d0 * Nrm(Par2, Par1)
        end do
      end function

      !###########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM times a vector using the Colleau
      !!         (2002, GSE 34: 409-421) method - adapted the code after
      !!         Aguilar et al. (2011, JABG 128: 422-428)
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 31, 2016
      !-------------------------------------------------------------------------
      pure function PedNrmTimesVector(RecPed, n, Inbreeding, Vector) result(Result)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3, 0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n                !< Number of individuals in pedigree
        real(real64), intent(in) :: Inbreeding(0:n)    !< Pedigree inbreeding coefficients; note Inbreeding(0) must be -1.0!
        real(real64), intent(in) :: Vector(0:n)        !< Vector to multiply NRM with
        real(real64) :: Result(0:n)                    !< @return PedNrm*Vector, i.e., Ax=b

        ! Other
        integer(int32) :: Ind, Par1, Par2
        real(real64) :: q(0:n), VarM, Tmp

        Result = 0.0d0
        q = 0.0d0

        do Ind = n, 1, -1
          q(Ind) = q(Ind) + Vector(Ind)
          Tmp = 0.5d0 * q(Ind)
          Par1 = min(RecPed(2, Ind), RecPed(3, Ind))
          Par2 = max(RecPed(3, Ind), RecPed(2, Ind))
          q(Par1) = q(Par1) + Tmp
          q(Par2) = q(Par2) + Tmp
        end do

        do Ind = 1, n
          Par1 = min(RecPed(2, Ind), RecPed(3, Ind))
          Par2 = max(RecPed(3, Ind), RecPed(2, Ind))
          ! Variance of founder effects and Mendelian sampling terms
          VarM = 0.5d0 - 0.25d0 * (Inbreeding(Par1) + Inbreeding(Par2))
          Tmp = 0.0d0
          Tmp = Tmp + Result(Par1) + Result(Par2)
          Result(Ind) = 0.5d0 * Tmp + VarM * q(Ind)
        end do
      end function

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM with an old NRM as input
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 31, 2016
      !-------------------------------------------------------------------------
      ! @todo: would be the best to write code such that we work with one Nrm
      !       and can then have say 1) base Nrm and this function kind of behaves
      !       like PedNrm, or 2) provide Nrm for the middle of pedigree and then
      !       this function gives Nrm for these individuals and the rest of pedigree
      ! pure function PedNrmForANewGeneration(RecPed, n, nNew, OldNrm, nOld, MinOldId, MaxOldId) result(Nrm)
      !   implicit none
      !
      !   ! Arguments
      !   integer(int32), intent(in) :: RecPed(1:3, 0:n)      !< Sorted and recoded pedigree array (unknown parents as 0)
      !   integer(int32), intent(in) :: n                     !< Number of individuals in pedigree
      !   integer(int32), intent(in) :: nNew                  !< Number of new generation individuals
      !   real(real64), intent(in) :: OldNrm(0:nOld, 0:nOld)  !< Old NRM
      !   integer(int32), intent(in) :: nOld                  !< Number of individuals in the old NRM
      !   integer(int32), intent(in) :: MinOldId              !< Minimal sequential id for the old NRM
      !   integer(int32), intent(in) :: MaxOldId              !< Maximal sequential id for the old NRM
      !   real(real64) :: Nrm(0:n, 0:n)                       !< @return Pedigree NRM for the new generation individuals
      !
      !   ! Other
      !   integer(int32) :: Ind1, Ind2, Par1, Par2
      !
      !   ! NOTE: The code assumes that the Nrm calculated pertains to a new generation
      !   ! of individuals, i.e., individuals in the OldNrm are ancestors of
      !   ! individuals in Nrm. It also assumes that ids in RecPed are sequential.
      !
      !   Nrm = 0.0d0
      !   do Ind1 = (MaxOldId + 1), n
      !       Par1 = max(RecPed(2, Ind1), RecPed(3, Ind1)) - MinOldId + 1
      !       Par2 = min(RecPed(3, Ind1), RecPed(2, Ind1)) - MinOldId + 1
      !       do Ind2 = 1, Ind1 - 1
      !           Nrm(Ind2, Ind1) = 0.5d0 * (Nrm(Ind2, Par1) + Nrm(Ind2, Par2))
      !           Nrm(Ind1, Ind2) = Nrm(Ind2, Ind1)
      !       end do
      !       Nrm(Ind1, Ind1) = 1.0d0 + 0.5d0 * Nrm(Par1, Par2)
      !   end do
      ! end function
      !
      ! @todo: Cleanup this old code
      !   integer(int32) :: i,j,k,l,m,n,s,d,div,MinId,MaxId,Start,Endin
      !     open(newunit=OldAMatUnit, file=OldAMatFile, status="unknown")
      !     allocate(OldAMatId(OldAMatNInd))
      !     do j = 1, OldAMatNInd
      !       read(OldAMatUnit, *) OldAMatId(j)
      !     end do
      !     rewind(OldAMatUnit)
      !     MinId = minval(OldAMatId)
      !     MaxId = maxval(OldAMatId)
      !     allocate(AMat(1:(OldAMatNInd+nAnisP-MaxId),&
      !                   1:(OldAMatNInd+nAnisP-MaxId)))
      !     do j = 1, OldAMatNInd
      !       read(OldAMatUnit, *) OldAMatId(j), AMat(1:OldAMatNInd,j)
      !       if (j > 1) then
      !         if (.not.(OldAMatId(j) > OldAMatId(j-1))) then
      !           print *, "Id are not sequential!"
      !           stop 1
      !         end if
      !       end if
      !     end do
      !     k = OldAMatNInd
      !     do i=MaxId+1,nAnisP
      !         k = k + 1
      !         s = RecPed(i,2) - MinId + 1
      !         d = RecPed(i,3) - MinId + 1
      !         l = k
      !         do j=1,k-1
      !             AMat(j,k)=(AMat(j,s)+AMat(j,d))/2.0d0
      !             AMat(k,j)=AMat(j,k)
      !             !print *,i,k,j,s,d,AMat(j,s),AMat(j,d),AMat(j,k)
      !         end do
      !         AMat(k,k)=1.0d0+AMat(s,d)/2.0d0
      !         !print *,i,k,s,d,AMat(s,d),AMat(k,k)
      !     end do

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM inverse
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk & John Hickey, john.hickey@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure function PedNrmInv(RecPed, n, Inbreeding) result(NrmInv)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3,0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n               !< Number of individuals in pedigree
        real(real64), intent(in) :: Inbreeding(0:n)   !< Pedigree inbreeding coefficients; note Inbreeding(0) must be -1.0!
        real(real64) :: NrmInv(0:n, 0:n)              !< @return Pedigree NRM inverse

        ! Other
        integer(int32) :: Ind, Par1, Par2
        real(real64) :: PreM

        NrmInv = 0.0d0
        do Ind = 1, n
          Par1 = RecPed(2, Ind)
          Par2 = RecPed(3, Ind)
          ! Precision (1/variance) of founder effects and Mendelian sampling terms
          ! PreM = 1.0d0 / (1.0d0 - 0.25d0 * (1.0d0 + Inbreeding(Par1)) - 0.25d0 * (1.0d0 + Inbreeding(Par2)))
          PreM = 1.0d0 / (0.5d0 - 0.25d0 * (Inbreeding(Par1) + Inbreeding(Par2)))
          ! Precision for the individual
          NrmInv(Ind,Ind) = PreM
          ! Add precision to the first parent and set the co-precision
          NrmInv(Par1, Par1) = NrmInv(Par1, Par1) + 0.25d0 * PreM
          NrmInv(Ind, Par1)  = NrmInv(Ind, Par1)  - 0.50d0 * PreM
          NrmInv(Par1, Ind)  = NrmInv(Ind, Par1)
          ! Add precision to the second parent and set the co-precision
          NrmInv(Par2, Par2) = NrmInv(Par2, Par2) + 0.25d0 * PreM
          NrmInv(Ind, Par2)  = NrmInv(Ind, Par2)  - 0.50d0 * PreM
          NrmInv(Par2, Ind)  = NrmInv(Ind, Par2)
          ! Add co-precision between the parents
          NrmInv(Par1, Par2) = NrmInv(Par1, Par2) + 0.25d0 * PreM
          NrmInv(Par2, Par1) = NrmInv(Par1, Par2)
        end do
        ! Reset the "margins"
        ! (the above algorithm does not need ifs for testing unknown parents as it
        !  relies on using the zeroth "margin" and Inbreeding(0)=-1. as placeholders;
        !  so should clear the "margin" now)
        NrmInv(0:n, 0) = 0.0d0
        NrmInv(0, 0:n) = 0.0d0
      end function

      !#########################################################################

    !###########################################################################

end module

!###############################################################################
