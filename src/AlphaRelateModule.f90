
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

  implicit none

  ! @todo: clean-up this code
  ! integer(int32) :: AllFreqSelCycle
  ! integer(int32) :: nGMat, OldAMatNInd
  ! integer(int32),allocatable :: RecodeGenotypeId(:),dooutput(:)
  ! integer(int32),allocatable :: OldAMatId(:)

  ! real(real64),allocatable :: Adiag(:)
  ! real(real64),allocatable :: tZMat(:,:),AMat(:,:),InvAMat(:,:)
  ! real(real64),allocatable :: GMat(:,:,:),InvGMat(:,:,:)

  private
  ! Types
  public :: AlphaRelateTitle, AlphaRelateSpec, AlphaRelateData, Inbreeding, Nrm, IndSet, GenotypeArray
  ! Functions
  public :: PedInbreeding, PedNrm, PedNrmTimesVector, PedNrmInv, MatchId

  !> @brief AlphaRelate specifications
  type AlphaRelateSpec
    character(len=FILELENGTH) :: SpecFile, PedigreeFile, GenotypeFile!, HaplotypeFile
    character(len=FILELENGTH) :: PedNrmSubsetFile!, PedNrmOldFile, LocusWeightFile, AlleleFreqFile
    character(len=SPECOPTIONLENGTH) :: OutputFormat, GenNrmType

    logical :: SpecPresent, PedigreePresent, GenotypePresent!, HaplotypePresent
    logical :: PedNrmSubsetPresent, PedNrmOldPresent!, LocusWeightPresent, AlleleFreqPresent, AlleleFreqFixed

    logical :: PedInbreeding, PedNrm, PedNrmIja, PedNrmInv, PedNrmInvIja
    logical :: GenInbreeding, GenNrm, GenNrmIja, GenNrmInv, GenNrmInvIja
    ! logical :: HapInbreeding, HapNrm, HapNrmIja, HapNrmInv, HapNrmInvIja
    ! logical :: FudgeGenNrmDiag, BlendGenNrm, FudgeHapNrmDiag, BlendHapNrm

    integer(int32):: nLoc!, nTrait, nGenMat

    !real(real64):: AlleleFreqAll
    !real(real64):: FudgeGenNrmDiagFactor, BlendGenNrmFactor, FudgeHapNrmDiagFactor, BlendHapNrmFactor
    contains
      procedure :: Init => InitAlphaRelateSpec
  end type

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
    real(real64), allocatable, dimension(:)            :: AlleleFreq
    contains
      procedure :: Init           => InitGenotypeArray
      procedure :: Destroy        => DestroyGenotypeArray
      procedure :: MatchId        => MatchIdGenotypeArray
      procedure :: Write          => WriteGenotypeArray
      procedure :: Read           => ReadGenotypeArray
      procedure :: CalcAlleleFreq => CalcAlleleFreqGenotypeArray
  end type

  !> @brief AlphaRelate data
  type AlphaRelateData
    ! Pedigree-based
    type(RecodedPedigreeArray) :: RecPed
    type(Inbreeding)           :: PedInbreeding
    type(Nrm)                  :: PedNrm
    type(IndSet)               :: PedNrmSubset
    type(Nrm)                  :: PedNrmInv
    ! Genotype-based
    type(GenotypeArray)        :: Gen
    ! type(Inbreeding)           :: GenInbreeding
    type(Nrm)                  :: GenNrm
    type(Nrm)                  :: GenNrmInv
    ! Haplotype-based
    ! type(Inbreeding)           :: HapInbreeding
    ! type(Nrm)                  :: HapNrm
    ! type(Nrm)                  :: HapNrmInv

    ! @todo: cleanup this code
    !integer(int32):: nAnisRawPedigree, nAnisP, nAnisG, nAnisH, nTrait
    !integer(int32), allocatable :: MapAnimal(:)

    !real(real64), allocatable :: AlleleFreq(:)
    !real(real64), allocatable :: ZMat(:,:), LocusWeight(:,:)

    !logical, allocatable :: MapToG(:), AnimalsInBoth(:)

    contains
      procedure :: Init              => InitAlphaRelateData
      procedure :: Destroy           => DestroyAlphaRelateData
      procedure :: Write             => WriteAlphaRelateData
      procedure :: CalcPedInbreeding
      procedure :: CalcPedNrm
      procedure :: CalcPedNrmInv
  end type

  contains

    !###########################################################################

    subroutine AlphaRelateTitle
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
        class(IndSet), intent(inout) :: This                                   !< @return IndSet holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(:)         !< Original Id of individuals in the      set
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
        if (present(OriginalId)) then
          This%OriginalId(0) = EMPTYID
          This%OriginalId(1:nInd) = OriginalId
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        if (present(OriginalIdSuperset)) then
          This%Id(0) = 0
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id = [(Ind, Ind=0, nInd)]
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
      subroutine WriteIndSet(This, File)
        implicit none
        class(IndSet), intent(in) :: This              !< IndSet holder
        character(len=*), intent(in), optional :: File !< File that will hold a set of Original Id and internal integer sequence

        integer(int32) :: Unit, Ind
        if (present(File)) then
          open(newunit=Unit, file=trim(File), status="unknown")
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
      subroutine ReadIndSet(This, File)
        implicit none
        class(IndSet), intent(out) :: This   !< @return IndSet holder
        character(len=*), intent(in) :: File !< File that holds a set of Original Id (internal integer sequence is not read)

        integer(int32) :: nInd, Ind, Unit

        nInd = CountLines(trim(File))
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
        class(Inbreeding), intent(inout) :: This                               !< @return Inbreeding holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(nInd)      !< Original Id of individuals in the      set
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip
        real(real64), intent(in), optional :: InbInput(nInd)                   !< Inbreeding coefficients

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
        if (present(OriginalId)) then
          This%OriginalId(0) = EMPTYID
          This%OriginalId(1:nInd) = OriginalId
        else
          This%OriginalId = EMPTYID
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        if (present(OriginalIdSuperset)) then
          This%Id(0) = 0
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id = [(Ind, Ind=0, nInd)]
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
      subroutine WriteInbreeding(This, File, OutputFormat)
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
          open(newunit=Unit, file=trim(File), status="unknown")
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
      subroutine ReadInbreeding(This, File)
        implicit none
        class(Inbreeding), intent(out) :: This !< @return Inbreeding holder
        character(len=*), intent(in) :: File   !< File that holds Original Id and pedigree inbreeding

        integer(int32) :: nInd, Ind, Unit

        nInd = CountLines(trim(File))
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
        class(Nrm), intent(inout) :: This                                      !< @return Nrm holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(nInd)      !< Original Id of individuals in the      set
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip
        real(real64), intent(in), optional :: NrmInput(nInd, nInd)             !< Relationship coefficients

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
        if (present(OriginalId)) then
          This%OriginalId(0) = EMPTYID
          This%OriginalId(1:nInd) = OriginalId
        else
          This%OriginalId = EMPTYID
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        if (present(OriginalIdSuperset)) then
          This%Id(0) = 0
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id = [(Ind, Ind=0, nInd)]
        end if

        if (allocated(This%Nrm)) then
          deallocate(This%Nrm)
        end if
        allocate(This%Nrm(0:nInd, 0:nInd))
        if (present(NrmInput)) then
          This%Nrm(0:nInd, 0) = 0.0d0
          This%Nrm(0, 0:nInd) = 0.0d0
          This%Nrm(1:nInd, 1:nInd) = NrmInput
        else
          This%Nrm = 0.0d0
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
      subroutine WriteNrm(This, File, Ija, OutputFormat)
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
          open(newunit=Unit, file=trim(File), status="unknown")
        else
          Unit = STDOUT
        end if
        if (IjaInternal) then
          ! Original Ids
          if (present(File)) then
            open(newunit=Unit2, file=trim(File)//"_IdMap", status="unknown")
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
      subroutine ReadNrm(This, File, Ija)
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
        class(GenotypeArray), intent(inout) :: This                            !< @return GenotypeArray holder
        integer(int32), intent(in) :: nInd                                     !< Number of individuals in the set
        integer(int32), intent(in) :: nLoc                                     !< Number of loci in the set
        character(len=IDLENGTH), intent(in), optional :: OriginalId(nInd)      !< Original Id of individuals in the      set
        character(len=IDLENGTH), intent(in), optional :: OriginalIdSuperset(:) !< Original Id of individuals in the superset
        integer(int32), intent(in), optional :: Skip                           !< How many elements of OriginalIdSuperset to skip
        integer(int8), intent(in), optional :: IntegerInput(nLoc, nInd)        !< Genotypes as an array of integers
        type(Genotype), intent(in), optional :: GenotypeInput(nInd)            !< Genotypes as genotype type

        ! Other
        integer(int32) :: Ind, Start
        integer(int8) :: TempInt(nLoc)
        type(Genotype) :: TempGeno

        if (present(Skip)) then
          Start = Skip + 1
        else
          Start = 1
        end if

        This%nInd = nInd
        This%nLoc = nLoc

        if (allocated(This%OriginalId)) then
          deallocate(This%OriginalId)
        end if
        allocate(This%OriginalId(0:nInd))
        if (present(OriginalId)) then
          This%OriginalId(0) = EMPTYID
          This%OriginalId(1:nInd) = OriginalId
        else
          This%OriginalId = EMPTYID
        end if

        if (allocated(This%Id)) then
          deallocate(This%Id)
        end if
        allocate(This%Id(0:nInd))
        if (present(OriginalIdSuperset)) then
          This%Id(0) = 0
          This%Id(1:nInd) = MatchId(IdSet=This%OriginalId(1:nInd),& ! to handle "0th margin"
                                    IdSuperset=OriginalIdSuperset(Start:size(OriginalIdSuperset))) ! to handle potential "0th margin"
        else
          This%Id = [(Ind, Ind=0, nInd)]
        end if

        if (allocated(This%Genotype)) then
          deallocate(This%Genotype)
        end if
        allocate(This%Genotype(0:nInd))
        TempInt = MISSINGGENOTYPECODE
        TempGeno = Genotype(Geno=TempInt)
        if (present(IntegerInput) .or. present(GenotypeInput)) then
          This%Genotype(0) = TempGeno
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
          do Ind = 0, nInd
            This%Genotype(Ind) = TempGeno
          end do
        end if

        if (allocated(This%AlleleFreq)) then
          deallocate(This%AlleleFreq)
        end if
        allocate(This%AlleleFreq(nLoc))
        This%AlleleFreq = 0.0d0
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
      !> @brief  Write GenotypeArray to a file or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine WriteGenotypeArray(This, File)
        implicit none
        class(GenotypeArray), intent(in) :: This       !< GenotypeArray Holder
        character(len=*), intent(in), optional :: File !< File that will hold Genotypes

        integer(int32) :: Unit, Ind
        character(len=:), allocatable :: Fmt

        if (present(File)) then
          open(newunit=Unit, file=trim(File), status="unknown")
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
      !> @brief  Read GenotypeArray from a file
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine ReadGenotypeArray(This, File, nLoc)
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
      !> @brief  Calculate allele frequencies on GenotypeArray
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine CalcAlleleFreqGenotypeArray(This)
        implicit none
        class(GenotypeArray), intent(inout) :: This !< @return GenotypeArray holder

        integer(int32) :: Ind, Loc
        integer(int32), allocatable, dimension(:) :: nObs
        real(real64), allocatable, dimension(:) :: Geno

        This%AlleleFreq = 0.0d0

        allocate(nObs(This%nLoc))
        allocate(Geno(This%nLoc))
        nObs = 0
        do Ind = 1, This%nInd
          Geno = dble(This%Genotype(Ind)%ToIntegerArray())
          do Loc = 1, This%nLoc
            if ((Geno(Loc) .ge. 0.0d0) .and. (Geno(Loc) .le. 2.0d0)) then
              This%AlleleFreq(Loc) = This%AlleleFreq(Loc) + Geno(Loc)
              nObs(Loc) = nObs(Loc) + 1
            end if
          end do
        end do
        do Loc = 1, This%nLoc
          if (nObs(Loc) .gt. 0) then
            This%AlleleFreq(Loc) = This%AlleleFreq(Loc) / (2.0d0 * dble(nObs(Loc)))
          else
            This%AlleleFreq(Loc) = 0.0d0
          end if
        end do
        deallocate(Geno)
        deallocate(nObs)
      end subroutine

      !#########################################################################

    !###########################################################################

    ! AlphaRelateSpec type methods

      !-------------------------------------------------------------------------
      !> @brief  AlphaRelateSpec constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine InitAlphaRelateSpec(This, SpecFile)
        implicit none

        ! Arguments
        class(AlphaRelateSpec), intent(inout) :: This      !< @return AlphaRelateSpec holder
        character(len=*), intent(in), optional :: SpecFile !< Spec file; when missing, a stub with defaults is created

        ! Other
        character(len=:), allocatable :: DumString
        character(len=SPECOPTIONLENGTH) :: Line
        character(len=SPECOPTIONLENGTH) :: First
        character(len=SPECOPTIONLENGTH), allocatable, dimension(:) :: Second

        integer(int32) :: SpecUnit, Stat

        ! Defaults
        This%SpecFile         = "None"
        This%PedigreeFile     = "None"
        This%PedNrmSubsetFile = "None"
        ! This%PedNrmOldFile    = "None"
        ! This%GenotypeFile     = "None"
        ! This%HaplotypeFile    = "None"
        ! This%LocusWeightFile  = "None"
        ! This%AlleleFreqFile   = "None"

        ! This%GenNrmType       = "None"
        This%OutputFormat  = "f"

        This%SpecPresent         = .false.
        This%PedigreePresent     = .false.
        This%PedNrmSubsetPresent = .false.
        This%PedNrmOldPresent    = .false.
        This%GenotypePresent     = .false.
        ! This%HaplotypePresent    = .false.
        ! This%LocusWeightPresent  = .false.
        ! This%AlleleFreqPresent   = .false.
        ! This%AlleleFreqFixed     = .false.

        This%PedInbreeding       = .false.
        This%PedNrm              = .false.
        This%PedNrmIja           = .false.
        This%PedNrmInv           = .false.
        This%PedNrmInvIja        = .false.

        ! This%GenInbreeding       = .false.
        ! This%GenNrm              = .false.
        ! This%GenNrmIja           = .false.
        ! This%GenNrmInv           = .false.
        ! This%GenNrmInvIja        = .false.
        ! This%FudgeGenNrmDiag     = .false.
        ! This%BlendGenNrm         = .false.

        ! This%HapInbreeding       = .false.
        ! This%HapNrm              = .false.
        ! This%HapNrmIja           = .false.
        ! This%HapNrmInv           = .false.
        ! This%HapNrmInvIja        = .false.
        ! This%FudgeHapNrmDiag     = .false.
        ! This%BlendHapNrm         = .false.

        This%nLoc             = 0
        ! This%nTrait           = 1
        ! This%nGenMat          = 0

        ! This%AlleleFreqAll         = 0.5d0
        ! This%FudgeHapNrmDiagFactor = 0.0d0
        ! This%BlendHapNrmFactor     = 0.0d0
        ! This%FudgeHapNrmDiagFactor = 0.0d0
        ! This%BlendHapNrmFactor     = 0.0d0

        if (present(SpecFile)) then

          This%SpecPresent = .true.
          This%SpecFile = SpecFile
          open(newunit=SpecUnit, file=trim(This%SpecFile), action="read", status="old")

          Stat = 0
          ReadSpec: do while (Stat == 0)
            read(SpecUnit, "(a)", iostat=Stat) Line
            if (len_trim(Line) == 0) then
              cycle
            end if
            call SplitLineIntoTwoParts(trim(Line), First, Second)
            DumString = ParseToFirstWhitespace(First)
            ! @todo why (len_trim(Line) == 0)? if we use (len_trim(Line) == 0) above
            if (First(1:1) == "=" .or. len_trim(Line) == 0) then
              cycle
            else
              select case (ToLower(trim(DumString)))

                case ("pedigreefile")
                  if (allocated(Second)) then
                    if (ToLower(trim(Second(1))) == "none") then
                      write(STDOUT, "(a)") " Not using pedigree file"
                    else
                      This%PedigreePresent = .true.
                      write(This%PedigreeFile, *) trim(Second(1))
                      write(STDOUT, "(2a)") " Using pedigree file: ", trim(This%PedigreeFile)
                    end if
                  else
                    write(STDERR, "(a)") " ERROR: Must specify a file for PedigreeFile, i.e., PedigreeFile, Pedigree.txt"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if

                case ("genotypefile")
                  if (allocated(Second)) then
                    if (ToLower(trim(Second(1))) == "none") then
                      write(STDOUT, "(a)") " Not using genotype file"
                    else
                      This%GenotypePresent = .true.
                      write(This%GenotypeFile, *) trim(Second(1))
                      write(STDOUT, "(2a)") " Using genotype file: ", trim(This%GenotypeFile)
                    end if
                  else
                    write(STDERR, "(a)") " ERROR: Must specify a file for GenotypeFile, i.e., GenotypeFile, Genotype.txt"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if

                ! case ("haplotypefile")
                !   if (allocated(Second)) then
                !     if (ToLower(trim(Second(1))) == "none") then
                !       write(STDOUT, "(a)") " Not using haplotype file"
                !     else
                !       This%HaplotypePresent = .true.
                !       write(This%HaplotypeFile, *) trim(Second(1))
                !       write(STDOUT, "(2a)") " Using haplotype file: ", trim(This%HaplotypeFile)
                !     end if
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a file for HaplotypeFile, i.e., HaplotypeFile, Haplotype.txt"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("locusweightfile")
                !   if (allocated(Second)) then
                !     if (ToLower(trim(Second(1))) == "none") then
                !       write(STDOUT, "(a)") " Not using locus weights file"
                !     else
                !       This%LocusWeightPresent = .true.
                !       write(This%LocusWeightFile, *) trim(Second(1))
                !       write(STDOUT, "(2a)") " Using locus weight file: ", trim(This%LocusWeightFile)
                !     end if
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a file for LocusWeightFile, i.e., LocusWeightFile, LocusWeight.txt"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("allelefreqfile")
                !   if (allocated(Second)) then
                !     if (ToLower(trim(Second(1))) == "none") then
                !       write(STDOUT, "(a)") " Not using precalculated/fixed allele frequencies file"
                !     else
                !       This%AlleleFreqPresent = .true.
                !       if (ToLower(trim(Second(1))) == "fixed") then
                !         This%AlleleFreqFixed = .true.
                !         if (size(Second) > 1) then
                !           This%AlleleFreqAll = Char2Double(trim(Second(2)), "(f20.16)")
                !         else
                !           This%AlleleFreqAll = 0.5d0
                !         end if
                !         write(STDOUT, "(2a)") " Using fixed allele frequency: ", Real2Char(This%AlleleFreqAll, "(f6.4)")
                !       else
                !         write(This%AlleleFreqFile, *) trim(Second(1))
                !         write(STDOUT, "(2a)") " Using allele frequencies file: ", trim(This%AlleleFreqFile)
                !       end if
                !     end if
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a file for AlleleFreqFile, i.e., AlleleFreqFile, AlleleFreq.txt"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                case ("numberofloci")
                  if (allocated(Second)) then
                    This%nLoc = Char2Int(trim(Second(1)))
                    write(STDOUT, "(a, i)") " Number of loci: ", This%nLoc
                  else
                    write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfLoci, i.e., NumberOfLoci, 10"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if

                ! case ("numberoftraits")
                !   if (allocated(Second)) then
                !     This%nTrait = Char2Int(trim(Second(1)))
                !     write(STDOUT, "(2a)") " Number of traits: ", trim(This%nTrait)
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfTraits, i.e., NumberOfTraits, 1"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("gennrmtype")
                !   if (allocated(Second)) then
                !     write(This%GenNrmType, *) ToLower(trim(Second(1)))
                !     if (trim(This%GenNrmType) /= "vanraden"        .and. &
                !         trim(This%GenNrmType) /= "vanraden1"       .and. &
                !         trim(This%GenNrmType) /= "vanraden2"       .and. &
                !         trim(This%GenNrmType) /= "yang"            .and. &
                !         trim(This%GenNrmType) /= "nejati-javaremi") then
                !         ! trim(This%GenNrmType) /= "day-williams") then
                !       write(STDERR, "(a)") " ERROR: GenNrmType must be either VanRaden=VanRaden1, VanRaden2, Yang, or Nejati-Javaremi"
                !       write(STDERR, "(a)") ""
                !       stop 1
                !     end if
                !     write(STDOUT, "(2a)") " Genotype NRM type: ", trim(This%GenNrmType)
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a method for GenNrmType, i.e., GenNrmType, VanRaden"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("fudgegennrmdiag")
                !   if (allocated(Second)) then
                !     This%FudgeGenNrmDiag = .true.
                !     This%FudgeGenNrmDiagFactor = Char2Double(trim(Second(1)), "(f20.16)")
                !     write(STDOUT, "(2a)") " Fudge genotype NRM diagonal: ", Real2Char(This%FudgeGenNrmDiagFactor, "(f6.4)")
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a number for FudgeGenNrmDiag, i.e., FudgeGenNrmDiag, 0.001"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("blendgennrm")
                !   if (allocated(Second)) then
                !     This%PedNrm = .true.
                !     This%BlendGenNrm = .true.
                !     This%BlendGenNrmFactor = Char2Double(trim(Second(1)), "(f20.16)")
                !     write(STDOUT, "(2a)") " Blend genotype NRM: ", Real2Char(This%BlendGenNrmFactor, "(f6.4)")
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a number for BlendGenNrm, i.e., BlendGenNrm, 0.95"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("fudgehapnrmdiag")
                !   if (allocated(Second)) then
                !     This%FudgeHapNrmDiag = .true.
                !     This%FudgeHapNrmDiagFactor = Char2Double(trim(Second(1)), "(f20.16)")
                !     write(STDOUT, "(2a)") " Fudge haplotype NRM diagonal: ", Real2Char(This%FudgeHapNrmDiagFactor, "(f6.4)")
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a number for FudgeHapNrmDiag, i.e., FudgeHapNrmDiag, 0.001"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                ! case ("blendhapnrm")
                !   if (allocated(Second)) then
                !     This%PedNrm = .true.
                !     This%BlendHapNrm = .true.
                !     This%BlendHapNrmFactor = Char2Double(trim(Second(1)), "(f20.16)")
                !     write(STDOUT, "(2a)") " Blend haplotype NRM: ", Real2Char(This%BlendHapNrmFactor, "(f6.4)")
                !   else
                !     write(STDERR, "(a)") " ERROR: Must specify a number for BlendHapNrm, i.e., BlendHapNrm, 0.95"
                !     write(STDERR, "(a)") ""
                !     stop 1
                !   end if

                case ("outputformat")
                  if (allocated(Second)) then
                    write(This%OutputFormat, *) trim(Second(1))
                    write(STDOUT, "(2a)") " Output precision: ", trim(This%OutputFormat)
                  else
                    write(STDERR, "(a)") " ERROR: Must specify Fortran format for OutputFormat, i.e., OutputFormat, f16.8"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if

                case ("pedinbreeding")
                  if (allocated(Second)) then
                    if (ToLower(trim(Second(1))) == "yes") then
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
                    if (ToLower(trim(Second(1))) == "yes") then
                      This%PedNrm = .true.
                      write(STDOUT, "(a)") " Calculate pedigree NRM: Yes"
                      if (size(Second) > 1) then
                        if (ToLower(trim(Second(2))) == "ija") then
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
                    if (ToLower(trim(Second(1))) == "none") then
                      write(STDOUT, "(a)") " Not using pedigree NRM set file"
                    else
                      This%PedNrmSubsetPresent = .true.
                      write(This%PedNrmSubsetFile, *) trim(Second(1))
                      write(STDOUT, "(2a)") " Using pedigree NRM set file: ", trim(This%PedNrmSubsetFile)
                    end if
                  else
                    write(STDERR, "(a)") " ERROR: Must specify a file for PedNrmSubsetFile, i.e., PedNrmSubsetFile, PedNrmSubset.txt"
                    write(STDERR, "(a)") ""
                    stop 1
                  end if

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
                    if (ToLower(trim(Second(1))) == "yes") then
                      This%PedNrmInv = .true.
                      write(STDOUT, "(a)") " Calculate pedigree NRM inverse: Yes"
                      if (size(Second) > 1) then
                        if (ToLower(trim(Second(2))) == "ija") then
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

                ! read(SpecUnit,*) DumC, Option
                ! This%GFullMat = trim(Option) == "Yes"

                ! read(SpecUnit,*) DumC, Option
                ! This%GIJA = trim(Option) == "Yes"

                ! if (This%GFullMat .or. This%GIJA) then
                !   This%MakeG = .true.
                ! end if

                ! read(SpecUnit,*) DumC, Option
                ! This%InvGFullMat = trim(Option) == "Yes"

                ! read(SpecUnit,*) DumC, Option
                ! This%InvGIJA = trim(Option) == "Yes"

                ! if (This%InvGFullMat .or. This%InvGIJA) then
                !   This%MakeInvG = .true.
                ! end if

                ! read(SpecUnit,*) DumC, Option
                ! This%HFullMat = trim(Option) == "Yes"

                ! read(SpecUnit,*) DumC, Option
                ! This%HIJA = trim(Option) == "Yes"

                ! if (This%HFullMat .or. This%HIJA) then
                !   This%MakeH = .true.
                !   This%MakeG = .true.
                !   This%MakeA = .true.
                ! end if

                ! read(SpecUnit,*) DumC, Option
                ! This%InvHFullMat = trim(Option) == "Yes"

                ! read(SpecUnit,*) DumC, Option
                ! This%InvHIJA = trim(Option) == "Yes"

                ! if (This%InvHFullMat .or. This%InvHIJA) then
                !   This%MakeInvH = .true.
                !   This%MakeG    = .true.
                !   This%MakeA    = .true.
                !   This%MakeInvA = .true.
                ! end if

                case default
                  write(STDOUT, "(3a)") " NOTE: Specification '", trim(Line), "' ignored"
                  write(STDOUT, "(a)") " "
              end select
            end if
          end do ReadSpec
          close(SpecUnit)

          if ((This%PedInbreeding .or. This%PedNrm .or. This%PedNrmInv .or. This%PedNrmSubsetPresent)&
              .and. .not. This%PedigreePresent) then
            write(STDERR, "(a)") " ERROR: Must provide pedigree file to calculate pedigree inbreeding, NRM, or NRM inverse"
            write(STDERR, "(a)") ""
            stop 1
          end if

          ! if ((This%GenInbreeding .or. This%GenNrm .or. This%GenNrmInv)&
          !     .and. .not. This%GenotypePresent) then
          !   write(STDERR, "(a)") " ERROR: Must provide genotype file to calculate genotype inbreeding, NRM, or NRM inverse"
          !   write(STDERR, "(a)") ""
          !   stop 1
          ! end if

          ! if ((This%HapInbreeding .or. This%HapNrm .or. This%HapNrmInv)&
          !     .and. .not. This%HaplotypePresent) then
          !   write(STDERR, "(a)") " ERROR: Must provide haplotype file to calculate haplotype inbreeding, NRM, or NRM inverse"
          !   write(STDERR, "(a)") ""
          !   stop 1
          ! end if

          ! if (This%BlendGenNrm .and. .not. This%PedigreePresent) then
          !   write(STDERR, "(a)") " ERROR: Must provide pedigree file to blend genotype NRM with pedigree NRM"
          !   write(STDERR, "(a)") ""
          !   stop 1
          ! end if

          ! if (This%BlendHapNrm .and. .not. This%PedigreePresent) then
          !   write(STDERR, "(a)") " ERROR: Must provide pedigree file to blend haplotype NRM with pedigree NRM"
          !   write(STDERR, "(a)") ""
          !   stop 1
          ! end if

          ! This%nGenMat=0
          ! do i = 1, This%nTrait
          !   do j = i, This%nTrait
          !     This%nGenMat = This%nGenMat + 1
          !   end do
          ! end do

          ! if ((This%MakeG .or. This%MakeInvG .or. This%MakeH .or. This%MakeInvH) .and. .not. This%GenotypePresent) then
          !   write(STDOUT, "(a)") " NOTE: To create G or H matrix, a genotype file must be given --> ommited G or H."
          !   write(STDOUT, "(a)") " "
          !   This%MakeG    = .false.
          !   This%MakeInvG = .false.
          !   This%MakeH    = .false.
          !   This%MakeInvH = .false.
          ! end if

          ! if ((This%MakeA .or. This%MakeInvA .or. This%MakeH .or. This%MakeInvH) .and. .not. This%PedigreePresent) then
          !   write(STDOUT, "(a)") " NOTE: To create A or H matrix, a pedigree file must be given --> ommited A or H."
          !   write(STDOUT, "(a)") " "
          !   This%MakeA    = .false.
          !   This%MakeInvA = .false.
          !   This%MakeH    = .false.
          !   This%MakeInvH = .false.
          ! end if

        end if
      end subroutine

    !###########################################################################

    ! AlphaRelateData type methods

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  AlphaRelateData constructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine InitAlphaRelateData(This, Spec)
        implicit none

        ! Arguments
        class(AlphaRelateData), intent(inout) :: This  !< @return AlphaRelateData holder
        type(AlphaRelateSpec), intent(in) :: Spec      !< Specifications

        ! Other
        type(PedigreeHolder) :: PedObj
        ! @todo cleanup this code
        ! integer(int32) :: Stat, nCols, GenoInPed
        ! integer(int32) :: PedNrmOldUnit, GenotypeUnit, AlleleFreqUnit, LocusWeightUnit
        ! character(len=IDLENGTH) :: DumC

        if (Spec%PedigreePresent) then
          ! Read in the pedigree
          PedObj = PedigreeHolder(Spec%PedigreeFile)

          ! Sort and recode pedigree
          call PedObj%MakeRecodedPedigreeArray(RecPed=This%RecPed)
          write(STDOUT, "(a1, i8, a)") " ", This%RecPed%nInd," individuals in the pedigree"
          call This%RecPed%Write(File=trim(Spec%PedigreeFile)//"_Recoded.txt")

          ! Free some memory
          call PedObj%DestroyPedigree

          ! Handle subset
          if (Spec%PedNrmSubsetPresent) then
            call This%PedNrmSubset%Read(File=Spec%PedNrmSubsetFile)
            write(STDOUT, "(a1, i8, a)") " ", This%PedNrmSubset%nInd," individuals in the pedigree NRM subset"
            call This%PedNrmSubset%MatchId(OriginalIdSuperset=This%RecPed%OriginalId, Skip=1) ! skip=1 because of the "0th margin" in This%RecPed%OriginalId
            block
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

        if (Spec%GenotypePresent) then

          call This%Gen%Read(File=Spec%GenotypeFile, nLoc=Spec%nLoc)
          write(STDOUT, "(a1, i8, a)") " ", This%Gen%nInd," individuals with genotypes"

          if (Spec%PedigreePresent) then
            call This%Gen%MatchId(OriginalIdSuperset=This%RecPed%OriginalId, Skip=1) ! skip=1 because of the "0th margin" in This%RecPed%OriginalId
            block
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

        !   This%nTrait = Spec%nTrait
        !   allocate(This%ZMat(This%nAnisG,This%nLoc))
        !
        !   ! Allele frequencies
        !   allocate(This%AlleleFreq(This%nLoc))
        !   if (.not. Spec%AlleleFreqPresent) then
        !     open(newunit=AlleleFreqUnit, file="AlleleFreq.txt", action=???, status="unknown")
        !     do j = 1, This%nLoc
        !       write(AlleleFreqUnit,*) j, This%AlleleFreq(j)
        !     end do
        !     close(AlleleFreqUnit)
        !   else
        !     if (trim(Spec%AlleleFreqFile) == "Fixed") then
        !       This%AlleleFreq(:) = Spec%AlleleFreqAll
        !       open(newunit=AlleleFreqUnit, file="AlleleFreq.txt", action=???, status="unknown")
        !       do j = 1, nLoc
        !         write(AlleleFreqUnit,*) j, This%AlleleFreq(j)
        !       end do
        !       close(AlleleFreqUnit)
        !     else
        !       ! Read allele frequencies from file.
        !       open(newunit=AlleleFreqUnit, file=trim(Spec%AlleleFreqFile), action="read", status="old")
        !       do i = 1, This%nLoc
        !         ! AlleleFrequencies are kept in second column to keep consistency with AlphaSim.
        !         read(AlleleFreqUnit, *, iostat=Stat) DumC, This%AlleleFreq(i)
        !         if (Stat /= 0) then
        !           write(STDERR, "(a)") " ERROR: Problems reading allele frequency file."
        !           write(STDERR, "(a)") " "
        !           stop 1
        !         end if
        !       end do
        !       close(AlleleFreqUnit)
        !     end if
        !   end if
        !
        !   ! LocusWeight
        !   allocate(This%LocusWeight(This%nLoc,This%nTrait))
        !   if (Spec%LocusWeightPresent) then
        !     open(newunit=LocusWeightUnit, file=trim(Spec%LocusWeightFile), action="read", status="old")
        !     do i = 1, This%nLoc
        !       read(LocusWeightUnit,*) DumC, This%LocusWeight(i,:)
        !     end do
        !     close(LocusWeightUnit)
        !   else
        !     This%LocusWeight(:,:) = 1.0d0
        !   end if
        end if

        ! if (.not. Spec%PedigreePresent .and. Spec%GenotypePresent) then
        !   allocate(This%RecPed(0:This%nAnisG,4))
        !   This%nAnisP = This%nAnisG
        !   This%RecPed(:,:) = 0
        !   This%RecPed(:,4) = 1
        !   do i = 1, This%nAnisP
        !     This%RecPed(i,1) = i
        !   end do
        ! end if

        ! if (Spec%PedigreePresent .and. Spec%GenotypePresent) then
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
      !> @brief  AlphaRelateData destructor
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine DestroyAlphaRelateData(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder

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
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Write AlphaRelateData to files or stdout
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      subroutine WriteAlphaRelateData(This, Basename)
        implicit none
        class(AlphaRelateData), intent(in) :: This         !< AlphaRelateData holder
        character(len=*), intent(in), optional :: Basename !< Basename for produced files

        if (allocated(This%RecPed%OriginalId)) then
          if (present(Basename)) then
            call This%RecPed%Write(File=trim(Basename)//"RecodedPedigree.txt")
          else
            write(STDOUT, "(a)") "Recoded pedigree"
            call This%RecPed%Write
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

        if (allocated(This%PedNrmSubset%OriginalId)) then
          if (present(Basename)) then
            call This%PedNrmSubset%Write(File=trim(Basename)//"PedNrmSubset.txt")
          else
            write(STDOUT, "(a)") "Pedigree NRM subset"
            call This%PedNrmSubset%Write
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
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree inbreeding on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcPedInbreeding(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder, note This%PedInbreeding(0) = -1.0!!!
        call This%PedInbreeding%Init(nInd=This%RecPed%nInd, OriginalId=This%RecPed%OriginalId)
        This%PedInbreeding%Inb = PedInbreeding(RecPed=This%RecPed%Id, n=This%PedInbreeding%nInd)
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree inbreeding using the Meuwissen and
      !!         Luo (1992, GSE 24: 305-313) method
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk & John Hickey, john.hickey@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure function PedInbreeding(RecPed, n) result(f)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3,0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n               !< Number of individuals in pedigree
        real(real64) :: f(0:n)                        !< @return Pedigree inbreeding, note PedInbreeding(0) = -1.0!!!

        ! Other
        integer(int32) :: i, is, id, j, k, ks, kd
        integer(int32) :: ped(3,0:n), point(0:n)
        real(real64) :: l(n), d(n), fi, r

        point = 0
        l = 0.0d0
        d = 0.0d0

        f = 0.0d0
        ped(1,:) = RecPed(1,:)
        ped(2,:) = RecPed(2,:)
        ped(3,:) = RecPed(3,:)

        f(0) = -1.0d0
        do i = 1, n
          is = RecPed(2,i)
          id = RecPed(3,i)
          ped(2,i) = max(is,id)
          ped(3,i) = min(is,id)
          d(i) = 0.5d0 - 0.25d0 * (f(is) + f(id))
          if (is .eq. 0 .or. id .eq. 0) then
            f(i) = 0.0d0
          else if ((ped(2,i-1) .eq. ped(2,i)) .and. (ped(3,i-1) .eq. ped(3,i))) then
            f(i) = f(i-1)
          else
            fi = -1.0d0
            l(i) = 1.0d0
            j = i

            do while (j .ne. 0)
              k = j
              r = 0.5d0 * l(k)
              ks = ped(2,k)
              kd = ped(3,k)
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
      !> @brief  Calculate pedigree NRM on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcPedNrm(This, Spec)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder
        type(AlphaRelateSpec), intent(in) :: Spec     !< Specifications

        if (Spec%PedNrmSubsetPresent) then
          ! Using the Colleau/Tier method to get Nrm for a subset of individuals
          call This%PedNrm%Init(nInd=This%PedNrmSubset%nInd,&
                                OriginalId=This%PedNrmSubset%OriginalId)
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
        else if (Spec%PedNrmOldPresent) then
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
          call This%PedNrm%Init(nInd=This%RecPed%nInd, OriginalId=This%RecPed%OriginalId)
          This%PedNrm%Nrm = PedNrm(RecPed=This%RecPed%Id, n=This%PedNrm%nInd)
        end if
      end subroutine

      !#########################################################################

      !-------------------------------------------------------------------------
      !> @brief  Calculate pedigree NRM
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk & John Hickey, john.hickey@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure function PedNrm(RecPed, n) result(Nrm)
        implicit none

        ! Arguments
        integer(int32), intent(in) :: RecPed(1:3,0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n               !< Number of individuals in pedigree
        real(real64) :: Nrm(0:n, 0:n)                 !< @return Pedigree NRM

        ! Other
        integer(int32) :: Ind1, Ind2, Par1, Par2

        Nrm = 0.0d0
        do Ind1 = 1, n
            Par1 = max(RecPed(2, Ind1), RecPed(3, Ind1))
            Par2 = min(RecPed(3, Ind1), RecPed(2, Ind1))
            do Ind2 = 1, Ind1 - 1
                Nrm(Ind2, Ind1) = (Nrm(Ind2, Par1) + Nrm(Ind2, Par2)) / 2.0d0
                Nrm(Ind1, Ind2) = Nrm(Ind2, Ind1)
            end do
            Nrm(Ind1, Ind1) = 1.0d0 + Nrm(Par1, Par2) / 2.0d0
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
        integer(int32), intent(in) :: RecPed(1:3,0:n) !< Sorted and recoded pedigree array (unknown parents as 0)
        integer(int32), intent(in) :: n               !< Number of individuals in pedigree
        real(real64), intent(in) :: Inbreeding(0:n)   !< Pedigree inbreeding coefficients; note Inbreeding(0) must be -1.0!
        real(real64), intent(in) :: Vector(0:n)       !< Vector to multiply NRM with
        real(real64) :: Result(0:n)                   !< @return PedNrm*Vector, i.e., Ax=b

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
      !   integer(int32), intent(in) :: RecPed(1:3,0:n)      !< Sorted and recoded pedigree array (unknown parents as 0)
      !   integer(int32), intent(in) :: n                    !< Number of individuals in pedigree
      !   integer(int32), intent(in) :: nNew                 !< Number of new generation individuals
      !   real(real64), intent(in) :: OldNrm(0:nOld, 0:nOld) !< Old NRM
      !   integer(int32), intent(in) :: nOld                 !< Number of individuals in the old NRM
      !   integer(int32), intent(in) :: MinOldId             !< Minimal sequential id for the old NRM
      !   integer(int32), intent(in) :: MaxOldId             !< Maximal sequential id for the old NRM
      !   real(real64) :: Nrm(0:n, 0:n)                      !< @return Pedigree NRM for the new generation individuals
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
      !           Nrm(Ind2, Ind1) = (Nrm(Ind2, Par1) + Nrm(Ind2, Par2)) / 2.0d0
      !           Nrm(Ind1, Ind2) = Nrm(Ind2, Ind1)
      !       end do
      !       Nrm(Ind1, Ind1) = 1.0d0 + Nrm(Par1, Par2) / 2.0d0
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
      !> @brief  Calculate pedigree NRM inverse on AlphaRelateData
      !> @author Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
      !> @date   December 22, 2016
      !-------------------------------------------------------------------------
      pure subroutine CalcPedNrmInv(This)
        implicit none
        class(AlphaRelateData), intent(inout) :: This !< @return AlphaRelateData holder
        if (.not. allocated(This%PedInbreeding%Inb)) then
          call This%CalcPedInbreeding
        end if
        call This%PedNrmInv%Init(nInd=This%RecPed%nInd, OriginalId=This%RecPed%OriginalId)
        This%PedNrmInv%Nrm = PedNrmInv(RecPed=This%RecPed%Id, n=This%PedNrmInv%nInd,&
                                       Inbreeding=This%PedInbreeding%Inb)
      end subroutine

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
          NrmInv(Par1, Par1) = NrmInv(Par1, Par1) + PreM / 4.0d0
          NrmInv(Ind, Par1)  = NrmInv(Ind, Par1)  - PreM / 2.0d0
          NrmInv(Par1, Ind)  = NrmInv(Ind, Par1)
          ! Add precision to the second parent and set the co-precision
          NrmInv(Par2, Par2) = NrmInv(Par2, Par2) + PreM / 4.0d0
          NrmInv(Ind, Par2)  = NrmInv(Ind, Par2)  - PreM / 2.0d0
          NrmInv(Par2, Ind)  = NrmInv(Ind, Par2)
          ! Add co-precision between the parents
          NrmInv(Par1, Par2) = NrmInv(Par1, Par2) + PreM / 4.0d0
          NrmInv(Par2, Par1) = NrmInv(Par1, Par2)
        end do
        ! Reset the "margins"
        ! (the above algorithm does not need ifs for testing unknown parents as it
        !  relies on using the zeroth "margin" and Inbreeding(0)=-1. as placeholders;
        !  so should clear the "margin" now)
        NrmInv(0:n, 0) = 0.0d0
        NrmInv(0, 0:n) = 0.0d0
      end function

      ! @todo: is this usefull when dealing with the single-step H matrix?
      !   if (InvAFullMat) then
      !     AnimToWrite = RecPed(1:nAnisP,4) == 1
      !     s = count(AnimToWrite)
      !     write(*,"(a40,i6,a11)") " Start writing A inverse full matrix for", s," individuals"
      !     write(nChar,*) s
      !     fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
      !     open(unit=202,file="InvAFullMatrix.txt",status="unknown")
      !     do m=1,nAnisP
      !       if (AnimToWrite(m)) then
      !         write(202,fmt) Id(m), pack(InvAMat(1:nAnisP,m), AnimToWrite)
      !       end if
      !     end do
      !     close(202)
      !     print*, "End writing A inverse full matrix"
      !   end if

      !#########################################################################

      !   if (AFullMat) then
      !     AnimToWrite = RecPed(1:nAnisP,4) == 1
      !     s = count(AnimToWrite)
      !     write(*,"(a32,i6,a11)") " Start writing A full matrix for", s," individuals"
      !     write(nChar,*) s
      !     fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
      !     open(unit=202,file="AFullMatrix.txt",status="unknown")
      !     if (.not.OldAMatPresent) then
      !       do m=1,nAnisP
      !         if (AnimToWrite(m)) then
      !           write(202,fmt) Id(m), pack(AMat(1:nAnisP,m), AnimToWrite)
      !         end if
      !       end do
      !     else
      !       Start = OldAMatNInd+1
      !       Endin = size(AMat,1)
      !       do m=Start,Endin
      !         !write(*,fmt)   Id(m+MinId-1), AMat(Start:Endin,m)
      !         write(202,fmt) Id(m+MinId-1), AMat(Start:Endin,m)
      !       end do
      !     end if
      !     close(202)
      !     print*, "End writing A full matrix"
      !   end if

      !   ! Record diagonals of animals in both A and G:
      !   if ((MakeH .or. MakeInvH) .and. ScaleGByRegression) then
      !     n = Count(AnimalsInBoth)
      !     allocate(Adiag(0:n))
      !     div = dble(n**2)
      !     AMatAvg = 0.0d0
      !     k = 0
      !     do i = 1,nAnisP
      !       if (.not. AnimalsInBoth(i)) then
      !         cycle
      !       end if
      !       k = k + 1
      !       Adiag(k) = AMat(i,i)
      !       do j=1,nAnisP
      !         if (AnimalsInBoth(j)) then
      !           AMatAvg=AMatAvg + AMat(j,i) * 2.0d0 / div
      !         end if
      !       end do
      !     end do
      !     Adiag(0) = AMatAvg
      !   end if
      ! end subroutine

      !#########################################################################

      ! subroutine MakeGAndInvGMatrix
      !   implicit none

      !   integer(int32) :: i,j,k,l,m,n,WhichMat

      !   real(real64) :: nLocD, DMatSum, Tmp, Tmp2, Tmp3
      !   real(real64), allocatable :: TmpZMat(:,:), DMat(:,:)

      !   character(len=1000) :: filout,nChar,fmt

      !   allocate(GMat(nAnisG,nAnisG,nGMat))
      !   allocate(tZMat(nLoc,nAnisG))
      !   allocate(TmpZMat(nAnisG,nLoc))
      !   if (LocusWeightPresent) then
      !     allocate(DMat(nLoc,nLoc))
      !     DMat(:,:)=0.0d0
      !   end if

      !   print*, "Start making G - ", trim(GType)

      !   nLocD = dble(nLoc)

      !   ! Center allele dosages (Z)
      !   if (trim(GType) == "VanRaden"  .or.&
      !       trim(GType) == "VanRaden1" .or.&
      !       trim(GType) == "VanRaden2" .or.&
      !       trim(GType) == "Yang") then
      !     do j=1,nLoc
      !       do i=1,nAnisG
      !         if ((Genos(i,j)>=0.0).and.(Genos(i,j)<=2.0)) then
      !           ZMat(i,j)=Genos(i,j)-2.0d0*AlleleFreq(j)
      !         else
      !           ZMat(i,j)=0.0d0
      !         end if
      !       end do
      !     end do
      !   end if
      !   if (trim(GType) == "Nejati-Javaremi" .or.&
      !       trim(GType) == "Day-Williams") then
      !     do j=1,nLoc
      !       do i=1,nAnisG
      !         if ((Genos(i,j)>=0.0).and.(Genos(i,j)<=2.0)) then
      !           ZMat(i,j)=Genos(i,j)-1.0d0
      !         else
      !           ZMat(i,j)=0.d00 ! @todo: is this OK?
      !         end if
      !       end do
      !     end do
      !   end if
      !   ! Scale centered allele dosages
      !   if (trim(GType) == "VanRaden2" .or. trim(GType) == "Yang") then
      !     do j=1,nLoc
      !       Tmp=2.0d0*AlleleFreq(j)*(1.0d0-AlleleFreq(j))
      !       if (Tmp > tiny(Tmp)) then
      !         ZMat(:,j)=ZMat(:,j)/sqrt(Tmp)
      !       end if
      !     end do
      !   end if

      !   ! Z'
      !   tZMat=transpose(ZMat)

      !   WhichMat=0
      !   do j=1,nTrait
      !     do i=j,nTrait
      !       WhichMat=WhichMat+1

      !       ! ZHZ'
      !       if (LocusWeightPresent) then
      !         DMatSum=0.0d0
      !         do k=1,nLoc
      !           DMat(k,k)=sqrt(LocusWeight(k,i))*sqrt(LocusWeight(k,j))
      !           DMatSum=DMatSum+DMat(k,k)
      !         end do
      !         ! @todo: use DGEMM equivalent for X * Diagonal
      !         TmpZMat=matmul(ZMat,DMat)
      !         ! @todo: use DGEMM
      !         GMat(:,:,WhichMat)=matmul(TmpZMat,tZMat)
      !       else
      !         ! @todo: use DGEMM
      !         GMat(:,:,WhichMat)=matmul(ZMat,tZMat)
      !       end if

      !       ! ZHZ'/Denom
      !       if (trim(GType) == "VanRaden" .or. trim(GType) == "VanRaden1") then
      !         GMat(:,:,WhichMat)=GMat(:,:,WhichMat)/(2.0d0*sum(AlleleFreq(:)*(1.0d0-AlleleFreq(:))))
      !       end if
      !       if (trim(GType) == "VanRaden2" .or.&
      !           trim(GType) == "Yang"      .or.&
      !           trim(GType) == "Nejati-Javaremi") then
      !         GMat(:,:,WhichMat)=GMat(:,:,WhichMat)/nLocD
      !       end if

      !       ! Put back scale from [-1,1] to [0,2]
      !       if (trim(GType) == "Nejati-Javaremi") then
      !         if (LocusWeightPresent) then
      !           Tmp=DMatSum/nLocD
      !         else
      !           Tmp=1.0d0
      !         end if
      !         GMat(:,:,WhichMat)=GMat(:,:,WhichMat)+Tmp
      !       end if

      !       ! @todo: needs testing (was getting some weird values)
      !       ! if (trim(GType) == "Day-Williams") then
      !       !   Tmp=0.0d0
      !       !   do k=1,nLoc
      !       !     Tmp=Tmp + AlleleFreq(k)*AlleleFreq(k) + (1.0d0-AlleleFreq(k))*(1.0d0-AlleleFreq(k))
      !       !   end do
      !       !   do k=1,nAnisG
      !       !     do l=1,nAnisG
      !       !       ! @todo: could do just lower triangle, but would have to jump around in memory, i.e., G(j,i)=G(i,j)
      !       !       !       which is faster?
      !       !       ! GMat(l,k,WhichMat)+nLoc is the total number of (observed) IBS matches, i.e., 2*e(i,j) in Day-Williams
      !       !       ! Multiply and divide by 2, because we are building covariance matrix instead of probability matrix
      !       !       GMat(l,k,WhichMat)=2.0d0*((GMat(l,k,WhichMat)+nLocD)/2.0d0-Tmp)/(nLocD-Tmp)
      !       !     end do
      !       !   end do
      !       ! end if

      !       ! Different diagonal for Yang altogether
      !       if (trim(GType) == "Yang") then
      !         do l=1,nAnisG
      !           GMat(l,l,WhichMat)=0.0d0
      !         end do
      !         do k=1,nLoc
      !           Tmp=2.0d0*AlleleFreq(k)*(1.0d0-AlleleFreq(k))
      !           if (Tmp > tiny(Tmp)) then
      !             if (LocusWeightPresent) then
      !               Tmp2=sqrt(LocusWeight(k,i))*sqrt(LocusWeight(k,j))
      !             else
      !               Tmp2=1.0d0
      !             end if
      !             do l=1,nAnisG
      !               Tmp3=Tmp2 * (1.0d0 + ((Genos(l,k)*Genos(l,k) - (1.0d0+2.0d0*AlleleFreq(k))*Genos(l,k) + 2.0d0*AlleleFreq(k)*AlleleFreq(k)) / Tmp))/nLocD
      !               GMat(l,l,WhichMat)=GMat(l,l,WhichMat)+Tmp3
      !             end do
      !           end if
      !         end do
      !       end if

      !       ! Fudge diagonal
      !       do l=1,nAnisG
      !         GMat(l,l,WhichMat)=GMat(l,l,WhichMat)+DiagFudge
      !       end do

      !       ! Export etc.
      !       if (GFullMat) then
      !         write(filout,'("GFullMatrix"i0,"-"i0".txt")')i,j
      !         write(nChar,*) nAnisG
      !         fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
      !         open(unit=202,file=trim(filout),status="unknown")
      !         do m=1,nAnisG
      !           write(202,fmt) IdGeno(m),GMat(:,m,WhichMat)
      !         end do
      !         close(202)
      !       end if

      !       if (GIJA) then
      !         fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
      !         write(filout,'("Gija"i0,"-"i0".txt")')i,j
      !         open(unit=202,file=trim(filout),status="unknown")
      !         do m=1,nAnisG
      !           do n=m,nAnisG
      !             ! No test for non-zero here as all elements are non-zero
      !             write(202,fmt) IdGeno(n),IdGeno(m),GMat(n,m,WhichMat)
      !           end do
      !         end do
      !         close(202)
      !       end if

      !       if (MakeInvG) then
      !         allocate(InvGMat(nAnisG,nAnisG,nGMat))

      !         print*, "Start inverting G - ", trim(GType)
      !         InvGMat(:,:,WhichMat)=GMat(:,:,WhichMat)
      !         call invert(InvGMat(:,:,WhichMat),nAnisG,.true., 1)

      !         print*, "Finished inverting G - ", trim(GType)

      !         if (InvGFullMat) then
      !           write(nChar,*) nAnisG
      !           fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
      !           write(filout,'("InvGFullMatrix"i0,"-"i0".txt")')i,j
      !           open(unit=202,file=trim(filout),status="unknown")
      !           do m=1,nAnisG
      !             write(202,fmt) IdGeno(m),InvGMat(:,m,WhichMat)
      !           end do
      !           close(202)
      !         end if

      !         if (InvGIJA) then
      !           write(filout,'("InvGija"i0,"-"i0".txt")')i,j
      !           fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
      !           open(unit=202,file=trim(filout),status="unknown")
      !           do m=1,nAnisG
      !             do n=m,nAnisG
      !               write(202,fmt) IdGeno(n),IdGeno(m),InvGMat(n,m,WhichMat)
      !             end do
      !           end do
      !           close(202)
      !         end if
      !       end if

      !     end do
      !   end do
      !   deallocate(tZMat)
      !   deallocate(TmpZMat)
      !   print*, "Finished making G - ", trim(GType)
      ! end subroutine

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

      !         if (HFullMat) then
      !           close(202)
      !         end if
      !         if (HIJA) then
      !           close(204)
      !         end if

      !         print *, "End writing H matrices"

      !       end if

      !       if (MakeInvH) then
      !         print *, "Start inverting scaled G matrix"
      !         call invert(G22, size(G22, 1), .true., 1)

      !         !print *, "Gw inverted"
      !         !write(fmt2, "(i0)") size(G22,1)
      !         !fmt1="(a8,"//trim(adjustl(fmt2))//"f8.4)"
      !         !do i=1,size(G22,1)
      !         ! write(*,fmt1) IdGeno(i), G22(i,:)
      !         !end do

      !         !print *, "A22 inverted"
      !         !do i=1,size(G22,1)
      !         ! write(*,fmt1) IdGeno(i), InvA22(i,:)
      !         !end do

      !         !print *, "InvA(22)"
      !         !do i=1,size(G22,1)
      !         ! j = i+10
      !         ! write(*,fmt1) Ids(j), InvAMat(j,11:25)
      !         !end do

      !         print *, "End inverting scaled G matrix"

      !         print *, "Start writing inverted H matrices (full and/or ija)"

      !         if (InvHFullMat) then
      !           write(filout,'("InvHFullMatrix"i0,"-"i0".txt")') t1,t2
      !           write(nChar,*) nAnisH
      !           fmt1="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
      !           open(unit=202,file=trim(filout),status="unknown")
      !         end if

      !         if (InvHIJA) then
      !           write(filout,'("InvHija"i0,"-"i0".txt")') t1,t2
      !           fmt2="(a,' ',a,' ',"//trim(adjustl(OutputFormat))//")"
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

      !         if (InvHFullMat) then
      !           close(202)
      !         end if
      !         if (InvHIJA) then
      !           close(204)
      !         end if
      !         print *, "End writing inverted H matrices (full and ija)"

      !       end if

      !       deallocate(Hrow)
      !       deallocate(G22)
      !     end do
      !   end do
      !   deallocate(Ids)
      ! end subroutine

      !###########################################################################

      ! subroutine invert(x,n,sym, method)

      !   ! Interface to call inverse subroutines from BLAS/LAPACK libraries

      !   ! x symmetric positive-definite matrix to be inverted
      !   ! n matrix dimension
      !   ! sym return lower-triangular (sym=.false) or full matrix (sym=.true.)
      !   ! method for inversion
      !   ! 0 -- Generalised solving using LU decomposition (dgetrs)
      !   ! 1 -- Cholesky decomposition

      !   implicit none
      !   integer(int32), intent(in) :: n,method
      !   integer(int32) :: i,j,info

      !   real(real64),intent(inout) :: x(n,n)
      !   real(real64),allocatable :: Iden(:,:)

      !   logical, intent(in) :: sym

      !   if (method == 0) then
      !     !Solves a general system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by SGETRF/CGETRF
      !     !http://physics.oregonstate.edu/~rubin/nacphy/lapack/routines/dgetrs.html

      !     allocate(Iden(n,n))
      !     ForAll(i = 1:n, j = 1:n) Iden(i,j) = (i/j)*(j/i)  !https://rosettacode.org/wiki/Identity_matrix#Notorious_trick

      !     !https://software.intel.com/en-us/node/468712
      !     !Solves a system of linear equations with an LU-factored square coefficient matrix, with multiple right-hand sides.
      !     ! dgetrs(trans,n,nrhs,A,b,lda,ldb,info)
      !     !Output: Solution overwrites `b`.
      !     call dgetrs("N",n,n,x,Iden,n,n,info)
      !     if (info /= 0) then
      !       print *, "Matrix not positive-definite - info",info
      !       stop 1
      !     end if

      !     x(:,:) = Iden(:,:)

      !   else if (method == 1) then

      !     ! Computes the Cholesky factorization of a symmetric positive definite matrix
      !     ! https://software.intel.com/en-us/node/468690
      !     call dpotrf("L",n,x,n,info)
      !     if (info /= 0) then
      !       print*,"Matrix not positive-definite - info",info
      !       stop 1
      !     end if

      !     ! Computes the inverse of a symmetric positive definite matrix,
      !     !   using the Cholesky factorization computed by dpotrf()
      !     ! https://software.intel.com/en-us/node/468824
      !     call dpotri("L",n,x,n,info)
      !     if (info /= 0) then
      !      print*,"Matrix not positive-definite - info",info
      !      stop 1
      !     end if

      !     ! Fills the upper triangle
      !     if (sym) then
      !       forall (i=1:n,j=1:n,j>i) x(i,j)=x(j,i)
      !     end if

      !   end if
      ! end subroutine

    !###########################################################################

end module

!###############################################################################
