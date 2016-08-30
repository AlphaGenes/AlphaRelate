#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#define BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX
#define DASH "/"
#define COPY "cp"
#define MKDIR "mkdir -p"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#else
#define DASH "\"
#define COPY "copy"
#define MKDIR "md"
#define RMDIR "rmdir /S"
#define RM "del"
#define RENAME "move /Y"
#endif

!TODO: write output subroutines and call them in different places

!#########################################################################

module AlphaAHGModule

  use ISO_Fortran_env
  use AlphaHouseMod, only : CountLines

  INTEGER(int32),PARAMETER :: LENGAN=20

  integer(int32) :: nTrait,nSnp,nAnisG,nAnisP,nAnisRawPedigree,AllFreqSelCycle,nCols,nAnisH
  integer(int32) :: nGMats, GlobalExtraAnimals, OldAmatNInd
  integer(int32) :: NRMmem, shell, shellmax, shellWarning
  integer(int32),allocatable :: MapAnimal(:),seqid(:),seqsire(:),seqdam(:),seqoutput(:)
  integer(int32),allocatable :: RecodeGenotypeId(:),passedorder(:),RecPed(:,:),dooutput(:)
  integer(int32),allocatable :: OldAmatId(:)

  real(real64) :: AlleleFreqAll,DiagFudge,ScaleGToA
  real(real64),allocatable :: AlleleFreq(:),Adiag(:)
  real(real64),allocatable :: Weights(:,:),Zmat(:,:),tpose(:,:),Genos(:,:),Amat(:,:),InvAmat(:,:)
  real(real64),allocatable :: WeightStand(:,:,:),Gmat(:,:,:),InvGmat(:,:,:)

  character(len=20) :: GType,OutputFormat
  character(len=1000) :: GenotypeFile,AlleleFreqFile,OldAmatFile
  character(len=LENGAN),allocatable :: ped(:,:),Id(:),sire(:),dam(:),IdGeno(:)

  logical :: PedigreePresent,WeightsPresent,ScaleGByRegression, OldAmatPresent
  logical :: Gmake, GInvMake, AMake, AInvMake, HMake, HInvMake
  logical :: GFullMat, GIJA, IGFullMat, IGIJA, AFullMat, AIJA, IAFullMat, IAIJA, HFullMat, HIJA, IHFullMat, IHIJA
  logical,allocatable :: MapToG(:), AnimalsInBoth(:)

  contains

  !#############################################################################

    subroutine ReadParam
      implicit none

      integer(int32) :: i,j,n,OutputPositions,OutputDigits

      character(len=1000) :: dumC,option,PedigreeFile,WeightFile
      character(len=200) :: OutputPositionsC,OutputDigitsC

      open(unit=11,file="AlphaAGHSpec.txt",status="old")

      ! Input parameters
      read(11,*) dumC
      read(11,*) dumC,GenotypeFile
      read(11,*) dumC,PedigreeFile
      read(11,*) dumC,WeightFile
      read(11,*) dumC,AlleleFreqFile
      if (trim(AlleleFreqFile) == "Fixed") then
        backspace(11)
        read(11,*) dumC,dumC,AlleleFreqAll
      end if
      read(11,*) dumC,nTrait
      read(11,*) dumC,nSnp
      read(11,*) dumC,DiagFudge
      read(11,*) dumC,GType
      if (trim(GenotypeFile) /= "None"     .and.&
          trim(GType) /= "VanRaden"        .and.&
          trim(GType) /= "VanRaden1"       .and.&
          trim(GType) /= "VanRaden2"       .and.&
          trim(GType) /= "Yang"            .and.&
          trim(GType) /= "Nejati-Javaremi" .and.&
          trim(GType) /= "Day-Williams") then
        print*, "TypeG must be either VanRaden=VanRaden1, VanRaden2, Yang, Nejati-Javaremi, or Day-Williams"
        print*, GType
        stop 1
      end if
      read(11,*) dumC,option
      if (trim(option) == 'Regression') then
        ScaleGByRegression = .true.
      else
        ScaleGByRegression = .false.
        read(option,*) ScaleGtoA
      end if

      ! Output options
      read(11,*) dumC

      read(11,*) dumC,OutputPositions,OutputDigits
      write(OutputPositionsC,*) OutputPositions
      write(OutputDigitsC,*)    OutputDigits
      OutputFormat="f"//trim(adjustl(OutputPositionsC))//"."//trim(adjustl(OutputDigitsC))

      ! Make G matrix?
      GFullMat = .false.
      GIJA = .false.
      GMake = .false.
      read(11,*) dumC, option
      GFullMat = trim(option) == 'Yes'
      if (GFullMat) then
        GMake = .true.
      end if

      read(11,*) dumC, option
      GIJA = trim(option) == 'Yes'
      if (GIJA) then
        GMake = .true.
      end if

      ! Make inverted G matrix ?
      IGFullMat = .false.
      IGIJA = .false.
      GInvMake = .false.
      read(11,*) dumC, option
      IGFullMat = trim(option) == 'Yes'
      if (IGFullMat) then
        GInvMake = .true.
      end if

      read(11,*) dumC, option
      IGIJA = trim(option) == 'Yes'
      if (IGIJA) then
        GInvMake = .true.
      end if

      ! Make A matrix?
      AFullMat = .false.
      AIJA = .false.
      AMake = .false.
      read(11,*) dumC, option
      AFullMat = trim(option) == 'Yes'
      if (AFullMat) then
        AMake = .true.
      end if

      read(11,*) dumC, option
      AIJA = trim(option) == 'Yes'
      if (AIJA) then
        AMake = .true.
      end if

      ! Make inverted A matrix?
      IAFullMat = .false.
      IAIJA = .false.
      AInvMake = .false.
      read(11,*) dumC, option
      IAFullMat = trim(option) == 'Yes'
      if (IAFullMat) then
        AInvMake = .true.
      end if

      read(11,*) dumC, option
      IAIJA = trim(option) == 'Yes'
      if (IAIJA) then
        AInvMake = .true.
      end if

      ! Make H matrix?
      HFullMat = .false.
      HIJA = .false.
      HMake = .false.
      read(11,*) dumC, option
      HFullMat = trim(option) == 'Yes'
      if (HFullMat) then
        HMake = .true.
      end if

      read(11,*) dumC, option
      HIJA = trim(option) == 'Yes'
      if (HIJA) then
        HMake = .true.
      end if

      ! Make inverted H matrix?
      IHFullMat = .false.
      IHIJA = .false.
      HInvMake = .false.
      read(11,*) dumC, option
      IHFullMat = trim(option) == 'Yes'
      if (IHFullMat) then
        HInvMake = .true.
      end if

      read(11,*) dumC, option
      IHIJA = trim(option) == 'Yes'
      if (IHIJA) then
        HInvMake = .true.
      end if

      if (HMake) then
        GMake = .true.
        AMake = .true.
      end if

      if (HInvMake) then
        GMake = .true.
        AMake = .true.
        AInvMake = .true.
      end if

      OldAmatPresent = .false.
      n=CountLines("AlphaAGHSpec.txt")
      if (n > 24) then
        print *, "This is experimental feature/hack = not well tested and might be removed !!!"
        print *, "  - It requires id of individuals to be numeric and sequential and no unknown parents"
        print *, "  - It requires the old A matrix between the parents of individuals whose A matrix will be built"
        print *, "  - It switches off creation of other matrices (exit after Amat is done)"
        read(11,*) dumC, OldAmatFile, OldAmatNInd
        OldAmatPresent = .true.
        open(unit=103, file=OldAmatFile, status="unknown")
      end if

      allocate(AlleleFreq(nSnp))
      allocate(WeightStand(nSnp,nTrait,nTrait))
      allocate(Weights(nSnp,nTrait))

      if (trim(GenotypeFile)/='None') then
        open(unit=101,file=trim(GenotypeFile),status="old")
      end if

      PedigreePresent=.false.
      if (trim(PedigreeFile)/='None') then
        open(unit=102,file=trim(PedigreeFile),status="old")
        PedigreePresent=.true.
      end if
      if (trim(WeightFile)/='None') then
        open(unit=103,file=trim(WeightFile),status="old")
        WeightsPresent = .true.
      else
        WeightsPresent = .false.
      end if

      nGMats=0
      do i=1,nTrait
        do j=i,nTrait
          nGMats=nGMats+1
        end do
      end do

      if (trim(GenotypeFile)=='None') then
        if (GMake .or. GInvMake .or. HMake .or. HInvMake) then
          print *, 'In order to create G or H matrices, a genotype file must be given.'
        end if
        GMake=.false.
        GInvMake=.false.
        HMake=.false.
        HInvMake=.false.
      end if

      if (trim(PedigreeFile)=='None') then
        if (AMake .or. AInvMake .or. HMake .or. HInvMake) then
          print *, 'In order to create A or H matrices, a pedigree must be given.'
        end if
        AMake=.false.
        AInvMake=.false.
        HMake=.false.
        HInvMake=.false.
      end if

      if (ScaleGByRegression) then
        Amake = .true.
      end if
    end subroutine

  !#############################################################################

    subroutine ReadData
      implicit none

      integer(int32) :: i,j,stat,GenoInPed,fourthColumn,nMissing

      character(len=1000) :: dumC
      character(len=LENGAN), dimension(1:3) :: pedline

      call CountInData

      allocate(Genos(nAnisG,nSnp))
      allocate(Zmat(nAnisG,nSnp))
      allocate(IdGeno(nAnisG))
      !allocate(RecodeIdGeno(nAnisG))

      if (.not. PedigreePresent) then
        allocate(RecPed(0:nAnisG,4))
        nAnisP=nAnisG
        RecPed(:,:)=0
        RecPed(:,4)=1
        do i=1,nAnisP
          RecPed(i,1)=i
        end do
      else
        !!! Attempt to magically detect whether there are three or four columns:
        read(102, '(a)', iostat=stat) dumC
        if (stat /= 0) then
          print *, 'Problems reading Pedigree file.'
          !print *, stat, dumC
          stop 1
        end if
        rewind(102)

        ! Test if the line contains three or four columns
        fourthColumn = -99
        read(dumC, *, iostat=stat) pedline, fourthColumn
        if (stat .eq. -1 .or. fourthColumn .eq. -99) then
          nCols = 3
        else
          nCols = 4
        end if
        !print *, nCols, pedline, fourthColumn
        !!! We now know whether there are three or four columns.

        allocate(Ped(nAnisRawPedigree,4))
        Ped(:,4) = '1'
        do i=1,nAnisRawPedigree
          read(102,*) ped(i,1:nCols)
        end do
        call PVseq(nAnisRawPedigree,nAnisP,0)

        allocate(RecPed(0:nAnisP,4))

        RecPed(0,:)=0
        RecPed(:,4)=1
        do i=1,nAnisP
          RecPed(i,1)=i
        end do
        RecPed(1:nAnisP,2)=seqsire(1:nAnisP)
        RecPed(1:nAnisP,3)=seqdam(1:nAnisP)
        RecPed(1:nAnisP,4)=dooutput(1:nAnisP)
      end if

      do i=1,nAnisG
        read(101,*) IdGeno(i),Genos(i,:)
      end do

      ! These three vectors uses the Pedigree animals as base,
      ! i.e. after reordering, the indices for the nth pedigree animal is n.
      allocate(MapAnimal(1:(nAnisP+nAnisG)))
      allocate(MapToG(1:(nAnisP+nAnisG)))
      allocate(AnimalsInBoth(1:nAnisP+nAnisG))
      MapAnimal = 0
      MapToG = .false.
      AnimalsInBoth = .false.
      nAnisH = nAnisP
      do i=1,nAnisP
        MapAnimal(i) = i
      end do

      AnimalsInBoth = .false.
      if (PedigreePresent) then
        !! Match Genotyped animals to pedigree animals
        do i=1,nAnisG
          GenoInPed=0
          do j=1,nAnisP
            if (trim(IdGeno(i))==trim(Id(j))) then
              MapToG(j) = .true.
              MapAnimal(j) = i
              AnimalsInBoth(j) = .true.
              GenoInPed=1
              exit
            end if
          end do
          if (GenoInPed==0) then
            nAnisH = nAnisH + 1
            MapAnimal(nAnisH) = i
            MapToG(nAnisH) = .true.
            print*, "Genotyped individual not in pedigree file - ",trim(IdGeno(i))
            !stop 1
          end if
        end do
      end if

      if (WeightsPresent) then
        do i=1,nSnp
          read(103,*) dumC,Weights(i,:)
        end do
      else
        Weights(:,:)=1
      end if

      !! Allele frequencies
      if (trim(GenotypeFile) /= "None") then
        if (trim(AlleleFreqFile) == 'None') then
          !Calculate Allele Freq
          AlleleFreq(:)=0.0d0
          do j=1,nSnp
            nMissing=0
            do i=1,nAnisG
              if ((Genos(i,j)>-0.1).and.(Genos(i,j)<2.1)) then
                AlleleFreq(j)=AlleleFreq(j)+Genos(i,j)
              else
                nMissing=nMissing+1
              end if
            end do
            ! Write the frequency of SNP j in array. If all SNPs are missing, then freq_j=0
            if (nAnisG>nMissing) then
              AlleleFreq(j)=AlleleFreq(j)/(2.0d0*dble((nAnisG-nMissing)))
            else
              AlleleFreq(j)=0.0d0
            end if
          end do
          open(unit=201, file='AlleleFreqTest.txt',status='unknown')
          do j=1,nSnp
            write(201,*) j,AlleleFreq(j)
          end do
          close(202)
        else
          if (trim(AlleleFreqFile) == "Fixed") then
            AlleleFreq(:) = AlleleFreqAll
            open(unit=201, file='AlleleFreqTest.txt',status='unknown')
            do j=1,nSnp
              write(201,*) j,AlleleFreq(j)
            end do
            close(202)
          else
            ! Read allele frequencies from file.
            open(unit=202, file=trim(AlleleFreqFile), status='OLD')
            do i=1,nSnp
              read(202,*,iostat=stat) dumC,AlleleFreq(i)  !AlleleFrequencies are kept in second column to keep consistency with AlphaSim.
              if (stat /= 0) then
                print*,"Problems reading allele frequency file."
                stop 1
              end if
            end do
            close(202)
          end if
        end if
      end if
    end subroutine

  !#############################################################################

    subroutine MakeInvAMatrix
      implicit none

      integer(int32) :: i,m,n,s,FId,MId

      real(real64) :: Inbreeding(0:nAnisP),Dii(nAnisP),InvDii

      character(len=1000) :: nChar,fmt

      logical :: AnimToWrite(nAnisP),FIdL,MIdL

      allocate(InvAmat(0:nAnisP,0:nAnisP))

      print*, "Start calculating inbreeding coefficients"
      call dinbreeding(RecPed(0:nAnisP,1),RecPed(0:nAnisP,2),RecPed(0:nAnisP,3),Inbreeding,nAnisP)
      open(unit=202,file="PedigreeBasedInbreeding.txt",status="unknown")
      print*, "End calculating inbreeding coefficients"
      do i=1,nAnisP
        write(202,'(a20,20000f10.5)') Id(i),Inbreeding(i)
      end do
      close(202)

      print*, "Start making A inverse"
      InvAmat=0.0d0
      do i=1,nAnisP
        FId=RecPed(i,2)
        FIdL=FId/=0
        MId=RecPed(i,3)
        MIdL=MId/=0
        ! Variance of founder effects and Mendelian sampling terms
        Dii(i)=1.0d0
        if (FIdL) then
          Dii(i)=Dii(i)-0.25d0*(1.0d0+Inbreeding(FId))
        end if
        if (MIdL) then
          Dii(i)=Dii(i)-0.25d0*(1.0d0+Inbreeding(MId))
        end if
        ! Precision for the individual
        InvDii=1.0d0/Dii(i)
        InvAmat(i,i)=InvDii
        ! Add precision to the father and set the co-precision
        if (FIdL) then
          InvAmat(FId,FId)=InvAmat(FId,FId)+InvDii/4.0d0
          InvAmat(i,FId)=InvAmat(i,FId)-InvDii/2.0d0
          InvAmat(FId,i)=InvAmat(i,FId)
        end if
        ! Add precision to the mother and set the co-precision
        if (MIdL) then
          InvAmat(MId,MId)=InvAmat(MId,MId)+InvDii/4.0d0
          InvAmat(i,MId)=InvAmat(i,MId)-InvDii/2.0d0
          InvAmat(MId,i)=InvAmat(i,MId)
        end if
        ! Add co-precision between the father and mother
        if (FIdL .and. MIdL) then
          InvAmat(FId,MId)=InvAmat(FId,MId)+InvDii/4.0d0
          InvAmat(MId,FId)=InvAmat(FId,MId)
        end if
      end do
      print*, "Finished making A inverse"

      if (IAFullMat) then
        AnimToWrite = RecPed(1:nAnisP,4) == 1
        s = count(AnimToWrite)
        write(*,'(a40,i6,a11)') " Start writing A inverse full matrix for", s," individuals"
        write(nChar,*) s
        fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="InvAFullMatrix.txt",status="unknown")
        do m=1,nAnisP
          if (AnimToWrite(m)) then
            write(202,fmt) Id(m), pack(InvAmat(1:nAnisP,m), AnimToWrite)
          end if
        end do
        close(202)
        print*, "End writing A inverse full matrix"
      end if

      if (IAIJA) then
        print*, "Start writing A inverse ija"
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="InvAija.txt",status="unknown")
        do m=1,nAnisP
          do n=m,nAnisP
            if (InvAmat(n,m) /= 0.0d0) then
              write(202,fmt) Id(n),Id(m),InvAmat(n,m)
            end if
          end do
        end do
        close(202)
        print*, "End writing A inverse ija"
      end if
    end subroutine

  !#############################################################################

    subroutine MakeAMatrix
      implicit none

      integer(int32) :: i,j,k,l,m,n,s,d,div,MinId,MaxId,Start,Endin

      real(real64) :: AmatAvg

      character(len=1000) :: nChar,fmt

      logical :: AnimToWrite(nAnisP)

      print*, "Start making A"
      if (OldAmatPresent) then
        allocate(OldAmatId(OldAmatNInd))
        do j = 1, OldAmatNInd
          read(103, *) OldAmatId(j)
        end do
        rewind(103)
        MinId = minval(OldAmatId)
        MaxId = maxval(OldAmatId)
        allocate(Amat(1:(OldAmatNInd+nAnisP-MaxId),&
                      1:(OldAmatNInd+nAnisP-MaxId)))
        Amat = 0.0d0
        do j = 1, OldAmatNInd
          read(103, *) OldAmatId(j), Amat(1:OldAmatNInd,j)
          if (j > 1) then
            if (.not.(OldAmatId(j) > OldAmatId(j-1))) then
              print *, "Id are not sequential!"
              stop 1
            end if
          end if
        end do
        k = OldAmatNInd
        do i=MaxId+1,nAnisP
            k = k + 1
            s = RecPed(i,2) - MinId + 1
            d = RecPed(i,3) - MinId + 1
            l = k
            do j=1,k-1
                Amat(j,k)=(Amat(j,s)+Amat(j,d))/2.0d0
                Amat(k,j)=Amat(j,k)
                !print *,i,k,j,s,d,Amat(j,s),Amat(j,d),Amat(j,k)
            end do
            Amat(k,k)=1.0d0+Amat(s,d)/2.0d0
            !print *,i,k,s,d,Amat(s,d),Amat(k,k)
        end do
        RecPed(1:nAnisP,4) = 0
        RecPed((MaxId+1):nAnisP,4) = 1
      else
        allocate(Amat(0:nAnisP,0:nAnisP))
        Amat=0.0d0
        do i=1,nAnisP
            do j=1,i-1
                Amat(j,i)=(Amat(j,RecPed(i,2))+Amat(j,RecPed(i,3)))/2.0d0
                Amat(i,j)=Amat(j,i)
            end do
            Amat(i,i)=1.0d0+Amat(RecPed(i,2),RecPed(i,3))/2.0d0
        end do
      end if
      print*, "Finished making A"

      if (AFullMat) then
        AnimToWrite = RecPed(1:nAnisP,4) == 1
        s = count(AnimToWrite)
        write(*,'(a32,i6,a11)') " Start writing A full matrix for", s," individuals"
        write(nChar,*) s
        fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="AFullMatrix.txt",status="unknown")
        if (.not.OldAmatPresent) then
          do m=1,nAnisP
            if (AnimToWrite(m)) then
              write(202,fmt) Id(m), pack(Amat(1:nAnisP,m), AnimToWrite)
            end if
          end do
        else
          Start = OldAmatNInd+1
          Endin = size(Amat,1)
          do m=Start,Endin
            !write(*,fmt)   Id(m+MinId-1), Amat(Start:Endin,m)
            write(202,fmt) Id(m+MinId-1), Amat(Start:Endin,m)
          end do
        end if
        close(202)
        print*, "End writing A full matrix"
      end if

      if (OldAmatPresent) then
        stop
      end if

      if (AIJA) then
        ! TODO: AnimToWrite is not being used here
        write(*,'(a24,i6,a11)') " Start writing A ija for", s," individuals"
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="Aija.txt",status="unknown")
        do m=1,nAnisP
          do n=m,nAnisP
            if (Amat(n,m) > 0.0d0)  then
              write(202,fmt) Id(n),Id(m),Amat(n,m)
            end if
          end do
        end do
        close(202)
        print*, "End writing A ija"
      end if

      ! Record diagonals of animals in both A and G:
      if ((Hmake .or. HInvMake) .and. ScaleGByRegression) then
        n = Count(AnimalsInBoth)
        allocate(Adiag(0:n))
        div = dble(n**2)
        AmatAvg = 0.0d0
        k = 0
        do i = 1,nAnisP
          if (.not. AnimalsInBoth(i)) then
            cycle
          end if
          k = k + 1
          Adiag(k) = Amat(i,i)
          do j=1,nAnisP
            if (AnimalsInBoth(j)) then
              AmatAvg=AmatAvg + Amat(j,i) * 2.0d0 / div
            end if
          end do
        end do
        Adiag(0) = AmatAvg
      end if
    end subroutine

  !#############################################################################

    subroutine MakeG
      implicit none

      integer(int32) :: i,j,k,l,m,n,WhichMat

      real(real64) :: Tmp,TmpWt(nSnp),TmpVal(1),Denom

      character(len=1000) :: filout,nChar,fmt

      allocate(Gmat(nAnisG,nAnisG,nGMats))
      allocate(tpose(nSnp,nAnisG))

      print*, "Start making G - ", trim(GType)

      !Standardise weights
      ! TODO: optimise use of indices to get good performance, perhaps these Weights
      !       should be ditched altogether
      if (trim(GType) == "VanRaden"  .or.&
          trim(GType) == "VanRaden1") then
        do i=1,nTrait
          do j=1,nTrait
            if (i==j) then
              TmpVal(1)=sum(Weights(:,j))
              do k=1,nSnp
                WeightStand(k,i,j)=Weights(k,j)/TmpVal(1)
              end do
            else
              do k=1,nSnp
                TmpWt(k)=Weights(k,i)*Weights(k,j)
              end do
              TmpVal(1)=sum(TmpWt(:))
              do k=1,nSnp
                WeightStand(k,i,j)=TmpWt(k)/TmpVal(1)
              end do
            end if
          end do
        end do
      end if
      if (trim(GType) == "VanRaden2"       .or.&
          trim(GType) == "Yang"            .or.&
          trim(GType) == "Nejati-Javaremi" .or.&
          trim(GType) == "Day-Williams") then
        if (WeightsPresent) then
          print*,"Weights not used with option ", trim(GType)
        end if
        WeightStand(:,:,:)=1.0d0
      end if

      !Make Z
      if (trim(GType) == "VanRaden"  .or.&
          trim(GType) == "VanRaden1" .or.&
          trim(GType) == "VanRaden2" .or.&
          trim(GType) == "Yang") then
        do j=1,nSnp
          do i=1,nAnisG
            if ((Genos(i,j)>=0.0).and.(Genos(i,j)<=2.0)) then
              Zmat(i,j)=Genos(i,j)-2.0d0*AlleleFreq(j)
            else
              Zmat(i,j)=0.0d0
            end if
          end do
        end do
      end if
      if (trim(GType) == "Nejati-Javaremi" .or.&
          trim(GType) == "Day-Williams") then
        do j=1,nSnp
          do i=1,nAnisG
            if ((Genos(i,j)>=0.0).and.(Genos(i,j)<=2.0)) then
              Zmat(i,j)=Genos(i,j)-1.0d0
            else
              Zmat(i,j)=0.d00 ! TODO: is this correct?
            end if
          end do
        end do
      end if

      if (trim(GType) == "VanRaden2" .or. trim(GType) == "Yang") then
        do j=1,nSnp
          Tmp=AlleleFreq(j)*(1.0d0-AlleleFreq(j))
          if (Tmp > tiny(Tmp)) then
            Zmat(:,j)=Zmat(:,j)/sqrt(2.0d0*Tmp)
          end if
        end do
      end if

      !Make G matrices
      WhichMat=0
      do j=1,nTrait
        do i=j,nTrait
          WhichMat=WhichMat+1

          ! Z'
          tpose=transpose(Zmat)

          ! ZH
          do k=1,nSnp
            Zmat(:,k)=Zmat(:,k)*WeightStand(k,i,j)
          end do

          ! ZHZ'
          GMat(:,:,WhichMat)=matmul(Zmat,tpose)

          ! ZHZ'/Denom
          Denom=1.0d0
          if (trim(GType) == "VanRaden" .or. trim(GType) == "VanRaden1") then
            Denom=0.0d0
            do k=1,nSnp
              Denom=Denom+(AlleleFreq(k)*(1.0d0-AlleleFreq(k))*WeightStand(k,i,j))
            end do
            Denom=2.0d0*Denom
          end if
          if (trim(GType) == "VanRaden2" .or.&
              trim(GType) == "Yang"      .or.&
              trim(GType) == "Nejati-Javaremi") then
            Denom = dble(nSnp)
          end if
          GMat(:,:,WhichMat)=GMat(:,:,WhichMat)/Denom

          if (trim(GType) == "Nejati-Javaremi") then
            GMat(:,:,WhichMat)=GMat(:,:,WhichMat)+1.0d0
          end if

          if (trim(GType) == "Day-Williams") then
            Tmp=0.0d0
            do k=1,nSnp
              Tmp=Tmp + AlleleFreq(k)*AlleleFreq(k) + (1.0d0-AlleleFreq(k))*(1.0d0-AlleleFreq(k))
            end do
            do k=1,nAnisG
              do l=1,nAnisG
                ! GMat(l,k,WhichMat)+nSnp is the observed number of IBS matches, i.e., e(i,j)
                GMat(l,k,WhichMat)=2.0d0*((GMat(l,k,WhichMat)+nSnp)/2.0d0-Tmp)/(nSnp-Tmp)
              end do
            end do
          end if

          ! Different diagonal for Yang altogether
          if (trim(GType) == "Yang") then
            do k=1,nAnisG
              GMat(k,k,WhichMat)=1.0d0
            end do
            do k=1,nSnp
              do l=1,nAnisG
                Tmp=AlleleFreq(k)*(1.0d0-AlleleFreq(k))
                if (Tmp > tiny(Tmp)) then
                  Tmp=(Genos(l,k)*Genos(l,k)-(1.0d0+2.0d0*AlleleFreq(k))*Genos(l,k)+2.0d0*AlleleFreq(k)*AlleleFreq(k))/(2.0d0*Tmp)
                  GMat(l,l,WhichMat)=GMat(l,l,WhichMat)+Tmp/Denom
                end if
              end do
            end do
          end if

          ! Fudge diagonal
          do l=1,nAnisG
            GMat(l,l,WhichMat)=GMat(l,l,WhichMat)+DiagFudge
          end do

          ! Export etc.
          if (GFullMat) then
            write(filout,'("GFullMatrix"i0,"-"i0".txt")')i,j
            write(nChar,*) nAnisG
            fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
            open(unit=202,file=trim(filout),status="unknown")
            do m=1,nAnisG
              write(202,fmt) IdGeno(m),GMat(:,m,WhichMat)
            end do
            close(202)
          end if

          if (GIJA) then
            fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
            write(filout,'("Gija"i0,"-"i0".txt")')i,j
            open(unit=202,file=trim(filout),status="unknown")
            do m=1,nAnisG
              do n=m,nAnisG
                ! No test for non-zero here as all elements are non-zero
                write(202,fmt) IdGeno(n),IdGeno(m),GMat(n,m,WhichMat)
              end do
            end do
            close(202)
          end if

          if (GInvMake) then
            allocate(InvGmat(nAnisG,nAnisG,nGMats))

            print*, "Start inverting G - ", trim(GType)
            InvGmat(:,:,WhichMat)=Gmat(:,:,WhichMat)
            call invert(InvGmat(:,:,WhichMat),nAnisG,.true., 1)

            print*, "Finished inverting G - ", trim(GType)

            if (IGFullMat) then
              write(nChar,*) nAnisG
              fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
              write(filout,'("InvGFullMatrix"i0,"-"i0".txt")')i,j
              open(unit=202,file=trim(filout),status="unknown")
              do m=1,nAnisG
                write(202,fmt) IdGeno(m),InvGMat(:,m,WhichMat)
              end do
              close(202)
            end if

            if (IGIJA) then
              write(filout,'("InvGija"i0,"-"i0".txt")')i,j
              fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
              open(unit=202,file=trim(filout),status="unknown")
              do m=1,nAnisG
                do n=m,nAnisG
                  write(202,fmt) IdGeno(n),IdGeno(m),InvGMat(n,m,WhichMat)
                end do
              end do
              close(202)
            end if
          end if

        end do
      end do
      deallocate(tpose)
      print*, "Finished making G - ", trim(GType)
    end subroutine

  !#############################################################################

    subroutine MakeHAndHInv
      ! Feature added by Stefan Hoj-Edwards, February 2016
      ! Making the Inverse H matrix ala Aguilar et al 201? and Christensen 2012

      ! Prerequisite and assumptions for this subroutine:
      ! There is given both a pedigree and genotype file, and there is an overlap
      ! of animals between two data sets.
      ! Diagonals of from both A have been collected during MakeA and MakeG,
      ! as well as average of A22.
      ! Gmat is already calculated and loaded in memory.
      ! Ainv is calculated an loaded into memory.
      !
      ! Further assumes that animals are ordered the same in both A and G.

      implicit none

      integer(int32) :: i,j,k,m,p,q,div,t1,t2,whichMat,nBoth
      integer(int32),allocatable :: MapToA11(:), MapToA22(:) !Gmap(:),

      real(real64) :: Gmatavg, nom, denom, slope, intercept, Gmean, Amean, Hii
      real(real64),allocatable :: Gdiag(:), Hrow(:), A22(:,:), A22Inv(:,:), G22(:,:), A11(:,:), A12(:,:), tmp(:,:), Gboth(:,:)

      character(len=1000) :: nChar,fmt1, fmt2,filout
      character(len=LENGAN),allocatable :: Ids(:)

      logical,allocatable :: AnimToWrite(:)

      nboth = count(AnimalsInBoth)
      ! Make H and/or Hinv
      allocate(Ids(1:nAnisH))
      allocate(AnimToWrite(1:nAnisH))

      do i=1,nAnisH
        if (MapToG(i)) then
          Ids(i) = IdGeno(MapAnimal(i))
          AnimToWrite(i) = .true.
        else
          Ids(i) = Id(MapAnimal(i))
          AnimToWrite(i) = RecPed(MapAnimal(i),4)
        end if
      end do

      allocate(A22Inv(nBoth,nBoth))
      allocate(MapToA22(nAnisH))
      if (HMake) then
        allocate(A22(nBoth,nBoth))
      end if

      k = 0
      do i=1,nAnisP
        if (.not. AnimalsInBoth(i)) then
          cycle
        end if
        k = k + 1
        MapToA22(i) = k
        m = 0
        do j=1,nAnisP
          if (.not. AnimalsInBoth(j)) then
            cycle
          end if
          m = m + 1
          A22Inv(m,k) = Amat(j,i)
        end do
      end do
      if (HMake) then
        A22 = A22Inv
      end if

      call invert(A22Inv,size(A22Inv,1),.true.,1)

      ! This is the G matrix in Legarra,
      ! Sadly, no genotypes where provided, instead the resulting G matrix was.
      if (.false.) then
        print *, 'Overwriting G matrix with example in Legarra 2008!'
        do i=1,nAnisG
          do j=1,nAnisG
            if (i==j) then
              Gmat(i,j,1) = 1
            else
              Gmat(i,j,1) = 0.7
            end if
          end do
        end do
      end if

      whichMat = 0
      do t1=1,nTrait
        do t2=t1,nTrait
          whichMat = whichMat + 1

          write(*, '(" Starting on H matrix "i0" - "i0)') t1, t2

          ! Collect G22
          allocate(G22(nAnisG,nAnisG))

          G22 = 0.0d0
          do j=1,nAnisG
            do i=1,nAnisG
              nom = Gmat(i,j,whichMat)
              if (i == j) then
                nom = nom - DiagFudge
              end if
              G22(i,j) = nom
            end do
          end do

          if (ScaleGByRegression) then
            allocate(Gdiag(0:nBoth))
            Gdiag=0.0d0
            Gmatavg=0.0d0
            div=dble(nBoth**2)
            !allocate(Gmap(nBoth))

            k = 0
            do i=1,nAnisH
              if (.not. AnimalsInBoth(i)) then
                cycle
              end if
              k = k+1
              Gdiag(k) = G22(MapAnimal(i),MapAnimal(i))
              do j=1,nAnisH
                if (.not. AnimalsInBoth(j)) then
                  cycle
                end if
                Gmatavg=Gmatavg + G22(MapAnimal(j),MapAnimal(i))/div
              end do
            end do
            Gdiag(0) = Gmatavg

            ! Now do simple linear regression
            nom = 0.0d0
            denom = 0.0d0
            Gmean = sum(Gdiag) / dble(size(Gdiag, 1))
            Amean = sum(Adiag) / dble(size(Adiag, 1))
            do i=0,ubound(Adiag, 1)
              nom = nom + (Adiag(i) - Amean) * (Gdiag(i) - Gmean)
              denom = denom + (Adiag(i) - Amean)**2
            end do
            slope = nom / denom
            intercept = Amean - slope * Gmean

            ! Scale G
            G22 = slope * G22 + intercept
            !do i=1,nAnisG
            ! G22(i,i) = G22(i,i) + DiagFudge
            !end do
            print *, 'Scaling of G:'
            write(*, '(a,f7.4,a,f7.4)'), " G* = G x ", slope, " + ", intercept
            deallocate(Gdiag)
          else
            do i=1,nAnisH
              if (.not. MapToG(i)) then
                cycle
              end if
              do j=1,nAnisH
                if (.not. MapToG(j)) then
                  cycle
                end if
                if (AnimalsInBoth(i) .and. AnimalsInBoth(j)) then
                  G22(MapAnimal(j),MapAnimal(i)) = ScaleGToA * G22(MapAnimal(j),MapAnimal(i)) + (1.0d0 - ScaleGToA) * Amat(j,i)
                end if
              end do
            end do
          end if

          do i=1,nAnisG
            G22(i,i) = G22(i,i) + DiagFudge
          end do

          allocate(Hrow(1:count(AnimToWrite)))

          if (HMake) then

            allocate(A11(nAnisP-nBoth, nAnisP-nBoth))
            allocate(A12(nAnisP-nBoth, nBoth))
            allocate(MapToA11(nAnisP))
            allocate(tmp(nAnisP-nBoth, nBoth))
            allocate(Gboth(nBoth,nBoth))

            MapToA11 = 0
            k = 0
            p = 0
            do i=1,nAnisP
              if (AnimalsInBoth(i)) then
                p = p + 1
                q = 0
                do j=1,nAnisP
                  if (.not. AnimalsInBoth(j)) then
                    cycle
                  end if
                  q = q + 1
                  Gboth(q,p) = G22(MapAnimal(j),MapAnimal(i))
                end do
              else
                k = k+1
                m = 0
                MapToA11(i) = k
                do j=1,nAnisP
                  if (AnimalsInBoth(j)) then
                    A12(k,MapAnimal(j)) = Amat(j,i)  !A12 is not symmetrical
                  else
                    m = m+1
                    A11(m,k) = Amat(j,i)
                  end if
                end do
              end if
            end do

            tmp = Matmul(A12, A22Inv)
            !tmp = Matmul(Matmul(tmp, (Gboth - A22)), transpose(tmp))
            tmp = Matmul(tmp, (Gboth-A22))
            tmp = Matmul(tmp, A22Inv)
            !tmp = Matmul(tmp, transpose(A12))

            A11 = A11 + Matmul(tmp, transpose(A12))
            A12 = matmul(matmul(A12, A22Inv), Gboth)

            deallocate(tmp)
            deallocate(Gboth)

            print *, 'Start writing H matrices (full and/or ija)'

            if (HFullMat) then
              write(filout,'("HFullMatrix"i0,"-"i0".txt")') t1,t2
              write(nChar,*) nAnisH
              fmt1="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
              open(unit=202,file=trim(filout),status="unknown")
            end if

            if (HIJA) then
              write(filout,'("Hija"i0,"-"i0".txt")') t1,t2
              fmt2="(a20,a20,"//trim(adjustl(OutputFormat))//")"
              open(unit=204,file=trim(filout),status="unknown")
            end if

            do i=1,nAnisH
              if (AnimToWrite(i) .eq. .false.) then
                cycle
              end if
              Hrow = 0
              k = 0
              do j=1,nAnisH
                if (AnimToWrite(j) .eq. .false.) then
                  cycle
                end if
                k = k + 1
                if (MapToG(i)) then
                  if (MapToG(j)) then
                    Hii = G22(MapAnimal(i),MapAnimal(j))
                  else
                    Hii = A12(MapToA11(j),MapAnimal(i)) ! Remember to transpose
                  end if
                else
                  if (MapToG(j)) then
                    Hii = A12(MapToA11(i),MapAnimal(j))
                  else
                    Hii = A11(MapToA11(i),MapToA11(j))
                  end if
                end if
                if (IHIJA .and. i .le. j .and. Hii /= 0.0d0) then
                  write(204,fmt2) Ids(i), Ids(j), Hii
                end if
                Hrow(k) = Hii
              end do
              if (HFullMat) then
                write(202,fmt1) Ids(i),Hrow(:)
              end if
            end do

            if (HFullMat) then
              close(202)
            end if
            if (HIJA) then
              close(204)
            end if

            print *, 'End writing H matrices'

          end if

          if (HInvMake) then
            print *, 'Start inverting scaled G matrix'
            call invert(G22, size(G22, 1), .true., 1)

            !print *, 'Gw inverted'
            !write(fmt2, '(i0)') size(G22,1)
            !fmt1="(a8,"//trim(adjustl(fmt2))//"f8.4)"
            !do i=1,size(G22,1)
            ! write(*,fmt1) IdGeno(i), G22(i,:)
            !end do

            !print *, 'A22 inverted'
            !do i=1,size(G22,1)
            ! write(*,fmt1) IdGeno(i), A22Inv(i,:)
            !end do

            !print *, 'Ainv(22)'
            !do i=1,size(G22,1)
            ! j = i+10
            ! write(*,fmt1) Ids(j), InvAmat(j,11:25)
            !end do

            print *, 'End inverting scaled G matrix'

            print *, 'Start writing inverted H matrices (full and/or ija)'

            if (IHFullMat) then
              write(filout,'("InvHFullMatrix"i0,"-"i0".txt")') t1,t2
              write(nChar,*) nAnisH
              fmt1="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
              open(unit=202,file=trim(filout),status="unknown")
            end if

            if (IHIJA) then
              write(filout,'("InvHija"i0,"-"i0".txt")') t1,t2
              fmt2="(a,' ',a,' ',"//trim(adjustl(OutputFormat))//")"
              open(unit=204,file=trim(filout),status="unknown")
            end if

            do i=1,nAnisH
              if (AnimToWrite(i) .eq. .false.) then
                cycle
              end if
              Hrow = 0
              k = 0
              do j=1,nAnisH
                if (AnimToWrite(j) .eq. .false.) then
                  cycle
                end if
                k = k + 1
                if (MapToG(i) .and. MapToG(j)) then
                  Hrow(k) = G22(MapAnimal(i),MapAnimal(j))
                  if (i <= nAnisP .and. j <= nAnisP) then
                    Hrow(k) = Hrow(k) + InvAmat(i,j) - A22Inv(MapToA22(i),MapToA22(j))
                  end if
                else if (i <= nAnisP .and. j <= nAnisP) then !if (MapToG(i) .eq. .false. .and. MapToG(j) .eq. .false.  ) then
                  Hrow(k) = InvAmat(i,j)
                end if
                if (IHIJA .and. i .le. j .and. Hrow(k) /= 0.0d0) then
                  write(204,fmt2) trim(Ids(i)), trim(Ids(j)), Hrow(k)
                end if
              end do
              if (IHFullMat) then
                write(202,fmt1) Ids(i),Hrow(:)
              end if
            end do

            if (IHFullMat) then
              close(202)
            end if
            if (IHIJA) then
              close(204)
            end if
            print *, 'End writing inverted H matrices (full and ija)'

          end if

          deallocate(Hrow)
          deallocate(G22)
        end do
      end do
      deallocate(Ids)
    end subroutine

  !#############################################################################

    subroutine CountInData
      implicit none

      integer(int32) :: k

      character(len=300) :: dumC

      nAnisG=0
      do
        read(101,*,iostat=k) dumC
        nAnisG=nAnisG+1
        if (k/=0) then
          nAnisG=nAnisG-1
          exit
        end if
      end do
      rewind(101)
      write(*,'(a2,i6,a33)') "   ",nAnisG," individuals in the genotype file"

      if (PedigreePresent) then
        nAnisRawPedigree=0
        do
          read(102,*,iostat=k) dumC
          nAnisRawPedigree=nAnisRawPedigree+1
          if (k/=0) then
            nAnisRawPedigree=nAnisRawPedigree-1
            exit
          end if
        end do
        rewind(102)
        write(*,'(a2,i6,a33)') "   ",nAnisRawPedigree," individuals in the pedigree file"
      end if
    end subroutine

  !#############################################################################

    subroutine dinbreeding(ianim,isire,idam,f,n)
      !c     inbreeding program from Meuwissen and Luo (1992)
      !c     GSE 24: 305-313
      !ARRAYS START FROM 0.......
      IMPLICIT NONE

      integer(int32) :: n,i,is,id,j,k,ks,kd
      integer(int32) :: ianim(0:n),isire(0:n),idam(0:n),ped(0:n,3),point(0:n)

      real(real64) :: f(0:n),l(n),d(n),fi,r

      point(:) = 0
      l(:) = 0.0d0
      d(:) = 0.0d0

      do i = 0,n
        f(i) = 0.0d0
        if (i.gt.0) then
          ped(i,1) = ianim(i)
          ped(i,2) = isire(i)
          ped(i,3) = idam(i)
        end if
      end do

      f(0) = -1.0d0
      do i = 1, n
        is = isire(i)
        id = idam(i)
        ped(i,2) = max(is,id)
        ped(i,3) = min(is,id)
        d(i) = 0.5d0 - 0.25d0 * (f(is) + f(id))
        if (is.eq.0.or.id.eq.0) then
          f(i) = 0.0d0
        else if ((ped(i-1,2).eq.ped(i,2)).and.(ped(i-1,3).eq.ped(i,3))) then
          f(i) = f(i-1)
        else
          fi = -1.0d0
          l(i) = 1.0d0
          j = i

          do while (j.ne.0)
            k = j
            r = 0.5d0 * l(k)
            ks = ped(k,2)
            kd = ped(k,3)
            if (ks.gt.0) then
              do while (point(k).gt.ks)
                k = point(k)
              end do
              l(ks) = l(ks) + r
              if (ks.ne.point(k)) then
                point(ks) = point(k)
                point(k) = ks
              end if
              if (kd.gt.0) then
                do while (point(k).gt.kd)
                  k = point(k)
                end do
                l(kd) = l(kd) + r
                if (kd.ne.point(k)) then
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
    end subroutine

  !#############################################################################

    subroutine PVseq(nObs,nAnisPedigree,mode)
      implicit none

      integer(int32) :: mode    ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
      integer(int32) :: i, j, newid, itth, itho, ihun, iten, iunit
      integer(int32) :: nsires, ndams, newsires, newdams, nbisexuals, flag
      integer(int32) :: iextra, oldnobs, kn, kb, oldkn, ks, kd
      integer(int32) :: Noffset, Limit, Switch, ihold, ipoint
      integer(int32) :: nObs,nAnisPedigree,verbose
      integer(int32),allocatable :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
      integer(int32),allocatable :: OldN(:), NewN(:), holdsire(:), holddam(:), holdoutput(:)

      character(len=LENGAN)              :: IDhold
      character(len=LENGAN) :: path
      character(len=LENGAN),allocatable :: holdsireid(:), holddamid(:)
      character(len=LENGAN),allocatable :: holdid(:), SortedId(:), SortedSire(:), SortedDam(:)

      allocate(id(0:nobs),sire(nobs),dam(nobs),dooutput(nobs),seqid(nobs),seqsire(nobs),seqdam(nobs),seqoutput(nobs))

      do i=1,nobs
          id(i)=ped(i,1)
          sire(i)=ped(i,2)
          dam(i)=ped(i,3)
          read(Ped(i,4), '(i)') dooutput(i)
      end do


      nAnisPedigree=nObs
      path=".\"

      Verbose=1

      do j = 1, nobs
          If (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') Then
              dam(j) = '0'
              seqdam(j)=0
          end if
          If (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') Then
              sire(j) = '0'
              seqsire(j)=0
          end if
      end do !j

      if(mode.eq.1) then
          !PRINT*,  ' Inserting dummy IDs ... '
          newid=0
          do j = 1, nobs
              if(((sire(j) == '0').and.(dam(j).ne.'0'))  .or. ((sire(j).ne.'0').and.(dam(j) == '0'))) then
                  newid=newid+1
                  if(newid.gt.99999) then
                      !         PRINT*, newid, ' ...'
                      print*,'too many dummy single parent IDs'
                      stop 1
                  end if
                  itth=int(newid/10000)
                  itho=int(newid/1000)-10*itth
                  ihun=int(newid/100)-10*itho-100*itth
                  iten=int(newid/10)-10*ihun-100*itho-1000*itth
                  iunit=newid-10*iten-100*ihun-1000*itho-10000*itth
                  if(sire(j) == '0') sire(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
                  if( dam(j) == '0')  dam(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
              end if
          end do
      end if

      !PRINT*,  ' Sorting Sires ... '
      allocate(SortedId(nobs), SortedIdIndex(nobs))
      SortedId(1:nobs) = Sire(1:nobs)
      Noffset = INT(nobs/2)
      DO WHILE (Noffset>0)
          Limit = nobs - Noffset
          switch=1
          DO WHILE (Switch.ne.0)
              Switch = 0
              do i = 1, Limit
                  IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
                      IDhold=SortedId(i)
                      SortedId(i)=SortedId(i + Noffset)
                      SortedId(i + Noffset)=IDhold
                      Switch = i
                  end if
              end do
              Limit = Switch - Noffset
          end do
          Noffset = INT(Noffset/2)
      end do

      ! Count number of unique sires
      nsires=0
      IF(SortedId(1) /= '0') nsires=1
      do i=2,nobs
          IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nsires=nsires+1
      end do

      ! Collect vector of unique sires in sorted order
      allocate  (SortedSire(0:nsires), SortedSireIndex(nsires))
      SortedSire(0) = '0'
      nsires=0

      IF(SortedId(1) /= '0') THEN
          nsires=1
          SortedSire(1) = SortedId(1)
      ENDIF

      do i=2,nobs
          IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
              nsires=nsires+1
              SortedSire(nsires) = SortedId(i)
          ENDIF
      end do

      !PRINT*,  ' Sorting Dams ... '
      SortedId(1:nobs) = Dam(1:nobs)
      Noffset = INT(nobs/2)
      DO WHILE (Noffset>0)
          Limit = nobs - Noffset
          switch=1
          DO WHILE (Switch.ne.0)
              Switch = 0
              do i = 1, Limit
                  IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
                      IDhold=SortedId(i)
                      SortedId(i)=SortedId(i + Noffset)
                      SortedId(i + Noffset)=IDhold
                      Switch = i
                  end if
              end do
              Limit = Switch - Noffset
          end do
          Noffset = INT(Noffset/2)
      end do

      nDams=0
      IF(SortedId(1) /= '0') nDams=1
      do i=2,nobs
          IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nDams=nDams+1
      end do

      allocate  (SortedDam(0:nDams), SortedDamIndex(ndams))
      SortedDam(0)='0'
      nDams=0
      IF(SortedId(1) /= '0') THEN
          nDams=1
          SortedDam(1) = SortedId(1)
      ENDIF

      do i=2,nobs
          IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
              nDams=nDams+1
              SortedDam(nDams) = SortedId(i)
          ENDIF
      end do

      !PRINT*,  ' Sorting IDs ... '
      SortedId(1:nobs) = ID(1:nobs)
      do i=1,nobs
          SortedIdIndex(i) = i
      end do
      Noffset = INT(nobs/2)
      DO WHILE (Noffset>0)
          Limit = nobs - Noffset
          switch=1
          DO WHILE (Switch.ne.0)
              Switch = 0
              do i = 1, Limit
                  IF (SortedId(i).gt.SortedId(i + Noffset)) THEN
                      IDhold=SortedId(i)
                      SortedId(i)=SortedId(i + Noffset)
                      SortedId(i + Noffset)=IDhold

                      ihold=SortedIdIndex(i)
                      SortedIdIndex(i)=SortedIdIndex(i + Noffset)
                      SortedIdIndex(i + Noffset)=ihold

                      Switch = i
                  end if
              end do
              Limit = Switch - Noffset
          end do
          Noffset = INT(Noffset/2)
      end do

      !PRINT*,  ' Check for duplicate IDs ... '
      flag = -1
      Do i = 2, nobs
          If (SortedID(i) == SortedID(i - 1)) Then
              If (flag == -1) Then
                  open (1,FILE='ID_err.txt',STATUS = 'unknown')
                  WRITE(1,*) 'Duplicated IDs ...'
                  flag = 0
              end if
              WRITE(1,*) SortedID(i)
              flag = flag + 1
          end if
      end do

      If (flag > -1) Then
          Close (1)
      ! PRINT*, flag,' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
      end if

      !PRINT*,  ' Males ... '
      !PRINT*,  '  Find or set sire indices ... '
      newsires = 0
      do j=1,nsires
          ! check if already listed as an individual
          ipoint=INT(nobs/2)
          Noffset = INT(ipoint/2)
          do while (Noffset>1)
              IF (SortedSire(j).lt.SortedId(ipoint)) THEN
                  ipoint = ipoint - Noffset
                  Noffset = INT(Noffset/2)
              else
                  ipoint = ipoint + Noffset
                  Noffset = INT(Noffset/2)
              end if
          end do

          kn=0
          if (SortedSire(j)==SortedId(ipoint)) kn=1
          do while (ipoint<nobs .and. kn==0 .and. SortedSire(j) > SortedId(ipoint))
              ipoint=ipoint+1
          end do
          if (SortedSire(j)==SortedId(ipoint)) kn=1
          do while (ipoint>1 .and. kn==0 .and. SortedSire(j) < SortedId(ipoint))
              ipoint=ipoint-1
          end do

          if (SortedSire(j)==SortedId(ipoint)) kn=1
          IF(kn==1) then
              SortedSireIndex(j) = SortedIdIndex(ipoint)
          else    ! sire is unlisted base sire
              newsires = newsires + 1
              SortedSireIndex(j) = nobs + newsires ! for now
          end if
      end do !j

      allocate  (holdsireid(newsires))
      kn=0
      do j=1,nsires
          if (SortedSireIndex(j) > nobs) then
              kn=kn+1
              holdsireid(SortedSireIndex(j)-nobs) = SortedSire(j)
          end if
      end do
      if (kn /= newsires) then
       print*,'newsires error'
       stop 1
      end if
      !PRINT*,  '  Find seqsire ... '
      do j = 1, nobs
          If (sire(j) == '0') Then
              seqsire(j)=0
          else
              ipoint=INT(nsires/2)
              Noffset = INT(ipoint/2)
              do while (Noffset>1)
                  IF (Sire(j).lt.SortedSire(ipoint)) THEN
                      ipoint = ipoint - Noffset
                      Noffset = INT(Noffset/2)
                  else
                      ipoint = ipoint + Noffset
                      Noffset = INT(Noffset/2)
                  end if
              end do
              kn=0
              if (Sire(j)==SortedSire(ipoint)) kn=1
              do while (ipoint<nsires .and. kn==0 .and. Sire(j) > SortedSire(ipoint))
                  ipoint=ipoint+1
              end do
              if (Sire(j)==SortedSire(ipoint)) kn=1
              do while (ipoint>1 .and. kn==0 .and. Sire(j) < SortedSire(ipoint))
                  ipoint=ipoint-1
              end do
              if (Sire(j)==SortedSire(ipoint)) kn=1
              IF(kn==1) then
                  seqsire(j) = SortedSireIndex(ipoint)
              else
                  PRINT*, ' Error: Sire missing: ', Sire(j)
                  stop 1
              end if
          end if
      ENDDO !j

      !PRINT*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'
      !PRINT*,  ' Females ... '
      !PRINT*,  '  Find or set dam indices ... '

      newdams = 0
      nbisexuals = 0

      do j=1,ndams
          ! check if already listed as an individual
          ipoint=INT(nobs/2)
          Noffset = INT(ipoint/2)
          do while (Noffset>1)
              IF (Sorteddam(j).lt.SortedId(ipoint)) THEN
                  ipoint = ipoint - Noffset
                  Noffset = INT(Noffset/2)
              else
                  ipoint = ipoint + Noffset
                  Noffset = INT(Noffset/2)
              end if
          end do
          kn=0
          if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint  ! store ipoint here as ipoint can change with bisexuals
          do while (ipoint<nobs .and. kn==0 .and. Sorteddam(j) > SortedId(ipoint))
              ipoint=ipoint+1
          end do

          if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
          do while (ipoint>1 .and. kn==0 .and. Sorteddam(j) < SortedId(ipoint))
              ipoint=ipoint-1
          end do

          if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
          ! check if already listed as a sire (and therefore bisexual)
          ipoint=INT(nsires/2)
          Noffset = INT(ipoint/2)
          do while (Noffset>1)
              IF (SortedDam(j).lt.SortedSire(ipoint)) THEN
                  ipoint = ipoint - Noffset
                  Noffset = INT(Noffset/2)
              else
                  ipoint = ipoint + Noffset
                  Noffset = INT(Noffset/2)
              end if
          end do

          kb=0
          if (SortedDam(j)==SortedSire(ipoint)) kb=1
          do while (ipoint<nsires .and. kb==0 .and. SortedDam(j) > SortedSire(ipoint))
              ipoint=ipoint+1
          end do

          if (SortedDam(j)==SortedSire(ipoint)) kb=1
          do while (ipoint>1 .and. kb==0 .and. SortedDam(j) < SortedSire(ipoint))
              ipoint=ipoint-1
          end do
          if (SortedDam(j)==SortedSire(ipoint)) kb=1

          IF(kb==1) then
              nbisexuals = nbisexuals + 1
              open (1,FILE='bisex.txt',position = 'append')
              WRITE(1,*) SortedDam(j)
              close(1)
          end if

          if (kb==1) then
              SorteddamIndex(j) = SortedSireIndex(ipoint)
          elseif (kn>=1) then
              SorteddamIndex(j) = SortedIdIndex(kn)
          else    ! dam is unlisted base dam
              newdams = newdams + 1
              SorteddamIndex(j) = nobs + newsires + newdams ! for now
          end if
      end do !j

      If (nbisexuals > 0)  then
        PRINT*, nbisexuals,' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'
      end if

      allocate(holddamid(newdams))

      kn=0
      do j=1,ndams
          if (SortedDamIndex(j) > nobs+newsires) then
              kn=kn+1
              holddamid(SortedDamIndex(j)-nobs-newsires) = SortedDam(j)
          end if
      end do

      if (kn /= newdams) then
        print*,'newdams error'
        stop 1
      end if

      !PRINT*,  '  Find seqdam ... '
      do j = 1, nobs
          If (dam(j) == '0') Then
              seqdam(j)=0
          else
              ipoint=INT(ndams/2)
              Noffset = INT(ipoint/2)
              do while (Noffset>1)
                  IF (dam(j).lt.Sorteddam(ipoint)) THEN
                      ipoint = ipoint - Noffset
                      Noffset = INT(Noffset/2)
                  else
                      ipoint = ipoint + Noffset
                      Noffset = INT(Noffset/2)
                  end if
              end do
              kn=0
              if (dam(j)==Sorteddam(ipoint)) kn=1
              do while (ipoint<ndams .and. kn==0 .and. dam(j) > Sorteddam(ipoint))
                  ipoint=ipoint+1
              end do
              if (dam(j)==Sorteddam(ipoint)) kn=1
              do while (ipoint>1 .and. kn==0 .and. dam(j) < Sorteddam(ipoint))
                  ipoint=ipoint-1
              end do
              if (dam(j)==Sorteddam(ipoint)) kn=1
              IF(kn==1) then
                  seqdam(j) = SorteddamIndex(ipoint)
              else
                  PRINT*, ' Error: dam missing: ', dam(j)
                  stop 1
              end if
          end if
      ENDDO !j

      !PRINT*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'
      !PRINT*,  ' Arranging unlisted base parents ... '

      iextra = newsires + newdams
      If (iextra > 0) then
              ! PRINT*, ' ', iextra, ' unlisted base parents found.'
          ! SortedId and SortedIdIndex just used as a holder while redimensioning
          SortedId(1:nobs)=id(1:nobs)
          deallocate (id)
          allocate(id(nobs+iextra))
          id(1+iextra:nobs+iextra)=SortedId(1:nobs)

          ! Do sires
          SortedId(1:nobs)=sire(1:nobs)
          deallocate (sire)
          allocate(sire(nobs+iextra))
          sire(1+iextra:nobs+iextra)=SortedId(1:nobs)
          SortedIdIndex(1:nobs)=seqsire(1:nobs)
          deallocate (seqsire)
          allocate(seqsire(nobs+iextra))
          seqsire(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)

          ! Do dams
          SortedId(1:nobs)=dam(1:nobs)
          deallocate (dam)
          allocate(dam(nobs+iextra))
          dam(1+iextra:nobs+iextra)=SortedId(1:nobs)
          SortedIdIndex(1:nobs)=seqdam(1:nobs)
          deallocate (seqdam)
          allocate(seqdam(nobs+iextra))
          seqdam(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)

          ! Do output switch, but use seqoutput as placeholder instead of SortedID
          !allocate(seqoutput(1:nobs))
          seqoutput(1:nobs) = dooutput(1:nobs)
          deallocate(dooutput)
          allocate(dooutput(nobs+iextra))
          dooutput(1+iextra:nobs+iextra)=seqoutput(1:nobs)
          if (nCols .eq. 4) dooutput(1:iextra) = 0
          if (nCols .eq. 3) dooutput(1:iextra) = 1
      end if

      !PRINT*, ' Inserting unlisted base parents ...'

      oldnobs = nobs
      nobs = nobs + iextra

      !PRINT*, ' Total number of animals = ',nobs

      allocate(passedorder(nobs))
      passedorder=0
      do i = 1+iextra, nobs
          passedorder(i)= i-iextra
          If (sire(i) == '0')then
              seqsire(i) = 0
          Else
              seqsire(i) = iextra + seqsire(i)
              If (seqsire(i) > nobs)  seqsire(i) = seqsire(i) - nobs  ! for unlisted sires
          end if

          If (dam(i) == '0') Then
              seqdam(i) = 0
          Else
              seqdam(i) = iextra + seqdam(i)
              If (seqdam(i) > nobs)  seqdam(i) = seqdam(i) - nobs
          end if
      ENDDO !i

      do i = 1, newsires
          ID(i) = holdsireid(i)
          passedorder(i)=0
          seqsire(i) = 0
          seqdam(i) = 0
      ENDDO !i

      do i = newsires + 1, newsires + newdams
          ID(i) = holddamid(i - newsires)
          passedorder(i)=0
          seqsire(i) = 0
          seqdam(i) = 0
      ENDDO !i

      deallocate(holdsireid, holddamid, SortedIdIndex, SortedId)

      flag = 0
      Do i = 1, nobs
          If (i <= seqsire(i) .Or. i <= seqdam(i) ) flag = 1
      end do !i
      !If (flag == 0) !PRINT*, 'not needed'!return

      !PRINT*, ' Re-Ordering pedigree ...'
      allocate( OldN(0:nobs), NewN(0:nobs) )
      allocate( holdid(0:nobs), holdsire(nobs), holddam(nobs), holdoutput(nobs) )

      OldN(0) = 0
      NewN=0
      !seqsire(0) = 0 !not needed !
      !seqdam(0) = 0

      holdid(1:nobs) = ID(1:nobs)
      holdsire = seqsire
      holddam = seqdam
      holdoutput(1:nobs) = dooutput(1:nobs)

      !Find base ancestors ...
      kn = 0
      do i = 1, nobs
          If (seqsire(i) == 0 .And. seqdam(i) == 0) Then
              kn = kn + 1
              NewN(i) = kn
              OldN(kn) = i
          end if
      ENDDO !i

      !Re-order pedigree ...
      NewN(0) = nobs + 1
      flag = 0
      Do While (kn < nobs)
          oldkn = kn
          do i = 1, nobs
              If (NewN(i) == 0) Then !And ID(i) <> 'UniqueNULL' Then
                  Ks = seqsire(i)
                  Kd = seqdam(i)
                  If (NewN(Ks) > 0 .And. NewN(Kd) > 0) Then
                      kn = kn + 1
                      NewN(i) = kn
                      OldN(kn) = i
                  end if
              end if
          end do !i
          ! to avoid hang on unexpected problem ...
          If (kn == oldkn) Then
              flag = flag + 1
          Else
              flag = 0
          end if

          If (flag > 10) Then
              open(1,file='ped_err.txt',status='unknown')
              write(1,*) 'Pedigree errors found involving two or more of the following relationships ...'
              write(1,*)
              write(1,*) '       Index numbers are followed by names.'
              write(1,*) '       Index number 0 means unknown, whence name is blank.'
              write(1,*)
              do i = 1, nobs
                  If (NewN(i) == 0) Then
                      write(1,*) 'Individual:',          i, ':  ', ID(i)
                      write(1,*) '    Father:', seqsire(i), ':  ', ID(seqsire(i))
                      write(1,*) '    Mother:',  seqdam(i), ':  ', ID(seqdam(i))
                      write(1,*)
                  end if
              ENDDO !i
              Close (1)
              PRINT*,  'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
              stop 1
          end if
      ENDDO !while

      NewN(0) = 0
      do i = 1, nobs
          ID(i) = holdid(OldN(i))
          dooutput(i) = holdoutput(OldN(i))
      end do

      do i = 1, nobs
          seqsire(i) = NewN(holdsire(OldN(i)))
          seqdam(i) = NewN(holddam(OldN(i)))
          If (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
              PRINT*,  'out of order'
              stop 1
          end if
      ENDDO !i

      DO i = 1, nobs
          holdsire(i) = passedorder(i)  ! holdsire just because it is free
      end do

      DO i = 1, nobs
          passedorder(i) = holdsire(OldN(i))
      end do

      deallocate ( OldN, NewN, holdid, holdsire, holddam) ! holdrec)
      !do i = 1, nobs
      ! PRINT'(3i5,2x,3a4,i5)', i, seqsire(i), seqdam(i), id(i), sire(i), dam(i), passedorder(i)
      !end do

      nAnisPedigree=nObs
      GlobalExtraAnimals=iextra   !Change John Hickey
    end subroutine

  !#############################################################################

    subroutine invert(x,n,sym, method)

      ! Interface to call inverse subroutines from BLAS/LAPACK libraries

      ! x symmetric positive-definite matrix to be inverted
      ! n matrix dimension
      ! sym return lower-triangular (sym=.false) or full matrix (sym=.true.)
      ! method for inversion
      ! 0 -- Generalised solving using LU decomposition (dgetrs)
      ! 1 -- Cholesky decomposition

      implicit none
      integer(int32), intent(in) :: n,method
      integer(int32) :: i,j,info

      real(real64),intent(inout) :: x(n,n)
      real(real64),allocatable :: Iden(:,:)

      logical, intent(in) :: sym

      if (method == 0) then
        !Solves a general system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by SGETRF/CGETRF
        !http://physics.oregonstate.edu/~rubin/nacphy/lapack/routines/dgetrs.html

        allocate(Iden(n,n))
        ForAll(i = 1:n, j = 1:n) Iden(i,j) = (i/j)*(j/i)  !https://rosettacode.org/wiki/Identity_matrix#Notorious_trick

        !https://software.intel.com/en-us/node/468712
        !Solves a system of linear equations with an LU-factored square coefficient matrix, with multiple right-hand sides.
        ! dgetrs(trans,n,nrhs,A,b,lda,ldb,info)
        !Output: Solution overwrites `b`.
        call dgetrs('N',n,n,x,Iden,n,n,info)
        if (info /= 0) then
          print *, 'Matrix not positive-definite - info',info
          stop 1
        end if

        x(:,:) = Iden(:,:)

      else if (method == 1) then

        ! Computes the Cholesky factorization of a symmetric positive definite matrix
        ! https://software.intel.com/en-us/node/468690
        call dpotrf('L',n,x,n,info)
        if (info /= 0) then
          print*,'Matrix not positive-definite - info',info
          stop 1
        end if

        ! Computes the inverse of a symmetric positive definite matrix,
        !   using the Cholesky factorization computed by dpotrf()
        ! https://software.intel.com/en-us/node/468824
        call dpotri('L',n,x,n,info)
        if (info /= 0) then
         print*,'Matrix not positive-definite - info',info
         stop 1
        end if

        ! Fills the upper triangle
        if (sym) then
          forall (i=1:n,j=1:n,j>i) x(i,j)=x(j,i)
        end if

      end if
    end subroutine

  !#############################################################################

    subroutine AlphaAGHTitles
       print*, ""
       write(*,'(a30,a,a30)') " ","**********************"," "
       write(*,'(a30,a,a30)') " ","*                    *"," "
       write(*,'(a30,a,a30)') " ","*      AlphaAGH      *"," "
       write(*,'(a30,a,a30)') " ","*                    *"," "
       write(*,'(a30,a,a30)') " ","**********************"
       write(*,'(a30,a,a30)') " ","VERSION:"//TOSTRING(VERS)," "
       print*, "              Software for bulding numerator relationship matrices       "
       print*, ""
       print*, "                                    No Liability              "
       print*, "                          Bugs to John.Hickey@roslin.ed.ac.uk"
       print*, ""
    end subroutine

  !#############################################################################
end module

!###############################################################################

program AlphaAGH

  use AlphaAHGModule

  implicit none

  include "mkl_vml.f90"

  real(real32) :: start, finish

  call cpu_time(start)
  call AlphaAGHTitles

  call ReadParam
  call ReadData

  if (AMake) then
    call MakeAMatrix
  end if

  if (AInvMake) then
    call MakeInvAMatrix
  end if

  if (GMake .or. GInvMake) then
    call MakeG
  end if

  if (HInvMake .or. HMake) then
    call MakeHAndHInv
  end if

  call cpu_time(finish)
  print *," "
  print '("  Time duration of AlphaAGH = ",f20.4," seconds.")',finish-start
  print *," "

end program

!###############################################################################
