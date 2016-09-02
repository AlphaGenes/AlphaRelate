
!TODO: write output subroutines and call them in different places

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

!#########################################################################

module AlphaRelateMod

  use ISO_Fortran_env
  use AlphaHouseMod, only : CountLines

  INTEGER(int32),PARAMETER :: LENGAN=20

  integer(int32) :: nTrait,nSnp,nAnisG,nAnisP,nAnisRawPedigree,AllFreqSelCycle,nCols,nAnisH
  integer(int32) :: nGMats, GlobalExtraAnimals, OldAMatNInd
  integer(int32) :: NRMmem, shell, shellmax, shellWarning
  integer(int32) :: PedigreeUnit,OldAMatUnit,GenotypeUnit,AlleleFreqUnit,WeightUnit
  integer(int32),allocatable :: MapAnimal(:),seqid(:),seqsire(:),seqdam(:),seqoutput(:)
  integer(int32),allocatable :: RecodeGenotypeId(:),passedorder(:),RecPed(:,:),dooutput(:)
  integer(int32),allocatable :: OldAMatId(:)

  real(real64) :: AlleleFreqAll,DiagFudge,ScaleGToA
  real(real64),allocatable :: AlleleFreq(:),Adiag(:)
  real(real64),allocatable :: Weights(:,:),Zmat(:,:),tpose(:,:),Genos(:,:),AMat(:,:),InvAMat(:,:)
  real(real64),allocatable :: WeightStand(:,:,:),GMat(:,:,:),InvGMat(:,:,:)

  character(len=20) :: GType,OutputFormat
  character(len=1000) :: PedigreeFile,WeightFile,GenotypeFile,AlleleFreqFile,OldAMatFile
  character(len=LENGAN),allocatable :: ped(:,:),Id(:),sire(:),dam(:),IdGeno(:)

  logical :: PedigreePresent,GenotypesPresent,AlleleFreqPresent,AlleleFreqFixed,WeightsPresent,OldAMatPresent,ScaleGByRegression
  logical :: MakeG, MakeInvG, MakeA, MakeInvA, MakeH, MakeInvH
  logical :: GFullMat, GIJA, InvGFullMat, InvGIJA, AFullMat, AIJA, InvAFullMat, InvAIJA, HFullMat, HIJA, InvHFullMat, InvHIJA
  logical,allocatable :: MapToG(:), AnimalsInBoth(:)

  contains

  !#############################################################################

    subroutine ReadParam
      implicit none

      integer(int32) :: i,j,n,OutputPositions,OutputDigits,SpecUnit

      character(len=1000) :: dumC,option
      character(len=200) :: OutputPositionsC,OutputDigitsC

      open(newunit=SpecUnit,file="AlphaRelateSpec.txt",status="old")

      ! Input parameters
      read(SpecUnit,*) dumC

      read(SpecUnit,*) dumC,GenotypeFile
      GenotypesPresent=.false.
      if (trim(GenotypeFile)/='None') then
        GenotypesPresent=.true.
      end if

      read(SpecUnit,*) dumC,PedigreeFile
      PedigreePresent=.false.
      if (trim(PedigreeFile)/='None') then
        PedigreePresent=.true.
      end if

      read(SpecUnit,*) dumC,WeightFile
      WeightsPresent = .false.
      if (trim(WeightFile)/='None') then
        WeightsPresent = .true.
      end if

      read(SpecUnit,*) dumC,AlleleFreqFile
      AlleleFreqPresent = .false.
      if (trim(AlleleFreqFile)/='None') then
        AlleleFreqPresent = .true.
      end if
      AlleleFreqFixed = .false.
      if (trim(AlleleFreqFile) == "Fixed") then
        AlleleFreqFixed = .true.
        backspace(SpecUnit)
        read(SpecUnit,*) dumC,dumC,AlleleFreqAll
      end if

      read(SpecUnit,*) dumC,nTrait
      read(SpecUnit,*) dumC,nSnp
      read(SpecUnit,*) dumC,DiagFudge

      read(SpecUnit,*) dumC,GType
      if (GenotypesPresent                 .and.&
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

      read(SpecUnit,*) dumC,option
      if (trim(option) == 'Regression') then
        ScaleGByRegression = .true.
      else
        ScaleGByRegression = .false.
        read(option,*) ScaleGtoA
      end if

      ! Output options
      read(SpecUnit,*) dumC

      read(SpecUnit,*) dumC,OutputPositions,OutputDigits
      write(OutputPositionsC,*) OutputPositions
      write(OutputDigitsC,*)    OutputDigits
      OutputFormat="f"//trim(adjustl(OutputPositionsC))//"."//trim(adjustl(OutputDigitsC))

      ! Make G matrix?
      GFullMat = .false.
      GIJA = .false.
      MakeG = .false.
      read(SpecUnit,*) dumC, option
      GFullMat = trim(option) == 'Yes'
      if (GFullMat) then
        MakeG = .true.
      end if

      read(SpecUnit,*) dumC, option
      GIJA = trim(option) == 'Yes'
      if (GIJA) then
        MakeG = .true.
      end if

      ! Make inverted G matrix ?
      InvGFullMat = .false.
      InvGIJA = .false.
      MakeInvG = .false.
      read(SpecUnit,*) dumC, option
      InvGFullMat = trim(option) == 'Yes'
      if (InvGFullMat) then
        MakeInvG = .true.
      end if

      read(SpecUnit,*) dumC, option
      InvGIJA = trim(option) == 'Yes'
      if (InvGIJA) then
        MakeInvG = .true.
      end if

      ! Make A matrix?
      AFullMat = .false.
      AIJA = .false.
      MakeA = .false.
      read(SpecUnit,*) dumC, option
      AFullMat = trim(option) == 'Yes'
      if (AFullMat) then
        MakeA = .true.
      end if

      read(SpecUnit,*) dumC, option
      AIJA = trim(option) == 'Yes'
      if (AIJA) then
        MakeA = .true.
      end if

      ! Make inverted A matrix?
      InvAFullMat = .false.
      InvAIJA = .false.
      MakeInvA = .false.
      read(SpecUnit,*) dumC, option
      InvAFullMat = trim(option) == 'Yes'
      if (InvAFullMat) then
        MakeInvA = .true.
      end if

      read(SpecUnit,*) dumC, option
      InvAIJA = trim(option) == 'Yes'
      if (InvAIJA) then
        MakeInvA = .true.
      end if

      ! Make H matrix?
      HFullMat = .false.
      HIJA = .false.
      MakeH = .false.
      read(SpecUnit,*) dumC, option
      HFullMat = trim(option) == 'Yes'
      if (HFullMat) then
        MakeH = .true.
      end if

      read(SpecUnit,*) dumC, option
      HIJA = trim(option) == 'Yes'
      if (HIJA) then
        MakeH = .true.
      end if

      ! Make inverted H matrix?
      InvHFullMat = .false.
      InvHIJA = .false.
      MakeInvH = .false.
      read(SpecUnit,*) dumC, option
      InvHFullMat = trim(option) == 'Yes'
      if (InvHFullMat) then
        MakeInvH = .true.
      end if

      read(SpecUnit,*) dumC, option
      InvHIJA = trim(option) == 'Yes'
      if (InvHIJA) then
        MakeInvH = .true.
      end if

      if (MakeH) then
        MakeG = .true.
        MakeA = .true.
      end if

      if (MakeInvH) then
        MakeG = .true.
        MakeA = .true.
        MakeInvA = .true.
      end if

      OldAMatPresent = .false.
      n=CountLines("AlphaRelateSpec.txt")
      if (n > 24) then
        print *, "This is experimental feature/hack = not well tested and might be removed !!!"
        print *, "  - It requires id of individuals to be numeric and sequential and no unknown parents"
        print *, "  - It requires the old A matrix between the parents of individuals whose A matrix will be built"
        print *, "  - It switches off creation of other matrices (exit after AMat is done)"
        read(SpecUnit,*) dumC, OldAMatFile, OldAMatNInd
        OldAMatPresent = .true.
        open(newunit=OldAMatUnit, file=OldAMatFile, status="unknown")
      end if

      nGMats=0
      do i=1,nTrait
        do j=i,nTrait
          nGMats=nGMats+1
        end do
      end do

      if (.not. GenotypesPresent) then
        if (MakeG .or. MakeInvG .or. MakeH .or. MakeInvH) then
          print *, 'In order to create G or H matrices, a genotype file must be given.'
        end if
        MakeG=.false.
        MakeInvG=.false.
        MakeH=.false.
        MakeInvH=.false.
      end if

      if (.not. PedigreePresent) then
        if (MakeA .or. MakeInvA .or. MakeH .or. MakeInvH) then
          print *, 'In order to create A or H matrices, a pedigree must be given.'
        end if
        MakeA=.false.
        MakeInvA=.false.
        MakeH=.false.
        MakeInvH=.false.
      end if

      if (ScaleGByRegression) then
        MakeA = .true.
      end if
    end subroutine

  !#############################################################################

    subroutine ReadData
      implicit none

      integer(int32) :: i,j,stat,GenoInPed,fourthColumn,nMissing

      character(len=1000) :: dumC
      character(len=LENGAN), dimension(1:3) :: pedline

      if (PedigreePresent) then
        nAnisRawPedigree=CountLines(trim(PedigreeFile))
        write(*,'(a2,i6,a33)') "   ",nAnisRawPedigree," individuals in the pedigree file"

        !!! Attempt to magically detect whether there are three or four columns:
        open(newunit=PedigreeUnit,file=trim(PedigreeFile),status="old")
        read(PedigreeUnit, '(a)', iostat=stat) dumC
        if (stat /= 0) then
          print *, 'Problems reading Pedigree file.'
          !print *, stat, dumC
          stop 1
        end if
        rewind(PedigreeUnit)

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
          read(PedigreeUnit,*) ped(i,1:nCols)
        end do
        close(PedigreeUnit)
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

      if (GenotypesPresent) then
        nAnisG=CountLines(trim(GenotypeFile))
        write(*,'(a2,i6,a33)') "   ",nAnisG," individuals in the genotype file"
        allocate(Genos(nAnisG,nSnp))
        allocate(Zmat(nAnisG,nSnp))
        allocate(IdGeno(nAnisG))
        !allocate(RecodeIdGeno(nAnisG))
        open(newunit=GenotypeUnit,file=trim(GenotypeFile),status="old")
        do i=1,nAnisG
          read(GenotypeUnit,*) IdGeno(i),Genos(i,:)
        end do
        close(GenotypeUnit)

        ! Allele frequencies
        allocate(AlleleFreq(nSnp))
        if (.not. AlleleFreqPresent) then
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
          open(newunit=AlleleFreqUnit, file='AlleleFreqTest.txt',status='unknown')
          do j=1,nSnp
            write(AlleleFreqUnit,*) j,AlleleFreq(j)
          end do
          close(AlleleFreqUnit)
        else
          if (trim(AlleleFreqFile) == "Fixed") then
            AlleleFreq(:) = AlleleFreqAll
            open(newunit=AlleleFreqUnit, file='AlleleFreqTest.txt',status='unknown')
            do j=1,nSnp
              write(AlleleFreqUnit,*) j,AlleleFreq(j)
            end do
            close(AlleleFreqUnit)
          else
            ! Read allele frequencies from file.
            open(newunit=AlleleFreqUnit, file=trim(AlleleFreqFile), status='OLD')
            do i=1,nSnp
              read(AlleleFreqUnit,*,iostat=stat) dumC,AlleleFreq(i)  !AlleleFrequencies are kept in second column to keep consistency with AlphaSim.
              if (stat /= 0) then
                print*,"Problems reading allele frequency file."
                stop 1
              end if
            end do
            close(AlleleFreqUnit)
          end if
        end if

        ! Weights
        allocate(Weights(nSnp,nTrait))
        allocate(WeightStand(nSnp,nTrait,nTrait))
        if (WeightsPresent) then
          open(newunit=WeightUnit,file=trim(WeightFile),status="old")
          do i=1,nSnp
            read(WeightUnit,*) dumC,Weights(i,:)
          end do
          close(WeightUnit)
        else
          Weights(:,:)=1
        end if
      end if

      if (.not. PedigreePresent .and. GenotypesPresent) then
        allocate(RecPed(0:nAnisG,4))
        nAnisP=nAnisG
        RecPed(:,:)=0
        RecPed(:,4)=1
        do i=1,nAnisP
          RecPed(i,1)=i
        end do
      end if

      if (PedigreePresent .and. GenotypesPresent) then
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
    end subroutine

  !#############################################################################

    subroutine MakeInvAMatrix
      implicit none

      integer(int32) :: i,m,n,s,FId,MId

      real(real64) :: Inbreeding(0:nAnisP),Dii(nAnisP),InvDii

      character(len=1000) :: nChar,fmt

      logical :: AnimToWrite(nAnisP),FIdL,MIdL

      allocate(InvAMat(0:nAnisP,0:nAnisP))

      print*, "Start calculating inbreeding coefficients"
      call dinbreeding(RecPed(0:nAnisP,1),RecPed(0:nAnisP,2),RecPed(0:nAnisP,3),Inbreeding,nAnisP)
      open(unit=202,file="PedigreeBasedInbreeding.txt",status="unknown")
      print*, "End calculating inbreeding coefficients"
      do i=1,nAnisP
        write(202,'(a20,20000f10.5)') Id(i),Inbreeding(i)
      end do
      close(202)

      print*, "Start making A inverse"
      InvAMat=0.0d0
      ! TODO: could remove the if statements bellow since InvAMat has zeroth row and column
      !       and could simply increment values in those positions - they are omitted on
      !       printout anyhow
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
        InvAMat(i,i)=InvDii
        ! Add precision to the father and set the co-precision
        if (FIdL) then
          InvAMat(FId,FId)=InvAMat(FId,FId)+InvDii/4.0d0
          InvAMat(i,FId)=InvAMat(i,FId)-InvDii/2.0d0
          InvAMat(FId,i)=InvAMat(i,FId)
        end if
        ! Add precision to the mother and set the co-precision
        if (MIdL) then
          InvAMat(MId,MId)=InvAMat(MId,MId)+InvDii/4.0d0
          InvAMat(i,MId)=InvAMat(i,MId)-InvDii/2.0d0
          InvAMat(MId,i)=InvAMat(i,MId)
        end if
        ! Add co-precision between the father and mother
        if (FIdL .and. MIdL) then
          InvAMat(FId,MId)=InvAMat(FId,MId)+InvDii/4.0d0
          InvAMat(MId,FId)=InvAMat(FId,MId)
        end if
      end do
      print*, "Finished making A inverse"

      if (InvAFullMat) then
        AnimToWrite = RecPed(1:nAnisP,4) == 1
        s = count(AnimToWrite)
        write(*,'(a40,i6,a11)') " Start writing A inverse full matrix for", s," individuals"
        write(nChar,*) s
        fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="InvAFullMatrix.txt",status="unknown")
        do m=1,nAnisP
          if (AnimToWrite(m)) then
            write(202,fmt) Id(m), pack(InvAMat(1:nAnisP,m), AnimToWrite)
          end if
        end do
        close(202)
        print*, "End writing A inverse full matrix"
      end if

      if (InvAIJA) then
        print*, "Start writing A inverse ija"
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="InvAija.txt",status="unknown")
        do m=1,nAnisP
          do n=m,nAnisP
            if (InvAMat(n,m) /= 0.0d0) then
              write(202,fmt) Id(n),Id(m),InvAMat(n,m)
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

      real(real64) :: AMatAvg

      character(len=1000) :: nChar,fmt

      logical :: AnimToWrite(nAnisP)

      print*, "Start making A"
      if (OldAMatPresent) then
        allocate(OldAMatId(OldAMatNInd))
        do j = 1, OldAMatNInd
          read(OldAMatUnit, *) OldAMatId(j)
        end do
        rewind(OldAMatUnit)
        MinId = minval(OldAMatId)
        MaxId = maxval(OldAMatId)
        allocate(AMat(1:(OldAMatNInd+nAnisP-MaxId),&
                      1:(OldAMatNInd+nAnisP-MaxId)))
        AMat = 0.0d0
        do j = 1, OldAMatNInd
          read(OldAMatUnit, *) OldAMatId(j), AMat(1:OldAMatNInd,j)
          if (j > 1) then
            if (.not.(OldAMatId(j) > OldAMatId(j-1))) then
              print *, "Id are not sequential!"
              stop 1
            end if
          end if
        end do
        k = OldAMatNInd
        do i=MaxId+1,nAnisP
            k = k + 1
            s = RecPed(i,2) - MinId + 1
            d = RecPed(i,3) - MinId + 1
            l = k
            do j=1,k-1
                AMat(j,k)=(AMat(j,s)+AMat(j,d))/2.0d0
                AMat(k,j)=AMat(j,k)
                !print *,i,k,j,s,d,AMat(j,s),AMat(j,d),AMat(j,k)
            end do
            AMat(k,k)=1.0d0+AMat(s,d)/2.0d0
            !print *,i,k,s,d,AMat(s,d),AMat(k,k)
        end do
        RecPed(1:nAnisP,4) = 0
        RecPed((MaxId+1):nAnisP,4) = 1
      else
        allocate(AMat(0:nAnisP,0:nAnisP))
        AMat=0.0d0
        do i=1,nAnisP
            do j=1,i-1
                AMat(j,i)=(AMat(j,RecPed(i,2))+AMat(j,RecPed(i,3)))/2.0d0
                AMat(i,j)=AMat(j,i)
            end do
            AMat(i,i)=1.0d0+AMat(RecPed(i,2),RecPed(i,3))/2.0d0
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
        if (.not.OldAMatPresent) then
          do m=1,nAnisP
            if (AnimToWrite(m)) then
              write(202,fmt) Id(m), pack(AMat(1:nAnisP,m), AnimToWrite)
            end if
          end do
        else
          Start = OldAMatNInd+1
          Endin = size(AMat,1)
          do m=Start,Endin
            !write(*,fmt)   Id(m+MinId-1), AMat(Start:Endin,m)
            write(202,fmt) Id(m+MinId-1), AMat(Start:Endin,m)
          end do
        end if
        close(202)
        print*, "End writing A full matrix"
      end if

      if (OldAMatPresent) then
        stop
      end if

      if (AIJA) then
        ! TODO: AnimToWrite is not being used here
        write(*,'(a24,i6,a11)') " Start writing A ija for", s," individuals"
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="Aija.txt",status="unknown")
        do m=1,nAnisP
          do n=m,nAnisP
            if (AMat(n,m) > 0.0d0)  then
              write(202,fmt) Id(n),Id(m),AMat(n,m)
            end if
          end do
        end do
        close(202)
        print*, "End writing A ija"
      end if

      ! Record diagonals of animals in both A and G:
      if ((MakeH .or. MakeInvH) .and. ScaleGByRegression) then
        n = Count(AnimalsInBoth)
        allocate(Adiag(0:n))
        div = dble(n**2)
        AMatAvg = 0.0d0
        k = 0
        do i = 1,nAnisP
          if (.not. AnimalsInBoth(i)) then
            cycle
          end if
          k = k + 1
          Adiag(k) = AMat(i,i)
          do j=1,nAnisP
            if (AnimalsInBoth(j)) then
              AMatAvg=AMatAvg + AMat(j,i) * 2.0d0 / div
            end if
          end do
        end do
        Adiag(0) = AMatAvg
      end if
    end subroutine

  !#############################################################################

    subroutine MakeGAndInvGMatrix
      implicit none

      integer(int32) :: i,j,k,l,m,n,WhichMat

      real(real64) :: Tmp,TmpWt(nSnp),TmpVal(1),Denom

      character(len=1000) :: filout,nChar,fmt

      allocate(GMat(nAnisG,nAnisG,nGMats))
      allocate(tpose(nSnp,nAnisG))

      print*, "Start making G - ", trim(GType)

      !Standardise weights
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

          if (MakeInvG) then
            allocate(InvGMat(nAnisG,nAnisG,nGMats))

            print*, "Start inverting G - ", trim(GType)
            InvGMat(:,:,WhichMat)=GMat(:,:,WhichMat)
            call invert(InvGMat(:,:,WhichMat),nAnisG,.true., 1)

            print*, "Finished inverting G - ", trim(GType)

            if (InvGFullMat) then
              write(nChar,*) nAnisG
              fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
              write(filout,'("InvGFullMatrix"i0,"-"i0".txt")')i,j
              open(unit=202,file=trim(filout),status="unknown")
              do m=1,nAnisG
                write(202,fmt) IdGeno(m),InvGMat(:,m,WhichMat)
              end do
              close(202)
            end if

            if (InvGIJA) then
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

    subroutine MakeHAndInvHMatrix
      ! Feature added by Stefan Hoj-Edwards, February 2016
      ! Making the Inverse H matrix ala Aguilar et al 201? and Christensen 2012

      ! Prerequisite and assumptions for this subroutine:
      ! There is given both a pedigree and genotype file, and there is an overlap
      ! of animals between two data sets.
      ! Diagonals of from both A have been collected during MakeA and MakeG,
      ! as well as average of A22.
      ! GMat is already calculated and loaded in memory.
      ! InvA is calculated an loaded into memory.
      !
      ! Further assumes that animals are ordered the same in both A and G.

      implicit none

      integer(int32) :: i,j,k,m,p,q,div,t1,t2,whichMat,nBoth
      integer(int32),allocatable :: MapToA11(:), MapToA22(:) !Gmap(:),

      real(real64) :: GMatavg, nom, denom, slope, intercept, Gmean, Amean, Hii
      real(real64),allocatable :: Gdiag(:), Hrow(:), A22(:,:), InvA22(:,:), G22(:,:), A11(:,:), A12(:,:), tmp(:,:), Gboth(:,:)

      character(len=1000) :: nChar,fmt1, fmt2,filout
      character(len=LENGAN),allocatable :: Ids(:)

      logical,allocatable :: AnimToWrite(:)

      nboth = count(AnimalsInBoth)
      ! Make H and/or InvH
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

      allocate(InvA22(nBoth,nBoth))
      allocate(MapToA22(nAnisH))
      if (MakeH) then
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
          InvA22(m,k) = AMat(j,i)
        end do
      end do
      if (MakeH) then
        A22 = InvA22
      end if

      call invert(InvA22,size(InvA22,1),.true.,1)

      ! This is the G matrix in Legarra,
      ! Sadly, no genotypes where provided, instead the resulting G matrix was.
      if (.false.) then
        print *, 'Overwriting G matrix with example in Legarra 2008!'
        do i=1,nAnisG
          do j=1,nAnisG
            if (i==j) then
              GMat(i,j,1) = 1
            else
              GMat(i,j,1) = 0.7
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
              nom = GMat(i,j,whichMat)
              if (i == j) then
                nom = nom - DiagFudge
              end if
              G22(i,j) = nom
            end do
          end do

          if (ScaleGByRegression) then
            allocate(Gdiag(0:nBoth))
            Gdiag=0.0d0
            GMatavg=0.0d0
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
                GMatavg=GMatavg + G22(MapAnimal(j),MapAnimal(i))/div
              end do
            end do
            Gdiag(0) = GMatavg

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
                  G22(MapAnimal(j),MapAnimal(i)) = ScaleGToA * G22(MapAnimal(j),MapAnimal(i)) + (1.0d0 - ScaleGToA) * AMat(j,i)
                end if
              end do
            end do
          end if

          do i=1,nAnisG
            G22(i,i) = G22(i,i) + DiagFudge
          end do

          allocate(Hrow(1:count(AnimToWrite)))

          if (MakeH) then

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
                    A12(k,MapAnimal(j)) = AMat(j,i)  !A12 is not symmetrical
                  else
                    m = m+1
                    A11(m,k) = AMat(j,i)
                  end if
                end do
              end if
            end do

            tmp = Matmul(A12, InvA22)
            !tmp = Matmul(Matmul(tmp, (Gboth - A22)), transpose(tmp))
            tmp = Matmul(tmp, (Gboth-A22))
            tmp = Matmul(tmp, InvA22)
            !tmp = Matmul(tmp, transpose(A12))

            A11 = A11 + Matmul(tmp, transpose(A12))
            A12 = matmul(matmul(A12, InvA22), Gboth)

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
                if (InvHIJA .and. i .le. j .and. Hii /= 0.0d0) then
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

          if (MakeInvH) then
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
            ! write(*,fmt1) IdGeno(i), InvA22(i,:)
            !end do

            !print *, 'InvA(22)'
            !do i=1,size(G22,1)
            ! j = i+10
            ! write(*,fmt1) Ids(j), InvAMat(j,11:25)
            !end do

            print *, 'End inverting scaled G matrix'

            print *, 'Start writing inverted H matrices (full and/or ija)'

            if (InvHFullMat) then
              write(filout,'("InvHFullMatrix"i0,"-"i0".txt")') t1,t2
              write(nChar,*) nAnisH
              fmt1="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
              open(unit=202,file=trim(filout),status="unknown")
            end if

            if (InvHIJA) then
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
                    Hrow(k) = Hrow(k) + InvAMat(i,j) - InvA22(MapToA22(i),MapToA22(j))
                  end if
                else if (i <= nAnisP .and. j <= nAnisP) then !if (MapToG(i) .eq. .false. .and. MapToG(j) .eq. .false.  ) then
                  Hrow(k) = InvAMat(i,j)
                end if
                if (InvHIJA .and. i .le. j .and. Hrow(k) /= 0.0d0) then
                  write(204,fmt2) trim(Ids(i)), trim(Ids(j)), Hrow(k)
                end if
              end do
              if (InvHFullMat) then
                write(202,fmt1) Ids(i),Hrow(:)
              end if
            end do

            if (InvHFullMat) then
              close(202)
            end if
            if (InvHIJA) then
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
      integer(int32) :: IdErrUnit, BisexUnit, PedErrUnit
      integer(int32),allocatable :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
      integer(int32),allocatable :: OldN(:), NewN(:), holdsire(:), holddam(:), holdoutput(:)

      character(len=LENGAN)              :: IDhold
      !character(len=LENGAN) :: path
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
      !path=".\"

      Verbose=1

      do j = 1, nobs
          if (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') Then
              dam(j) = '0'
              seqdam(j)=0
          end if
          if (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') Then
              sire(j) = '0'
              seqsire(j)=0
          end if
      end do !j

      if(mode.eq.1) then
          !print*,  ' Inserting dummy IDs ... '
          newid=0
          do j = 1, nobs
              if(((sire(j) == '0').and.(dam(j).ne.'0'))  .or. ((sire(j).ne.'0').and.(dam(j) == '0'))) then
                  newid=newid+1
                  if(newid.gt.99999) then
                      !         print*, newid, ' ...'
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

      !print*,  ' Sorting Sires ... '
      allocate(SortedId(nobs), SortedIdIndex(nobs))
      SortedId(1:nobs) = Sire(1:nobs)
      Noffset = INT(nobs/2)
      do while (Noffset>0)
          Limit = nobs - Noffset
          switch=1
          do while (Switch.ne.0)
              Switch = 0
              do i = 1, Limit
                  IF (SortedId(i).gt.SortedId(i + Noffset)) then
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
      if (SortedId(1) /= '0') nsires=1
      do i=2,nobs
          if (SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nsires=nsires+1
      end do

      ! Collect vector of unique sires in sorted order
      allocate  (SortedSire(0:nsires), SortedSireIndex(nsires))
      SortedSire(0) = '0'
      nsires=0

      if (SortedId(1) /= '0') then
          nsires=1
          SortedSire(1) = SortedId(1)
      end if

      do i=2,nobs
          if (SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
              nsires=nsires+1
              SortedSire(nsires) = SortedId(i)
          end if
      end do

      !print*,  ' Sorting Dams ... '
      SortedId(1:nobs) = Dam(1:nobs)
      Noffset = INT(nobs/2)
      do while (Noffset>0)
          Limit = nobs - Noffset
          switch=1
          do while (Switch.ne.0)
              Switch = 0
              do i = 1, Limit
                  IF (SortedId(i).gt.SortedId(i + Noffset)) then
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
      if (SortedId(1) /= '0') nDams=1
      do i=2,nobs
          if (SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nDams=nDams+1
      end do

      allocate  (SortedDam(0:nDams), SortedDamIndex(ndams))
      SortedDam(0)='0'
      nDams=0
      if (SortedId(1) /= '0') then
          nDams=1
          SortedDam(1) = SortedId(1)
      end if

      do i=2,nobs
          if (SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') then
              nDams=nDams+1
              SortedDam(nDams) = SortedId(i)
          end if
      end do

      !print*,  ' Sorting IDs ... '
      SortedId(1:nobs) = ID(1:nobs)
      do i=1,nobs
          SortedIdIndex(i) = i
      end do
      Noffset = INT(nobs/2)
      do while (Noffset>0)
          Limit = nobs - Noffset
          switch=1
          do while (Switch.ne.0)
              Switch = 0
              do i = 1, Limit
                  IF (SortedId(i).gt.SortedId(i + Noffset)) then
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

      !print*,  ' Check for duplicate IDs ... '
      flag = -1
      Do i = 2, nobs
          if (SortedID(i) == SortedID(i - 1)) Then
              if (flag == -1) Then
                  open(IdErrUnit,FILE='ID_err.txt',STATUS = 'unknown')
                  write(IdErrUnit,*) 'Duplicated IDs ...'
                  flag = 0
              end if
              write(IdErrUnit,*) SortedID(i)
              flag = flag + 1
          end if
      end do

      if (flag > -1) Then
          close(IdErrUnit)
      ! print*, flag,' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
      end if

      !print*,  ' Males ... '
      !print*,  '  Find or set sire indices ... '
      newsires = 0
      do j=1,nsires
          ! check if already listed as an individual
          ipoint=INT(nobs/2)
          Noffset = INT(ipoint/2)
          do while (Noffset>1)
              IF (SortedSire(j).lt.SortedId(ipoint)) then
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
          if (kn==1) then
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
      !print*,  '  Find seqsire ... '
      do j = 1, nobs
          if (sire(j) == '0') Then
              seqsire(j)=0
          else
              ipoint=INT(nsires/2)
              Noffset = INT(ipoint/2)
              do while (Noffset>1)
                  IF (Sire(j).lt.SortedSire(ipoint)) then
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
              if (kn==1) then
                  seqsire(j) = SortedSireIndex(ipoint)
              else
                  print*, ' Error: Sire missing: ', Sire(j)
                  stop 1
              end if
          end if
      end do !j

      !print*,  '  Sires: ',newsires,' unlisted, ',nsires,' in total'
      !print*,  ' Females ... '
      !print*,  '  Find or set dam indices ... '

      newdams = 0
      nbisexuals = 0

      do j=1,ndams
          ! check if already listed as an individual
          ipoint=INT(nobs/2)
          Noffset = INT(ipoint/2)
          do while (Noffset>1)
              IF (Sorteddam(j).lt.SortedId(ipoint)) then
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
              IF (SortedDam(j).lt.SortedSire(ipoint)) then
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

          if (kb==1) then
              nbisexuals = nbisexuals + 1
              open(BisexUnit,FILE='bisex.txt',position = 'append')
              write(BisexUnit,*) SortedDam(j)
              close(BisexUnit)
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

      if (nbisexuals > 0)  then
        print*, nbisexuals,' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'
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

      !print*,  '  Find seqdam ... '
      do j = 1, nobs
          if (dam(j) == '0') Then
              seqdam(j)=0
          else
              ipoint=INT(ndams/2)
              Noffset = INT(ipoint/2)
              do while (Noffset>1)
                  IF (dam(j).lt.Sorteddam(ipoint)) then
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
              if (kn==1) then
                  seqdam(j) = SorteddamIndex(ipoint)
              else
                  print*, ' Error: dam missing: ', dam(j)
                  stop 1
              end if
          end if
      end do !j

      !print*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'
      !print*,  ' Arranging unlisted base parents ... '

      iextra = newsires + newdams
      if (iextra > 0) then
              ! print*, ' ', iextra, ' unlisted base parents found.'
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

      !print*, ' Inserting unlisted base parents ...'

      oldnobs = nobs
      nobs = nobs + iextra

      !print*, ' Total number of animals = ',nobs

      allocate(passedorder(nobs))
      passedorder=0
      do i = 1+iextra, nobs
          passedorder(i)= i-iextra
          if (sire(i) == '0')then
              seqsire(i) = 0
          else
              seqsire(i) = iextra + seqsire(i)
              if (seqsire(i) > nobs)  seqsire(i) = seqsire(i) - nobs  ! for unlisted sires
          end if

          if (dam(i) == '0') Then
              seqdam(i) = 0
          else
              seqdam(i) = iextra + seqdam(i)
              if (seqdam(i) > nobs)  seqdam(i) = seqdam(i) - nobs
          end if
      end do !i

      do i = 1, newsires
          ID(i) = holdsireid(i)
          passedorder(i)=0
          seqsire(i) = 0
          seqdam(i) = 0
      end do !i

      do i = newsires + 1, newsires + newdams
          ID(i) = holddamid(i - newsires)
          passedorder(i)=0
          seqsire(i) = 0
          seqdam(i) = 0
      end do !i

      deallocate(holdsireid, holddamid, SortedIdIndex, SortedId)

      flag = 0
      Do i = 1, nobs
          if (i <= seqsire(i) .Or. i <= seqdam(i) ) flag = 1
      end do !i
      !if (flag == 0) !print*, 'not needed'!return

      !print*, ' Re-Ordering pedigree ...'
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
          if (seqsire(i) == 0 .And. seqdam(i) == 0) Then
              kn = kn + 1
              NewN(i) = kn
              OldN(kn) = i
          end if
      end do !i

      !Re-order pedigree ...
      NewN(0) = nobs + 1
      flag = 0
      Do While (kn < nobs)
          oldkn = kn
          do i = 1, nobs
              if (NewN(i) == 0) Then !And ID(i) <> 'UniqueNULL' Then
                  Ks = seqsire(i)
                  Kd = seqdam(i)
                  if (NewN(Ks) > 0 .And. NewN(Kd) > 0) Then
                      kn = kn + 1
                      NewN(i) = kn
                      OldN(kn) = i
                  end if
              end if
          end do !i
          ! to avoid hang on unexpected problem ...
          if (kn == oldkn) Then
              flag = flag + 1
          else
              flag = 0
          end if

          if (flag > 10) Then
              open(PedErrUnit,file='ped_err.txt',status='unknown')
              write(PedErrUnit,*) 'Pedigree errors found involving two or more of the following relationships ...'
              write(PedErrUnit,*)
              write(PedErrUnit,*) '       Index numbers are followed by names.'
              write(PedErrUnit,*) '       Index number 0 means unknown, whence name is blank.'
              write(PedErrUnit,*)
              do i = 1, nobs
                  if (NewN(i) == 0) Then
                      write(PedErrUnit,*) 'Individual:',          i, ':  ', ID(i)
                      write(PedErrUnit,*) '    Father:', seqsire(i), ':  ', ID(seqsire(i))
                      write(PedErrUnit,*) '    Mother:',  seqdam(i), ':  ', ID(seqdam(i))
                      write(PedErrUnit,*)
                  end if
              end do !i
              close(PedErrUnit)
              print*,  'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
              stop 1
          end if
      end do !while

      NewN(0) = 0
      do i = 1, nobs
          ID(i) = holdid(OldN(i))
          dooutput(i) = holdoutput(OldN(i))
      end do

      do i = 1, nobs
          seqsire(i) = NewN(holdsire(OldN(i)))
          seqdam(i) = NewN(holddam(OldN(i)))
          if (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
              print*,  'out of order'
              stop 1
          end if
      end do !i

      do i = 1, nobs
          holdsire(i) = passedorder(i)  ! holdsire just because it is free
      end do

      do i = 1, nobs
          passedorder(i) = holdsire(OldN(i))
      end do

      deallocate ( OldN, NewN, holdid, holdsire, holddam) ! holdrec)
      !do i = 1, nobs
      ! print'(3i5,2x,3a4,i5)', i, seqsire(i), seqdam(i), id(i), sire(i), dam(i), passedorder(i)
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

    subroutine AlphaRelateTitles
       print*, ""
       write(*,'(a30,a,a30)') " ","**********************"," "
       write(*,'(a30,a,a30)') " ","*                    *"," "
       write(*,'(a30,a,a30)') " ","*    AlphaRelate     *"," "
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
