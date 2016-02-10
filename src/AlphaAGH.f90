#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#DEFINE BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"

#else
#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#endif

!TODO: write output subroutines and call them in different places

!#########################################################################

module GlobalAGH

    integer,parameter :: lengan=20

    integer :: nTrait,nSnp,nAnisG,nAnisP,nAnisRawPedigree,AllFreqSelCycle,nCols,nAnisH
    integer :: nGMats !,PedigreePresent,WeightYes1No0,
    integer :: GlobalExtraAnimals		!Change John Hickey
		real :: ScaleGtoA  ! For H matrix
    double precision :: DiagFudge
    character(len=20) :: OutputFormat
    character(len=1000) :: AlleleFreqFile

    double precision,allocatable,dimension(:) :: AlleleFreq,Pmat,Adiag
    double precision,allocatable,dimension(:,:) :: Weights,Zmat,tpose,Genos,Amat,InvAmat
    double precision,allocatable,dimension(:,:,:) :: WeightStand,Gmat,InvGmat

		logical :: PedigreePresent,WeightsPresent,ScaleGByRegression
		logical :: Gmake, GInvMake, AMake, AInvMake, HMake, HInvMake  ! What to make
		logical :: GFullMat, GIJA, IGFullMat, IGIJA, AFullMat, AIJA, IAFullMat, IAIJA, HFullMat, HIJA, IHFullMat, IHIJA
		integer :: GType ! InvOpt -- inversion option for G
    !integer :: InvOut,GFullMat,GIJA,IGFullMat,IGIJA,AFullMat,AIJA,IAFullMat,IAIJA,HFullMat,HIJA,IHFullMat,IHIJA
    !integer :: GMake,GInvMake,AMake,AInvMake,HMake,HInvMake,GType

    real(kind=4),allocatable :: xnumrelmatHold(:)
    integer :: NRMmem, shell, shellmax, shellWarning
    integer,allocatable :: seqid(:),seqsire(:),seqdam(:),seqoutput(:),RecodeGenotypeId(:),passedorder(:),RecPed(:,:),dooutput(:)
    character*(lengan),allocatable :: ped(:,:),Id(:),sire(:),dam(:),IdGeno(:)

		! Variable for mapping between Pedigree and genotype animals.
		integer, allocatable :: MapAnimal(:)
		logical, allocatable :: MapToG(:), AnimalsInBoth(:)

end module GlobalAGH

!#########################################################################

program AlphaG

    use GlobalAGH

    implicit none

    include "mkl_vml.f90"

    real :: start, finish

    call cpu_time(start)
    call Titles

    call ReadParam
    call ReadData

    if (AMake) then
        call MakeAMatrix
    endif

    if (AInvMake) then
        call MakeInvAMatrix
    endif

    if (GMake .or. GInvMake) then
        if (GType==1) then
            call MakeGVanRaden
        endif
        if (GType==2) then
            call MakeGNejatiJavaremi
        endif
    endif
    
    if (HInvMake .or. HMake) then
    	call MakeH
    endif


    call cpu_time(finish)
    print *," "
    print '("  Time duration of AlphaAGH = ",f20.4," seconds.")',finish-start
    print *," "

end program AlphaG

!#########################################################################

subroutine ReadParam
    use GlobalAGH
    implicit none

    integer :: i,j,OutputPositions,OutputDigits

    character(len=1000) :: dumC,option,GenotypeFile,PedigreeFile,WeightFile! ,scaleC
    !character(len=1000) :: GFullMatC,GIJAC,InverseGFullMatC,InverseGIJA,AmatFullMat,AmatIJA,InverseAmatFullMat,InverseAmatIJA,HmatFullMat,HmatIJA,InverseHmatFullMat,InverseHmatIJA
    character(len=200) :: OutputPositionsC,OutputDigitsC

    open(unit=11,file="AlphaAGHSpec.txt",status="old")

    read(11,*) dumC			! Input parameters
    read(11,*) dumC,GenotypeFile
    read(11,*) dumC,PedigreeFile
    read(11,*) dumC,WeightFile
    read(11,*) dumC,AlleleFreqFile
    read(11,*) dumC,nTrait
    read(11,*) dumC,nSnp
    read(11,*) dumC,DiagFudge
		read(11,*) dumC,option
    	if (trim(option)=="VanRaden")        GType=1
    	if (trim(option)=="Nejati-Javaremi") GType=2
		read(11,*) dumC,option
		if (trim(option) == 'Regression') then
			ScaleGByRegression = .true.
		else
			ScaleGByRegression = .false.
			read(option,*) ScaleGtoA
		endif
		
    read(11,*) dumC			! Output options

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
		if (GFullMat) GMake = .true.

		read(11,*) dumC, option
		GIJA = trim(option) == 'Yes'
		if (GIJA) GMake = .true.

		! Make inverted G matrix ?
		IGFullMat = .false.
		IGIJA = .false.
		GInvMake = .false.
		read(11,*) dumC, option
		IGFullMat = trim(option) == 'Yes'
		if (IGFullMat) GInvMake = .true.
		
		read(11,*) dumC, option
		IGIJA = trim(option) == 'Yes'
		if (IGIJA) GInvMake = .true.

		! Make A matrix?
		AFullMat = .false.
		AIJA = .false.
		AMake = .false.
		read(11,*) dumC, option
		AFullMat = trim(option) == 'Yes'
		if (AFullMat) AMake = .true.

		read(11,*) dumC, option
		AIJA = trim(option) == 'Yes'
		if (AIJA) AMake = .true.

		! Make inverted A matrix?
		IAFullMat = .false.
		IAIJA = .false.
		AInvMake = .false.
		read(11,*) dumC, option
		IAFullMat = trim(option) == 'Yes'
		if (IAFullMat) AInvMake = .true.

		read(11,*) dumC, option
		IAIJA = trim(option) == 'Yes'
		if (IAIJA) AInvMake = .true.

		! Make H matrix?
		HFullMat = .false.
		HIJA = .false.
		HMake = .false.
		read(11,*) dumC, option
		HFullMat = trim(option) == 'Yes'
		if (HFullMat) HMake = .true.

		read(11,*) dumC, option
		HIJA = trim(option) == 'Yes'
		if (HIJA) HMake = .true.

		! Make inverted H matrix?
		IHFullMat = .false.
		IHIJA = .false.
		HInvMake = .false.
		read(11,*) dumC, option
		IHFullMat = trim(option) == 'Yes'
		if (IHFullMat) HInvMake = .true.

		read(11,*) dumC, option
		IHIJA = trim(option) == 'Yes'
		if (IHIJA) HInvMake = .true.

		if (HMake) then
			GMake = .true.
			AMake = .true.
		endif
		
		if (HInvMake) then
			GMake = .true.
			AMake = .true.
			AInvMake = .true.
		endif



    allocate(AlleleFreq(nSnp))
    allocate(Pmat(nSnp))
    allocate(WeightStand(nSnp,nTrait,nTrait))
    allocate(Weights(nSnp,nTrait))

    if (trim(GenotypeFile)/='None') then
        open(unit=101,file=trim(GenotypeFile),status="old")
    endif

    PedigreePresent=.false.
    if (trim(PedigreeFile)/='None') then
        open(unit=102,file=trim(PedigreeFile),status="old")
        PedigreePresent=.true.
    endif
    if (trim(WeightFile)/='None') then
        open(unit=103,file=trim(WeightFile),status="old")
        WeightsPresent = .true.
    else
        WeightsPresent = .false.
    endif

    nGMats=0
    do i=1,nTrait
        do j=i,nTrait
            nGMats=nGMats+1
        enddo
    enddo

    if (trim(GenotypeFile)=='None') then
    		if (GMake .or. GInvMake .or. HMake .or. HInvMake) print *, 'In order to create G or H matrices, a genotype file must be given.'
        GMake=.false.
        GInvMake=.false.
        HMake=.false.
        HInvMake=.false.
    endif

    if (trim(PedigreeFile)=='None') then
    		if (AMake .or. AInvMake .or. HMake .or. HInvMake) print *, 'In order to create A or H matrices, a pedigree must be given.'
        AMake=.false.
        AInvMake=.false.
        HMake=.false.
        HInvMake=.false.     
    endif

		if (GType==2 .and. trim(AlleleFreqFile) /= 'None') print *, 'The  Nejati-Javaremi approach does not utilise allele frequencies.'
		
		if (ScaleGByRegression) Amake = .true.

end subroutine ReadParam

!#########################################################################

subroutine ReadData
    use GlobalAGH
    implicit none

    integer :: i,j,stat,GenoInPed,fourthColumn
    character(len=1000) :: dumC
    character(lengan), dimension(1:3) :: pedline

    call CountInData

    allocate(Genos(nAnisG,nSnp))
    allocate(Zmat(nAnisG,nSnp))
    allocate(IdGeno(nAnisG))
    !allocate(RecodeIdGeno(nAnisG))

    if (.not. PedigreePresent) then
        allocate(RecPed(0:nAnisG,4))
        nAnisP=nAnisG
        RecPed(:,:)=0
        do i=1,nAnisP
            RecPed(i,1)=i
            RecPed(i,4)=1
        enddo
    else

				!!! Attempt to magically detect whether there are three or four columns:
				read(102, '(a)', iostat=stat) dumC
				if (stat /= 0) then
					print *, 'Problems reading Pedigree file.'
					!print *, stat, dumC
					stop 1
				endif
				rewind(102)
				
				! Test if the line contains three or four columns
				fourthColumn = -99
				read(dumC, *, iostat=stat) pedline, fourthColumn
				if (stat .eq. -1 .or. fourthColumn .eq. -99) then
					nCols = 3
				else
					nCols = 4
				endif
				!print *, nCols, pedline, fourthColumn
				!!! We now know whether there are three or four columns.
				
				
        allocate(Ped(nAnisRawPedigree,4))
        Ped(:,4) = '1'
        do i=1,nAnisRawPedigree
            read(102,*) ped(i,1:nCols)
        enddo
        call PVseq(nAnisRawPedigree,nAnisP)
        
        allocate(RecPed(0:nAnisP,4))

        RecPed(0,:)=0
        RecPed(:,4)=1
        do i=1,nAnisP
            RecPed(i,1)=i
        enddo
        RecPed(1:nAnisP,2)=seqsire(1:nAnisP)
        RecPed(1:nAnisP,3)=seqdam(1:nAnisP)
        RecPed(1:nAnisP,4)=dooutput(1:nAnisP)        
    endif

    do i=1,nAnisG
        read(101,*) IdGeno(i),Genos(i,:)
    enddo

		! These three vectors uses the Pedigree animals as base,
		! i.e. after reordering, the indices for the nth pedigree animal is n.
		allocate(MapAnimal(1:(nAnisP+nAnisG)))
		allocate(MapToG(1:(nAnisP+nAnisG)))
		allocate(AnimalsInBoth(1:nAnisP+nAnisG))
		MapAnimal = 0
		MapToG = .false.
		AnimalsInBoth = .false.
		nAnisH = nAnisP
		!MapAnimal(1:nAnisP) = ( i, i=1, nAnisP ) ! Why doesn't this work??
		do i=1,nAnisP
			MapAnimal(i) = i
		enddo
		
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
                endif
            enddo
            if (GenoInPed==0) then
            		nAnisH = nAnisH + 1
            		MapAnimal(nAnisH) = i
            		MapToG(nAnisH) = .true.
                print*, "Genotyped individual not in pedigree file - ",trim(IdGeno(i))
                !stop
            endif
        enddo
    endif

    if (WeightsPresent) then
        do i=1,nSnp
            read(103,*) dumC,Weights(i,:)
        enddo
    else
        Weights(:,:)=1
    endif

end subroutine ReadData

!#########################################################################

subroutine MakeInvAMatrix
    use GlobalAGH
    implicit none

    integer :: i,m,n,s
    double precision :: Inbreeding(0:nAnisP),Dii(0:nAnisP)
    character(len=1000) :: nChar,fmt
    logical :: AnimToWrite(nAnisP)

    allocate(InvAmat(0:nAnisP,0:nAnisP))

    print*, "Start calculating inbreeding coefficients"
    call dinbreeding(RecPed(0:nAnisP,1),RecPed(0:nAnisP,2),RecPed(0:nAnisP,3),Inbreeding,nAnisP)
    open(unit=202,file="PedigreeBasedInbreeding.txt",status="unknown")
    print*, "End calculating inbreeding coefficients"
    do i=1,nAnisP
        write(202,'(a20,20000f10.5)') Id(i),Inbreeding(i)
    enddo
    close(202)
    print*, "Start making A inverse"

    !Inverse Ordinary Nrm
    Dii=0
    do i=1,nAnisP
        if (RecPed(i,2)==0) then
            Dii(i)=1
        else
            Dii(i)=0.5-(0.25*(Inbreeding(RecPed(i,2))+Inbreeding(RecPed(i,3))))
        endif
    enddo

    InvAmat=0
    do i=1,nAnisP
        InvAmat(i,i)=1/Dii(i)
        if (RecPed(i,2)/=0) then
            InvAmat(i,RecPed(i,2))=InvAmat(i,RecPed(i,2))+(-1/Dii(i))/2
            InvAmat(i,RecPed(i,3))=InvAmat(i,RecPed(i,3))+(-1/Dii(i))/2
            InvAmat(RecPed(i,2),i)=InvAmat(RecPed(i,2),i)+(-1/Dii(i))/2
            InvAmat(RecPed(i,3),i)=InvAmat(RecPed(i,3),i)+(-1/Dii(i))/2

            InvAmat(RecPed(i,2),RecPed(i,2))=InvAmat(RecPed(i,2),RecPed(i,2))+((1/Dii(i))/4)
            InvAmat(RecPed(i,3),RecPed(i,3))=InvAmat(RecPed(i,3),RecPed(i,3))+((1/Dii(i))/4)
            InvAmat(RecPed(i,2),RecPed(i,3))=InvAmat(RecPed(i,2),RecPed(i,3))+((1/Dii(i))/4)
            InvAmat(RecPed(i,3),RecPed(i,2))=InvAmat(RecPed(i,3),RecPed(i,2))+((1/Dii(i))/4)
        endif
    enddo
    print*, "Finished making A inverse"

    if (IAFullMat) then
    		AnimToWrite = RecPed(1:nAnisP,4) == 1
    		s = COUNT(AnimToWrite)
    		
    		write(*,'(a40,i6,a11)') " Start writing A inverse full matrix for", s," individuals"
        write(nChar,*) s
        fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="InvAFullMatrix.txt",status="unknown")
        do m=1,nAnisP
        		if (AnimToWrite(m)) then
        			write(202,fmt) Id(m), pack(InvAmat(m,1:nAnisP), AnimToWrite)
	          endif
        enddo
        close(202)
        print*, "End writing A inverse full matrix"
    endif

    if (IAIJA) then
        print*, "Start writing A inverse ija"
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="InvAija.txt",status="unknown")
        do m=1,nAnisP
            do n=1,m
                if (RecPed(m,4) == 1 .and. RecPed(n,4) == 1 .and. InvAmat(m,n) /= 0) write(202,fmt) Id(m),Id(n),InvAmat(m,n)
            enddo
        enddo
        close(202)
        print*, "End writing A inverse ija"
    endif

		!if (HInvMake) deallocate(InvAmat)

end subroutine MakeInvAMatrix

!#########################################################################

subroutine MakeAMatrix
    use GlobalAGH
    implicit none

    integer :: i,j,k,m,n,s,div
    double precision :: Amatavg
    character(len=1000) :: nChar,fmt
    logical :: AnimToWrite(nAnisP)
		
    allocate(Amat(0:nAnisP,0:nAnisP))

    print*, "Start making A"
    Amat=0
    do i=1,nAnisP
        Amat(i,i)=1+Amat(RecPed(i,2),RecPed(i,3))/2
        do j=i+1,nAnisP
            Amat(i,j)=(Amat(i,RecPed(j,2))+Amat(i,RecPed(j,3)))/2
            Amat(j,i)=Amat(i,j)
        enddo
    enddo
    print*, "Finished making A"

		! Record diagonals of animals in both A and G:
		if (ScaleGByRegression) then
			allocate(Adiag(0:Count(AnimalsInBoth)))
			div = Count(AnimalsInBoth)**2
			Amatavg = 0
			k = 0
			do i = 1,nAnisP
				if (.not. AnimalsInBoth(i)) cycle
				k = k + 1
				Adiag(k) = Amat(i,i)
				do j=1,nAnisP
					if (AnimalsInBoth(j)) Amatavg=Amatavg + Amat(i,j) * 2 / div
				enddo
			enddo
			Adiag(0) = Amatavg
		endif

    if (AFullMat) then
    		AnimToWrite = RecPed(1:nAnisP,4) == 1
    		s = COUNT(AnimToWrite)
    		
    		write(*,'(a32,i6,a11)') " Start writing A full matrix for", s," individuals"
        write(nChar,*) s
        fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="AFullMatrix.txt",status="unknown")
        do m=1,nAnisP
        		if (AnimToWrite(m)) then
        			write(202,fmt) Id(m), pack(Amat(m,1:nAnisP), AnimToWrite)
	          endif
        enddo
        close(202)
        print*, "End writing A full matrix"
    endif

    if (AIJA) then
        write(*,'(a24,i6,a11)') " Start writing A ija for", s," individuals"
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="Aija.txt",status="unknown")
        do m=1,nAnisP
            do n=1,m
            	  if (RecPed(m,4) == 1 .and. RecPed(n,4) == 1 .and. Amat(m,n) /= 0)  write(202,fmt) Id(m),Id(n),Amat(m,n)
            enddo
        enddo
        close(202)
        print*, "End writing A ija"
    endif
    
    !if (HMake /= 1 .and. HInvMake /= 1)  deallocate(Amat)

end subroutine MakeAMatrix

!#########################################################################

subroutine MakeGVanRaden
    use GlobalAGH
    implicit none

    integer :: i,j,k,l,m,n,nMissing,dumI,WhichMat,stat!,errorflag
    double precision :: TmpWt(nSnp),TmpVal(1),Denom
    character(len=1000) :: filout,nChar,fmt

    allocate(Gmat(nAnisG,nAnisG,nGMats))
    allocate(tpose(nSnp,nAnisG))    


    print*, "Start making G - VanRaden"
    
    !! Allele frequencies
    if (trim(AlleleFreqFile) == 'None') then
    	!Calculate Allele Freq
			AlleleFreq(:)=0.0
			do j=1,nSnp
					nMissing=0
					do i=1,nAnisG
							if ((Genos(i,j)>-0.1).and.(Genos(i,j)<2.1)) then
									AlleleFreq(j)=AlleleFreq(j)+Genos(i,j)
							else
									nMissing=nMissing+1
							endif
					enddo
					! Write the frequency of SNP j in array. If all SNPs are missing, then freq_j=0
					if (nAnisG>nMissing) then
							AlleleFreq(j)=AlleleFreq(j)/(2*(nAnisG-nMissing))
					else
							AlleleFreq(j)=0.0
					endif
			enddo
			open(unit=201, file='AlleleFreqTest.txt',status='unknown')
			do j=1,nSnp
					write(201,*) j,AlleleFreq(j)
			enddo
			close(202)			

		else
			! Read allele frequencies from file.
			open(unit=202, file=trim(AlleleFreqFile), status='OLD')
			do i=1,nSnp
				read(202,*,iostat=stat) dumI,AlleleFreq(i)  !AlleleFrequencies are kept in second column to keep consistency with AlphaSim.
				if (stat /= 0) stop "Problems reading allele frequency file."
			enddo
			close(202)
		endif

    !Create Pmat
    do i=1,nSnp
        Pmat(i)=2*(AlleleFreq(i)-0.5)
    enddo

    !Standardise weights
    do i=1,nTrait
        do j=1,nTrait
            if (i==j) then
                TmpVal(1)=sum(Weights(:,j))
                do k=1,nSnp
                    WeightStand(k,i,j)=Weights(k,j)/TmpVal(1)
                enddo
            else
                do k=1,nSnp
                    TmpWt(k)=Weights(k,i)*Weights(k,j)
                enddo
                TmpVal(1)=sum(TmpWt(:))
                do k=1,nSnp
                    WeightStand(k,i,j)=TmpWt(k)/TmpVal(1)
                enddo

            endif
        enddo
    enddo

    !Make Z
    do i=1,nAnisG
        do j=1,nSnp
            if ((Genos(i,j)>-0.1).and.(Genos(i,j)<2.1)) then
                Zmat(i,j)=(Genos(i,j)-1.0)-(Pmat(j))
            else
                Zmat(i,j)=((2*AlleleFreq(j))-1.0)-(Pmat(j))
            endif
        enddo
    enddo


    !Make G matrices
    WhichMat=0
    do i=1,nTrait
        do j=i,nTrait
            WhichMat=WhichMat+1

            !Get Denom
            Denom = 0.000000001
            do k=1,nSnp
                Denom=Denom+(AlleleFreq(k)*(1-AlleleFreq(k))*WeightStand(k,i,j))
            enddo
            Denom = 2*Denom

            tpose=transpose(Zmat)

            !Make ZH
            do k=1,nSnp
                Zmat(:,k)= Zmat(:,k)*WeightStand(k,i,j)
            enddo

            GMat(:,:,WhichMat)=matmul(Zmat,tpose)

            do l=1,nAnisG
                do m=1,l
                    GMat(l,m,WhichMat)=GMat(l,m,WhichMat)/Denom
                    GMat(m,l,WhichMat)=GMat(l,m,WhichMat)
                enddo
            enddo
            do l=1,nAnisG
                GMat(l,l,WhichMat)=GMat(l,l,WhichMat)+DiagFudge
            enddo

            if (GFullMat) then
                write(filout,'("GFullMatrix"i0,"-"i0".txt")')i,j
                write(nChar,*) nAnisG
                fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
                open(unit=202,file=trim(filout),status="unknown")
                do m=1,nAnisG
                    write(202,fmt) IdGeno(m),GMat(m,:,WhichMat)
                enddo
                close(202)
            endif

            if (GIJA) then
                fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
                write(filout,'("Gija"i0,"-"i0".txt")')i,j
                open(unit=202,file=trim(filout),status="unknown")
                do m=1,nAnisG
                    do n=1,m
                        write(202,fmt) IdGeno(m),IdGeno(n),GMat(m,n,WhichMat)
                    enddo
                enddo
                close(202)
            endif

            if (GInvMake) then
                allocate(InvGmat(nAnisG,nAnisG,nGMats))

                print*, "Start inverting G - VanRaden"
								InvGmat(:,:,WhichMat)=Gmat(:,:,WhichMat)
								call invert(InvGmat(:,:,WhichMat),nAnisG,.true., 1)

                print*, "Finished inverting G - VanRaden"

                if (IGFullMat) then
                    write(nChar,*) nAnisG
                    fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
                    write(filout,'("InvGFullMatrix"i0,"-"i0".txt")')i,j
                    open(unit=202,file=trim(filout),status="unknown")
                    do m=1,nAnisG
                        write(202,fmt) IdGeno(m),InvGMat(m,:,WhichMat)
                    enddo
                    close(202)
                endif

                if (IGIJA) then
                    write(filout,'("InvGija"i0,"-"i0".txt")')i,j
                    fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
                    open(unit=202,file=trim(filout),status="unknown")
                    do m=1,nAnisG
                        do n=1,m
                            write(202,fmt) IdGeno(m),IdGeno(n),InvGMat(m,n,WhichMat)
                        enddo
                    enddo
                    close(202)
                endif
            endif

        enddo
    enddo
    deallocate(tpose)

    print*, "Finished making G - VanRaden"

end subroutine MakeGVanRaden

!#########################################################################

subroutine MakeGNejatiJavaremi
    use GlobalAGH
    implicit none

    integer :: i,j,l,m,n!,errorflag
    character(len=1000) :: nChar,fmt

    allocate(Gmat(nAnisG,nAnisG,1))
    allocate(tpose(nSnp,nAnisG))


    print*, "Start making G - Nejati-Javaremi"

    !Make Z
    do i=1,nAnisG
        do j=1,nSnp
            if ((Genos(i,j)>=0).and.(Genos(i,j)<=2)) then
                Zmat(i,j)=Genos(i,j)-1.0
            else
                Zmat(i,j)=0
            endif
        enddo
    enddo

    !Make G matrices
    tpose=transpose(Zmat)
    GMat(:,:,1)=matmul(Zmat,tpose)
    do l=1,nAnisG
        do m=1,l
            GMat(l,m,1)=GMat(l,m,1)/nSnp + 1.0
            GMat(m,l,1)=GMat(l,m,1)
        enddo
        Gmat(l,l,1)=Gmat(l,l,1)+DiagFudge
    enddo

    if (GFullMat) then
        write(nChar,*) nAnisG
        fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="GFullMatrix1-1.txt",status="unknown")
        do m=1,nAnisG
            write(202,fmt) IdGeno(m),GMat(m,:,1)
        enddo
        close(202)
    endif

    if (GIJA) then
        fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
        open(unit=202,file="Gija1-1.txt",status="unknown")
        do m=1,nAnisG
            do n=1,m
                write(202,fmt) IdGeno(m),IdGeno(n),GMat(m,n,1)
            enddo
        enddo
        close(202)
    endif

    if (GInvMake) then
        allocate(InvGmat(nAnisG,nAnisG,1))

        print*, "Start inverting G - Nejati-Javaremi"
				InvGmat(:,:,1)=Gmat(:,:,1)
				call invert(InvGmat(:,:,1),nAnisG,.true.,1)
        print*, "Finished inverting G - Nejati-Javaremi"

        if (IGFullMat) then
            write(nChar,*) nAnisG
            fmt="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
            open(unit=202,file="InvGFullMatrix1-1.txt",status="unknown")
            do m=1,nAnisG
                write(202,fmt) IdGeno(m),InvGMat(m,:,1)
            enddo
            close(202)
        endif

        if (IGIJA) then
            fmt="(2a20,"//trim(adjustl(OutputFormat))//")"
            open(unit=202,file="InvGija1-1.txt",status="unknown")
            do m=1,nAnisG
                do n=1,m
                    write(202,fmt) IdGeno(m),IdGeno(n),InvGMat(m,n,1)
                enddo
            enddo
            close(202)
        endif
    endif
    deallocate(tpose)

    print*, "Finished making G - Nejati-Javaremi"

end subroutine MakeGNejatiJavaremi

!#########################################################################

subroutine MakeH  ! Both H and Hinv
! Feature added by Stefan Hoj-Edwards, aka The Handsome One, February 2016
! Making the Inverse H matrix ala Christensen 2012 requires:
! Scaling G to A22 (subset of A that is covered by G) by linear regression.
! Inverting the scaled G.
! Replacing subset of inverse A with inverted, scaled G.

! Prerequisite and assumptions for this subroutine:
! There is given both a pedigree and genotype file, and there is an overlap
! of animals between two data sets.
! Diagonals of from both A have been collected during MakeA and MakeG,
! as well as average of A22 .
! Gmat is already calculated and loaded in memory.
! Ainv is calculated an loaded into memory.
!
! Further assumes that animals are ordered the same in both A and G.

	use GlobalAGH
	implicit none
	
	integer :: i,j,k,m,p,q,div,t1,t2,whichMat,nboth
	double precision :: Gmatavg, nom, denom, slope, intercept, Gmean, Amean, Hii
	character(len=1000) :: nChar,fmt1, fmt2,filout
  double precision,allocatable :: Gdiag(:), Hrow(:), A22(:,:), A22Inv(:,:), G22(:,:), A11(:,:), A12(:,:), tmp(:,:), Gboth(:,:)
  character*(lengan),allocatable,dimension(:) :: Ids
  logical,allocatable,dimension(:) :: AnimToWrite !, AhasG
  integer,allocatable :: MapToA11(:), MapToA22(:) !Gmap(:), 
   
   
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
		endif
	enddo

	allocate(A22Inv(nBoth,nBoth))
	allocate(MapToA22(nAnisH))
	if (HMake) allocate(A22(nBoth,nBoth))
	
	k = 0
	do i=1,nAnisP
		if (.not. AnimalsInBoth(i)) cycle
		k = k + 1
		MapToA22(i) = k
		m = 0
		do j=1,nAnisP
			if (.not. AnimalsInBoth(j)) cycle
			m = m + 1
			A22Inv(k,m) = Amat(i,j)
		enddo
	enddo
	if (HMake) A22 = A22Inv
	
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
				endif
			enddo
		enddo
	endif


  
  whichMat = 0
  do t1=1,nTrait
  	do t2=t1,nTrait
  		whichMat = whichMat + 1
  		
			write(*, '(" Starting on H matrix "i0" - "i0)') t1, t2
			
			! Collect G22
			allocate(G22(nAnisG,nAnisG))		

			G22 = 0
			do i=1,nAnisG
				do j=1,nAnisG
					nom = Gmat(i,j,whichMat)
					if (i == j) nom = nom - DiagFudge
					G22(i,j) = nom
				enddo
			enddo
			
			if (ScaleGByRegression) then
				allocate(Gdiag(0:nBoth))
				Gdiag=0
				Gmatavg=0
				div=nBoth**2
				!allocate(Gmap(nBoth))
			
				k = 0
				do i=1,nAnisH
					if (.not. AnimalsInBoth(i)) cycle
					k = k+1
					Gdiag(k) = G22(MapAnimal(i),MapAnimal(i))
					do j=1,nAnisH
						if (.not. AnimalsInBoth(j)) cycle
						Gmatavg=Gmatavg + G22(MapAnimal(i),MapAnimal(j))/div
					enddo
				enddo
				Gdiag(0) = Gmatavg	
	
				! Now do simple linear regression
				nom = 0
				denom = 0
				Gmean = sum(Gdiag) / size(Gdiag, 1)
				Amean = sum(Adiag) / size(Adiag, 1)
				do i=0,ubound(Adiag, 1)
					nom = nom + (Adiag(i) - Amean) * (Gdiag(i) - Gmean)
					denom = denom + (Adiag(i) - Amean) ** 2
				enddo
				slope = nom / denom
				intercept = Amean - slope * Gmean
				
				! Scale G
				G22 = slope * G22 + intercept
				!do i=1,nAnisG
				!	G22(i,i) = G22(i,i) + DiagFudge
				!enddo
				print *, 'Scaling of G:'
				write(*, '(a,f7.4,a,f7.4)'), " G* = G x ", slope, " + ", intercept
				deallocate(Gdiag)
			else
				do i=1,nAnisH
					if (.not. MapToG(i)) cycle
					do j=1,nAnisH
						if (.not. MapToG(j)) cycle
						if (AnimalsInBoth(i) .and. AnimalsInBoth(j)) then
							G22(MapAnimal(i),MapAnimal(j)) = ScaleGToA * G22(MapAnimal(i),MapAnimal(j)) + (1 - ScaleGToA) * Amat(i,j)
						endif
					enddo
				enddo
			endif
			
			do i=1,nAnisG
					G22(i,i) = G22(i,i) + DiagFudge
			enddo
	
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
							if (.not. AnimalsInBoth(j)) cycle
							q = q + 1
							Gboth(p,q) = G22(MapAnimal(i),MapAnimal(j))
						enddo
					else
						k = k+1
						m = 0
						MapToA11(i) = k
						do j=1,nAnisP
							if (AnimalsInBoth(j)) then
								A12(k,MapAnimal(j)) = Amat(i,j)
							else
								m = m+1
								A11(k,m) = Amat(i,j)
							endif
						enddo
					endif
				enddo
				

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
				endif
	
				if (HIJA) then
						write(filout,'("Hija"i0,"-"i0".txt")') t1,t2
						fmt2="(a20,a20,"//trim(adjustl(OutputFormat))//")"
						open(unit=204,file=trim(filout),status="unknown")  
				endif
		
				do i=1,nAnisH
					if (AnimToWrite(i) .eq. .false.) cycle
					Hrow = 0
					k = 0
					do j=1,nAnisH
						if (AnimToWrite(j) .eq. .false.) cycle
						k = k + 1			  
						if (MapToG(i)) then
							if (MapToG(j)) then
								Hii = G22(MapAnimal(i),MapAnimal(j))
							else
								Hii = A12(MapToA11(j),MapAnimal(i)) ! Remember to transpose
							endif
						else
							if (MapToG(j)) then
								Hii = A12(MapToA11(i),MapAnimal(j))
							else
								Hii = A11(MapToA11(i),MapToA11(j))
							endif
						endif
						if (IHIJA .and. i .le. j .and. Hii /= 0) write(204,fmt2) Ids(i), Ids(j), Hii
						Hrow(k) = Hii
					enddo
					if (HFullMat) write(202,fmt1) Ids(i),Hrow(:)
				enddo 
		
	
				if (HFullMat) close(202)
				if (HIJA) close(204)  
		
		
				print *, 'End writing H matrices'
				
			endif
	
			if (HInvMake) then 
				print *, 'Start inverting scaled G matrix'
				call invert(G22, size(G22, 1), .true., 1)

				!print *, 'Gw inverted'				
				!write(fmt2, '(i0)') size(G22,1)
				!fmt1="(a8,"//trim(adjustl(fmt2))//"f8.4)"
				!do i=1,size(G22,1)
				!	write(*,fmt1) IdGeno(i), G22(i,:)
				!enddo
				
				!print *, 'A22 inverted'
				!do i=1,size(G22,1)
				!	write(*,fmt1) IdGeno(i), A22Inv(i,:)
				!enddo
				
				!print *, 'Ainv(22)'
				!do i=1,size(G22,1)
				!	j = i+10
				!	write(*,fmt1) Ids(j), InvAmat(j,11:25)
				!enddo
								
				
				print *, 'End inverting scaled G matrix'

				print *, 'Start writing inverted H matrices (full and/or ija)'

				if (IHFullMat) then
						write(filout,'("InvHFullMatrix"i0,"-"i0".txt")') t1,t2
						write(nChar,*) nAnisH
						fmt1="(a20,"//trim(adjustl(nChar))//trim(adjustl(OutputFormat))//")"
						open(unit=202,file=trim(filout),status="unknown")
				endif
	
				if (IHIJA) then
						write(filout,'("InvHija"i0,"-"i0".txt")') t1,t2
						fmt2="(a,' ',a,' ',"//trim(adjustl(OutputFormat))//")"
						open(unit=204,file=trim(filout),status="unknown")  
				endif
	
				do i=1,nAnisH
					if (AnimToWrite(i) .eq. .false.) cycle
					Hrow = 0
					k = 0
					do j=1,nAnisH
						if (AnimToWrite(j) .eq. .false.) cycle
						k = k + 1				
						if (MapToG(i) .and. MapToG(j)) then
							Hrow(k) = G22(MapAnimal(i),MapAnimal(j))
							if (i <= nAnisP .and. j <= nAnisP) Hrow(k) = Hrow(k) + InvAmat(i,j) - A22Inv(MapToA22(i),MapToA22(j))
						elseif (i <= nAnisP .and. j <= nAnisP) then !if (MapToG(i) .eq. .false. .and. MapToG(j) .eq. .false.	) then
							Hrow(k) = InvAmat(i,j)
						endif
						if (IHIJA .and. i .le. j .and. Hrow(k) /= 0) write(204,fmt2) trim(Ids(i)), trim(Ids(j)), Hrow(k) 
					enddo
					if (IHFullMat) write(202,fmt1) Ids(i),Hrow(:)
				enddo 
	
				if (IHFullMat) close(202)
				if (IHIJA) close(204)
				print *, 'End writing inverted H matrices (full and ija)'
 
			endif
			
			deallocate(Hrow)
			deallocate(G22)
		enddo
	enddo
  
  deallocate(Ids)
  
end subroutine MakeH

!#########################################################################

subroutine CountInData
    use GlobalAGH
    implicit none

    integer :: k
    character(len=300) :: dumC

    nAnisG=0
    do
        read(101,*,iostat=k) dumC
        nAnisG=nAnisG+1
        if (k/=0) then
            nAnisG=nAnisG-1
            exit
        endif
    enddo
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
            endif
        enddo
        rewind(102)
        write(*,'(a2,i6,a33)') "   ",nAnisRawPedigree," individuals in the pedigree file"
    endif

end subroutine CountInData

!#########################################################################

subroutine dinbreeding(ianim,isire,idam,f,n)
    !c     inbreeding program from Meuwissen and Luo (1992)
    !c     GSE 24: 305-313
    !ARRAYS START FROM 0.......
    IMPLICIT NONE
    integer :: n,i,is,id,j,k,ks,kd
    integer :: ianim(0:n),isire(0:n),idam(0:n),ped(0:n,3),point(0:n)
    DOUBLE PRECISION :: f(0:n),l(n),d(n),fi,r

    point(:) = 0
    l(:) = 0
    d(:) = 0

    do i = 0,n
        f(i) = 0.
        if (i.gt.0) then
            ped(i,1) = ianim(i)
            ped(i,2) = isire(i)
            ped(i,3) = idam(i)
        endif
    enddo

    f(0) = -1.0
    do i = 1, n
        is = isire(i)
        id = idam(i)
        ped(i,2) = max(is,id)
        ped(i,3) = min(is,id)
        d(i) = 0.5 - 0.25 * (f(is) + f(id))
        if (is.eq.0.or.id.eq.0) then
            f(i) = 0.0
        else if((ped(i-1,2).eq.ped(i,2)).and.(ped(i-1,3).eq.ped(i,3))) then
            f(i) = f(i-1)
        else
            fi = -1.0
            l(i) = 1.0
            j = i

            do while (j.ne.0)
                k = j
                r = 0.5 * l(k)
                ks = ped(k,2)
                kd = ped(k,3)
                if (ks.gt.0) then
                    do while (point(k).gt.ks)
                        k = point(k)
                    enddo
                    l(ks) = l(ks) + r
                    if (ks.ne.point(k)) then
                        point(ks) = point(k)
                        point(k) = ks
                    endif
                    if (kd.gt.0) then
                        do while (point(k).gt.kd)
                            k = point(k)
                        enddo
                        l(kd) = l(kd) + r
                        if (kd.ne.point(k)) then
                            point(kd) = point(k)
                            point(k) = kd
                        endif
                    endif
                endif
                fi = fi + l(j) * l(j) * d(j)
                l(j) = 0.0
                k = j
                j = point(j)
                point(k) = 0
            enddo

            f(i) = fi
        endif
    enddo

end subroutine dinbreeding

!#########################################################################

subroutine PVseq(nObs,nAnisPedigree)
    use GlobalAGH
    implicit none

    character(LEN=lengan), ALLOCATABLE :: holdsireid(:), holddamid(:)
    character(LEN=lengan), ALLOCATABLE :: holdid(:), SortedId(:), SortedSire(:), SortedDam(:)
    character(LEN=lengan)              :: IDhold
    integer, ALLOCATABLE               :: SortedIdIndex(:), SortedSireIndex(:), SortedDamIndex(:)
    integer, ALLOCATABLE               :: OldN(:), NewN(:), holdsire(:), holddam(:), holdoutput(:)
    INTEGER :: mode    ! mode=1 to generate dummy ids where one parent known.  Geneprob->1  Matesel->0
    INTEGER :: i, j, newid, itth, itho, ihun, iten, iunit
    integer :: nsires, ndams, newsires, newdams, nbisexuals, flag
    INTEGER :: iextra, oldnobs, kn, kb, oldkn, ks, kd
    INTEGER :: Noffset, Limit, Switch, ihold, ipoint
    integer :: nObs,nAnisPedigree,verbose
    character(LEN=lengan) :: path

    mode=1

    allocate(id(0:nobs),sire(nobs),dam(nobs),dooutput(nobs),seqid(nobs),seqsire(nobs),seqdam(nobs),seqoutput(nobs))

    do i=1,nobs
        id(i)=ped(i,1)
        sire(i)=ped(i,2)
        dam(i)=ped(i,3)
        read(Ped(i,4), '(i)') dooutput(i) 
    enddo      
      

    nAnisPedigree=nObs
    path=".\"

    Verbose=1

    do j = 1, nobs
        If (dam(j) == ''.or. dam(j) == '0'.or. dam(j) == '#'.or. dam(j) == '*' .or. dam(j) == '.') Then
            dam(j) = '0'
            seqdam(j)=0
        endif
        If (sire(j) == ''.or.sire(j) == '0'.or.sire(j) == '#'.or.sire(j) == '*'.or.sire(j) == '.') Then
            sire(j) = '0'
            seqsire(j)=0
        endif
    enddo !j

    if(mode.eq.1) then
        !PRINT*,  ' Inserting dummy IDs ... '
        newid=0
        do j = 1, nobs
            if(((sire(j) == '0').and.(dam(j).ne.'0'))  .or. ((sire(j).ne.'0').and.(dam(j) == '0'))) then
                newid=newid+1
                if(newid.gt.99999) then
                    !         PRINT*, newid, ' ...'
                    stop 'too many dummy single parent IDs'
                endif
                itth=int(newid/10000)
                itho=int(newid/1000)-10*itth
                ihun=int(newid/100)-10*itho-100*itth
                iten=int(newid/10)-10*ihun-100*itho-1000*itth
                iunit=newid-10*iten-100*ihun-1000*itho-10000*itth
                if(sire(j) == '0') sire(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
                if( dam(j) == '0')  dam(j)='dum'//achar(48+itth)//achar(48+itho)//achar(48+ihun)//achar(48+iten)//achar(48+iunit)
            endif
        enddo
    endif

    !PRINT*,  ' Sorting Sires ... '
    ALLOCATE(SortedId(nobs), SortedIdIndex(nobs))
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
                endif
            enddo
            Limit = Switch - Noffset
        enddo
        Noffset = INT(Noffset/2)
    enddo
		
		! Count number of unique sires
    nsires=0
    IF(SortedId(1) /= '0') nsires=1
    do i=2,nobs
        IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nsires=nsires+1
    enddo

		! Collect vector of unique sires in sorted order
    ALLOCATE  (SortedSire(0:nsires), SortedSireIndex(nsires))
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
    enddo

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
                endif
            enddo
            Limit = Switch - Noffset
        enddo
        Noffset = INT(Noffset/2)
    enddo

    nDams=0
    IF(SortedId(1) /= '0') nDams=1
    do i=2,nobs
        IF(SortedId(i) /= SortedId(i-1) .and. SortedId(i) /= '0') nDams=nDams+1
    enddo

    ALLOCATE  (SortedDam(0:nDams), SortedDamIndex(ndams))
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
    enddo

    !PRINT*,  ' Sorting IDs ... '
    SortedId(1:nobs) = ID(1:nobs)
    do i=1,nobs
        SortedIdIndex(i) = i
    enddo
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
                endif
            enddo
            Limit = Switch - Noffset
        enddo
        Noffset = INT(Noffset/2)
    enddo

    !PRINT*,  ' Check for duplicate IDs ... '
    flag = -1
    Do i = 2, nobs
        If (SortedID(i) == SortedID(i - 1)) Then
            If (flag == -1) Then
                open (1,FILE='ID_err.txt',STATUS = 'unknown')
                WRITE(1,*) 'Duplicated IDs ...'
                flag = 0
            endif
            WRITE(1,*) SortedID(i)
            flag = flag + 1
        endif
    enddo

    If (flag > -1) Then
        Close (1)
    ! PRINT*, flag,' case(s) of duplicate ID. See ID_ERR.TXT        <------------ WARNING !!!'
    endif

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
            endif
        enddo

        kn=0
        if (SortedSire(j)==SortedId(ipoint)) kn=1
        do while (ipoint<nobs .and. kn==0 .and. SortedSire(j) > SortedId(ipoint))
            ipoint=ipoint+1
        enddo
        if (SortedSire(j)==SortedId(ipoint)) kn=1
        do while (ipoint>1 .and. kn==0 .and. SortedSire(j) < SortedId(ipoint))
            ipoint=ipoint-1
        enddo

        if (SortedSire(j)==SortedId(ipoint)) kn=1
        IF(kn==1) then
            SortedSireIndex(j) = SortedIdIndex(ipoint)
        else    ! sire is unlisted base sire
            newsires = newsires + 1
            SortedSireIndex(j) = nobs + newsires ! for now
        endif
    enddo !j

    ALLOCATE  (holdsireid(newsires))
    kn=0
    do j=1,nsires
        if (SortedSireIndex(j) > nobs) then
            kn=kn+1
            holdsireid(SortedSireIndex(j)-nobs) = SortedSire(j)
        endif
    enddo
    IF(kn /= newsires) stop'newsires error'

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
                endif
            enddo
            kn=0
            if (Sire(j)==SortedSire(ipoint)) kn=1
            do while (ipoint<nsires .and. kn==0 .and. Sire(j) > SortedSire(ipoint))
                ipoint=ipoint+1
            enddo
            if (Sire(j)==SortedSire(ipoint)) kn=1
            do while (ipoint>1 .and. kn==0 .and. Sire(j) < SortedSire(ipoint))
                ipoint=ipoint-1
            enddo
            if (Sire(j)==SortedSire(ipoint)) kn=1
            IF(kn==1) then
                seqsire(j) = SortedSireIndex(ipoint)
            else
                !PRINT*, ' Error: Sire missing: ', Sire(j)
                stop
            endif
        endif
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
            endif
        enddo
        kn=0
        if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint  ! store ipoint here as ipoint can change with bisexuals
        do while (ipoint<nobs .and. kn==0 .and. Sorteddam(j) > SortedId(ipoint))
            ipoint=ipoint+1
        enddo

        if (Sorteddam(j)==SortedId(ipoint)) kn=ipoint
        do while (ipoint>1 .and. kn==0 .and. Sorteddam(j) < SortedId(ipoint))
            ipoint=ipoint-1
        enddo

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
            endif
        enddo

        kb=0
        if (SortedDam(j)==SortedSire(ipoint)) kb=1
        do while (ipoint<nsires .and. kb==0 .and. SortedDam(j) > SortedSire(ipoint))
            ipoint=ipoint+1
        enddo

        if (SortedDam(j)==SortedSire(ipoint)) kb=1
        do while (ipoint>1 .and. kb==0 .and. SortedDam(j) < SortedSire(ipoint))
            ipoint=ipoint-1
        enddo
        if (SortedDam(j)==SortedSire(ipoint)) kb=1

        IF(kb==1) then
            nbisexuals = nbisexuals + 1
            open (1,FILE='bisex.txt',position = 'append')
            WRITE(1,*) SortedDam(j)
            close(1)
        endif

        if (kb==1) then
            SorteddamIndex(j) = SortedSireIndex(ipoint)
        elseif (kn>=1) then
            SorteddamIndex(j) = SortedIdIndex(kn)
        else    ! dam is unlisted base dam
            newdams = newdams + 1
            SorteddamIndex(j) = nobs + newsires + newdams ! for now
        endif
    enddo !j

    If (nbisexuals > 0)  PRINT*, nbisexuals,' bisexual parent(s) found. See file bisex.txt.  <------------ WARNING !!!'

    ALLOCATE  (holddamid(newdams))

    kn=0
    do j=1,ndams
        if (SortedDamIndex(j) > nobs+newsires) then
            kn=kn+1
            holddamid(SortedDamIndex(j)-nobs-newsires) = SortedDam(j)
        endif
    enddo

    IF(kn /= newdams) stop'newdams error'

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
                endif
            enddo
            kn=0
            if (dam(j)==Sorteddam(ipoint)) kn=1
            do while (ipoint<ndams .and. kn==0 .and. dam(j) > Sorteddam(ipoint))
                ipoint=ipoint+1
            enddo
            if (dam(j)==Sorteddam(ipoint)) kn=1
            do while (ipoint>1 .and. kn==0 .and. dam(j) < Sorteddam(ipoint))
                ipoint=ipoint-1
            enddo
            if (dam(j)==Sorteddam(ipoint)) kn=1
            IF(kn==1) then
                seqdam(j) = SorteddamIndex(ipoint)
            else
                ! PRINT*, ' Error: dam missing: ', dam(j)
                stop
            endif
        endif
    ENDDO !j

    !PRINT*,  '  Dams: ',newdams,' unlisted, ',ndams,' in total'
    !PRINT*,  ' Arranging unlisted base parents ... '

    iextra = newsires + newdams
    If (iextra > 0) then
            ! PRINT*, ' ', iextra, ' unlisted base parents found.'
        ! SortedId and SortedIdIndex just used as a holder while redimensioning
        SortedId(1:nobs)=id(1:nobs)
        deallocate (id)
        ALLOCATE(id(nobs+iextra))
        id(1+iextra:nobs+iextra)=SortedId(1:nobs)
        
        ! Do sires
        SortedId(1:nobs)=sire(1:nobs)
        deallocate (sire)
        ALLOCATE(sire(nobs+iextra))
        sire(1+iextra:nobs+iextra)=SortedId(1:nobs)
        SortedIdIndex(1:nobs)=seqsire(1:nobs)
        deallocate (seqsire)
        ALLOCATE(seqsire(nobs+iextra))
        seqsire(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)
				
        ! Do dams
        SortedId(1:nobs)=dam(1:nobs)
        deallocate (dam)
        ALLOCATE(dam(nobs+iextra))
        dam(1+iextra:nobs+iextra)=SortedId(1:nobs)
        SortedIdIndex(1:nobs)=seqdam(1:nobs)
        deallocate (seqdam)
        ALLOCATE(seqdam(nobs+iextra))
        seqdam(1+iextra:nobs+iextra)=SortedIdIndex(1:nobs)
        
        ! Do output switch, but use seqoutput as placeholder instead of SortedID
        !allocate(seqoutput(1:nobs))
        seqoutput(1:nobs) = dooutput(1:nobs)
        deallocate(dooutput)
        allocate(dooutput(nobs+iextra))
        dooutput(1+iextra:nobs+iextra)=seqoutput(1:nobs)
        if (nCols .eq. 4) dooutput(1:iextra) = 0
        if (nCols .eq. 3) dooutput(1:iextra) = 1
    endif

    !PRINT*, ' Inserting unlisted base parents ...'

    oldnobs = nobs
    nobs = nobs + iextra

    !PRINT*, ' Total number of animals = ',nobs

    ALLOCATE (passedorder(nobs))
    passedorder=0
    do i = 1+iextra, nobs
        passedorder(i)= i-iextra
        If (sire(i) == '0')then
            seqsire(i) = 0
        Else
            seqsire(i) = iextra + seqsire(i)
            If (seqsire(i) > nobs)  seqsire(i) = seqsire(i) - nobs  ! for unlisted sires
        endif

        If (dam(i) == '0') Then
            seqdam(i) = 0
        Else
            seqdam(i) = iextra + seqdam(i)
            If (seqdam(i) > nobs)  seqdam(i) = seqdam(i) - nobs
        endif
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

    DEALLOCATE(holdsireid, holddamid, SortedIdIndex, SortedId)

    flag = 0
    Do i = 1, nobs
        If (i <= seqsire(i) .Or. i <= seqdam(i) ) flag = 1
    enddo !i
    !If (flag == 0) !PRINT*, 'not needed'!return

    !PRINT*, ' Re-Ordering pedigree ...'
    Allocate ( OldN(0:nobs), NewN(0:nobs) )
    ALLOCATE ( holdid(0:nobs), holdsire(nobs), holddam(nobs), holdoutput(nobs) )

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
        endif
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
                endif
            endif
        enddo !i
        ! to avoid hang on unexpected problem ...
        If (kn == oldkn) Then
            flag = flag + 1
        Else
            flag = 0
        endif

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
                endif
            ENDDO !i
            Close (1)
            PRINT*,  'Logical error when re-ordering pedigree - see details in file PED_ERR.TXT'
            stop
        endif
    ENDDO !while

    NewN(0) = 0
    do i = 1, nobs
        ID(i) = holdid(OldN(i))
        dooutput(i) = holdoutput(OldN(i))
    enddo

    do i = 1, nobs
        seqsire(i) = NewN(holdsire(OldN(i)))
        seqdam(i) = NewN(holddam(OldN(i)))
        If (i <= NewN(holdsire(OldN(i))) .Or. i <= NewN(holddam(OldN(i)))) then
            !PRINT*,  'out of order'
            stop
        endif
    ENDDO !i

    DO i = 1, nobs
        holdsire(i) = passedorder(i)  ! holdsire just because it is free
    enddo

    DO i = 1, nobs
        passedorder(i) = holdsire(OldN(i))
    enddo

    deallocate ( OldN, NewN, holdid, holdsire, holddam) ! holdrec)
    !do i = 1, nobs
    ! PRINT'(3i5,2x,3a4,i5)', i, seqsire(i), seqdam(i), id(i), sire(i), dam(i), passedorder(i)
    !enddo

    nAnisPedigree=nObs
    GlobalExtraAnimals=iextra 	!Change John Hickey

end subroutine PVseq

!#############################################################################

!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
    IMPLICIT NONE
    !Declarations
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
    REAL*8, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
    REAL*8, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

    LOGICAL :: FLAG = .TRUE.
    INTEGER :: i, j, k
    REAL*8 :: m
    REAL*8, DIMENSION(n,2*n) :: augmatrix !augmented matrix

    !Augment input matrix with an identity matrix
    DO i = 1, n
        DO j = 1, 2*n
            IF (j <= n ) THEN
                augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
                augmatrix(i,j) = 1
            Else
                augmatrix(i,j) = 0
            ENDIF
        enddo
    enddo

    !Reduce augmented matrix to upper traingular form
    DO k =1, n-1
        IF (augmatrix(k,k) == 0) THEN
            FLAG = .FALSE.
            DO i = k+1, n
                IF (augmatrix(i,k) /= 0) THEN
                    DO j = 1,2*n
                        augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                    enddo
                    FLAG = .TRUE.
                    EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
                    PRINT*, "Matrix is non - invertible"
                    inverse = 0
                    errorflag = -1
                    return
                ENDIF
            enddo
        ENDIF
        DO j = k+1, n
            m = augmatrix(j,k)/augmatrix(k,k)
            DO i = k, 2*n
                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            enddo
        enddo
    enddo

    !Test for invertibility
    DO i = 1, n
        IF (augmatrix(i,i) == 0) THEN
            PRINT*, "Matrix is non - invertible"
            inverse = 0
            errorflag = -1
            return
        ENDIF
    enddo

    !Make diagonal elements as 1
    DO i = 1 , n
        m = augmatrix(i,i)
        DO j = i , (2 * n)
            augmatrix(i,j) = (augmatrix(i,j) / m)
        enddo
    enddo

    !Reduced right side half of augmented matrix to identity matrix
    DO k = n-1, 1, -1
        DO i =1, k
            m = augmatrix(i,k+1)
            DO j = k, (2*n)
                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            enddo
        enddo
    enddo

    !store answer
    DO i =1, n
        DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
        enddo
    enddo
    errorflag = 0
END SUBROUTINE FINDinv

!#########################################################################

subroutine invert(x,n,sym, method)

  ! Interface to call inverse subroutines from BLAS/LAPACK libraries

  ! x symmetric positive-definite matrix to be inverted
  ! n matrix dimension
  ! sym return lower-triangular (sym=.false) or full matrix (sym=.true.)
  ! method for inversion
  ! 0 -- Generalised solving using LU decomposition (dgetrs)
  ! 1 -- Cholesky decomposition

  implicit none
  integer, intent(in) :: n,method
  logical, intent(in) :: sym
  integer :: i,j,info
  double precision,intent(inout) :: x(n,n)
  double precision,allocatable :: Iden(:,:)

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
			stop
		endif
		
		x(:,:) = Iden(:,:)
		
	else if (method == 1) then

		! Computes the Cholesky factorization of a symmetric positive definite matrix
		! https://software.intel.com/en-us/node/468690
		call dpotrf('L',n,x,n,info)
		if (info /= 0) then
			print*,'Matrix not positive-definite - info',info
			stop
		endif

		! Computes the inverse of a symmetric positive definite matrix,
		!   using the Cholesky factorization computed by dpotrf()
		! https://software.intel.com/en-us/node/468824
		call dpotri('L',n,x,n,info)
		if (info /= 0) then
		 print*,'Matrix not positive-definite - info',info
		 stop
		endif

		! Fills the upper triangle
		if (sym) then
			forall (i=1:n,j=1:n,j>i) x(i,j)=x(j,i)
		endif

	endif

end subroutine invert

!#########################################################################

subroutine Titles
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
end subroutine Titles

!#########################################################################
