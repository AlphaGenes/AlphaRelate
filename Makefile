# General variables
NAME:=AlphaRelate
VERSION:=$(shell git rev-parse --short HEAD)
SUBVERSION:=0
PROGRAM:=$(NAME)$(VERSION).$(SUBVERSION)
ALPHAHOUSEDIR:=../AlphaHouse/

# Set the default compiler to iFort
FC:=ifort
FFLAGS:=-O3 -DVERS=""commit-$(VERSION)""

# If -D WEB is specified, stops will be put into the binary

# MS Windows
ifeq ($(OS), Windows_NT)
	SRCDIR:=src/
	BUILDDIR:=
	TARGETDIR:=
	MAKEDIR:=
	TESTDIR:=
	PFUNIT:=
	OSFLAG:="OS_WIN"
	# TODO: What is the MKL path on Windoze?
	MKLROOT:=#???
	MKLLIB:=-L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
	MKLINC:=-I$(MKLROOT)/include
	## see also https://software.intel.com/en-us/compiler_winapp_f (2014-12-03)
	FFLAGS:=$(FFLAGS) /fpp /Qmkl /Qopenmp /static /Qopenmp-link:static /Qlocation,link,"${VCINSTALLDIR}/bin" /D $(OSFLAG) $(MKLINC) $(MKLLIB)
	obj:=.obj
	exe:=.exe
	CC:=cl
	CFLAGS:=/EHsc
	DEL:=del
else
	# Linux or Mac OSX
	SRCDIR:=src/
	BUILDDIR:=objs/
	TARGETDIR:=bin/
	MAKEDIR:=@mkdir -p
	TESTDIR:=tests/
	PFUNIT:=/usr/local/pFUnit_serial
	OSFLAG:="OS_UNIX"
	# On Mac
	MKLROOT:=/opt/intel/mkl
	# On Eddie2
	# MKLROOT:=/exports/applications/apps/intel/ClusterStudio2013/mkl
	# On Eddie3
	# MKLROOT:=/exports/applications/apps/SL7/intel/parallel_studio_xe_2016/mkl
	MKLLIB:=-L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
	MKLINC:=-I$(MKLROOT)/include
	obj:=.o
	exe:=
	DEL:=rm -rf
	FFLAGS:=$(FFLAGS) -fpp -mkl -static-intel -module $(BUILDDIR) -D $(OSFLAG) $(MKLINC) $(MKLLIB)
	OMPFLAGS:=-openmp -openmp-link=static
	uname:=$(shell uname)
	# Linux only
	ifeq ($(uname), Linux)
		FFLAGS:=$(FFLAGS) -static -static-libgcc -static-libstdc++ $(OMPFLAGS)
	endif
endif

# Required modules
MODS:=$(ALPHAHOUSEDIR)AlphaHouseMod.f90

FILES:=$(SRCDIR)AlphaRelateMod.f90 $(SRCDIR)AlphaRelate.f90

#Get all the pf files

PFTESTS:=$(wildcard $(TESTDIR)*.pf)

F90TESTS:=$(PFTESTS:.pf=.F90)

#Get the current date and time (for git)
BUILDDATE:=$(shell date +%Y%m%d-%H:%M:%S)

# Compile everything
all: directories $(TARGETDIR)$(NAME)$(exe) preprocess_test test_clean

all_git: directories $(TARGETDIR)$(NAME)$(exe) git preprocess_test test_clean

directories:
	$(MAKEDIR) $(TARGETDIR)
	$(MAKEDIR) $(BUILDDIR)

# Compilation options for debugging
# With warnings about not used variables
debuglong: FFLAGS:=$(FFLAGS) -traceback -g -debug all -ftrapuv -fpe0 -warn -check all

debuglong: all

# With memory checks
debug: FFLAGS:=$(FFLAGS) -traceback -g -debug all -warn -check bounds -check format \
		-check output_conversion -check pointers -check uninit

debug: all

web: FFLAGS:=$(FFLAGS) -D "WEB"

web: all

# If binary is made, intermediate files will be binary
binary: FFLAGS:=$(FFLAGS) -D "BINARY"

binary: all

# Compile
$(TARGETDIR)$(NAME)$(exe): Makefile $(MODS) $(SRCDIR)$(NAME).f90
	@echo "Compiling $(NAME)..."
	$(FC) $(MODS) $(FILES) $(FFLAGS) -o $(TARGETDIR)$(NAME)$(exe)
	@echo

# Cleaning
sparklinglyclean: veryclean
	$(DEL) TARGETDIR

veryclean: clean
	$(DEL) $(TARGETDIR)$(NAME)$(exe)

clean:
	$(DEL) $(BUILDDIR) *$(obj) *.mod *.dwarf *.i90 *__genmod* *~

preprocess_test: $(F90TESTS)
	cd $(BUILDDIR)
	cp $(TES)testSuites.inc .
	$(FC) -o tests.x $^ -I$(THIS_DIR)$(BUILDDIR) $(PFUNIT)/include/driver.f90 -I$(PFUNIT)/mod -lpfunit -module $(BUILDDIR)
	./tests.x

git: $(FILES) $(PFTESTS)
	git add $^
	git commit --allow-empty -m "Successful compilation $(BUILDDATE)"

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@

test_clean:
	$(RM) $(THIS_DIR)$(TESTDIR)*.F90
	$(RM) $(BUILDDIR)tests.x
	$(RM) $(BUILDDIR)testSuites.inc

.PHONY: make veryclean all
