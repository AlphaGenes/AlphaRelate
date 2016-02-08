# General variables
NAME:=AlphaAGH
VERSION:= $(shell git rev-parse --short HEAD)
SUBVERSION:=0
PROGRAM:=$(NAME)$(VERSION).$(SUBVERSION)

# Set the default compiler to iFort
FC:=ifort
FFLAGS:=-O3 -DVERS=""commit-$(VERSION)""

# Hello friends,
# If compiling with MKL fails -- or succeeds and running the binary fails --
# try running 
# source /opt/intel/composer_xe_2015.3.187/mkl/bin/mklvars.sh intel64
# with the path replace as necessary...

#  If -D WEB is specified, stops will be put into AlphaAGH.

# MS Windows
ifeq ($(OS), Windows_NT)
	SRCDIR      := src/
	BUILDDIR    :=
	TARGETDIR   :=
	OSFLAG := "OS_WIN"
	# TODO: What is the MKL path on Windoze?
	MKLROOT := #???
	MKLLIB := -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
	MKLINC := -I$(MKLROOT)/include
	FFLAGS := $(FFLAGS) /static /i8 /fpp  /Qmkl /D $(OSFLAG) $(MKLLIB) $(MKLINC)
	ABOPT := -static  -Qmkl # not used!
	obj := .obj

	MAKEDIR :=
	exe := .exe
	CC := cl
	CFLAGS := /EHsc

	DEL := del
else
	# Linux or Mac OSX
	SRCDIR      := src/
	BUILDDIR    := objs/
	TARGETDIR   := bin/
	obj := .o
	OSFLAG := "OS_UNIX"
	# TODO: can we make this generic?
	MKLROOT := /opt/intel/mkl
	# On Eddie
	# MKLROOT:=/exports/applications/apps/intel/ClusterStudio2013/mkl
	MKLLIB := -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm
	MKLINC := -I$(MKLROOT)/include
	ABOPT := -mkl -static-intel -openmp-link=static # not used!
	exe :=
	FFLAGS:= $(FFLAGS) -mkl -i8 -static-intel -fpp -openmp-link=static  -module $(BUILDDIR) -D $(OSFLAG) $(MKLLIB) $(MKLINC)
	uname := $(shell uname)
	MAKEDIR := @mkdir -p
	DEL := rm -rf
	# Linux only
	ifeq ($(uname), Linux)
		FFLAGS := $(FFLAGS) -static -static-libgcc -static-libstdc++
	endif
endif

# Compile everything
all: directories $(TARGETDIR)$(NAME)$(exe) $(TARGETDIR)AlphaAGH$(exe)

directories:
	$(MAKEDIR)  $(TARGETDIR)
	$(MAKEDIR)  $(BUILDDIR)

# Compilation options for debugging
# With warnings about not used variables
debuglong: FFLAGS:= -i8 -traceback -g -debug all -fpp -ftrapuv -module $(BUILDDIR) -fpe0 -warn -check all -D $(OSFLAG)

debuglong: all

# With memory checks
debug: FFLAGS:= $(FFLAGS) -i8 -traceback -g -D VERS=""commit-$(VERSION)"" -D $(OSFLAG) -debug all -warn -check bounds -check format \
		-check output_conversion -check pointers -check uninit -fpp -module $(BUILDDIR)

debug: all

web: FFLAGS := $(FFLAGS) -D "WEB"

web: all

# If binary is made, intermediate files will be binary
binary: FFLAGS := $(FFLAGS) -D "BINARY"

binary: all

# Compile AlphaAGH
$(TARGETDIR)AlphaAGH$(exe): $(SRCDIR)AlphaAGH.f90
	@echo "Compiling AlphaAGH..."
	$(FC) $(SRCDIR)AlphaAGH.f90 $(FFLAGS) -o $(TARGETDIR)AlphaAGH$(exe)
	@echo

# Cleaning
sparklinglyclean: veryclean
	$(DEL) TARGETDIR

veryclean: clean
	$(DEL) $(TARGETDIR)AlphaAGH$(exe)

clean:
	$(DEL) $(BUILDDIR) *$(obj) *.mod *.dwarf *.i90 *__genmod* *~

.PHONY: make veryclean all
