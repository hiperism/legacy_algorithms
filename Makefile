#! /bin/csh -f
## Copyright 2019 HiPERiSM Consulting, LLC, 08-01-2019
#
##+++++++++ The following is for the Intel 32-core  ++++++++++
##++++++++++ with Intel ifort compiler 17.6    +++++++++++++++++++++
## ======================== Sophie Build Script ================= 
# Usage: make > make.log 2>&1 
#
# Requirements: and Intel Fortran compiler
#               
#          now using ifort v17.6 on node20
 FC    = /opt/intel17/compilers_and_libraries_2017.6.256/linux/bin/intel64/ifort
#CC    = /opt/intel17/compilers_and_libraries_2017.6.256/linux/bin/intel64/icc
#FP   = $FC
#
## /opt/intel17/compilers_and_libraries_2017.6.256/linux/mpi/intel64/bin/mpiifort
## /opt/intel17/compilers_and_libraries_2017.6.256/linux/bin/intel64/icpc
## /opt/intel17/compilers_and_libraries_2017.6.256/linux/mpi/intel64/bin/mpiicc
## /opt/intel17/compilers_and_libraries_2017.6.256/linux/mpi/intel64/bin/mpiicpc
## /opt/intel17/compilers_and_libraries_2017.6.256/linux/mpi/intel64/bin/mpirun
## vendor library is Intel BLAS.
##
## for sequential compilation use  the compiler option "-lmkl=sequential" and link for LP64 (ILP64 requires "-i8" for large integers)

 MKLROOT         = /opt/intel17/compilers_and_libraries_2017.6.256/linux/mkl/
 MKLPATH         = $(MKLROOT)/lib/intel64_lin/
 MKLINCLUDE      = $(MKLROOT)/include
#INTEL_BLAS      = -L$(MKLPATH) -I$(MKLINCLUDE) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 INTEL_BLAS      = -L$(MKLPATH) -I$(MKLINCLUDE) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

## in fortran dynamic linking and sequential verison of Intel MKL supporting LP64 interface:
## ifort myprog.f -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

BLASLIB         = $(INTEL_BLAS)

#
# The following is commentary on setting compiler flags for Intel ifort
#
#=>> for Intel some optimization flags are as follows:
#
#  -xHost               [optimize for IA-64 processor on host]
#  -O[n]                [optimization levels n=0, 1, 2, 3]
#  -inline-level=2      [Enables inlining of  any  function
#                             at the compiler's discretion]
#  -fast                [enables -ipo, -O3, and -static
#                        use is avoided because of -static flag ]
#  -ftz                 [On systems using Intel 64 architecture, the compiler
#                               denoraml results are flushed to zero]
#  -ip                  [Enables ipo for single file compilation]
#  -ipo                 [Enables ipo between files ]
#  -openmp, -parallel,  [thread parallel mode]
#  -prof-gen, -prof-use [profile guided optimization - requires runtime data]
#  -fnsplit             [function splitting is only enabled with -prof-use]
#
# compiler reporting options are as follows
#
# -opt-report-file=name [compiler generates an optimization report to a file]
# -opt-report-phase=all [report on all optimization phases]
#
#=>> select standard options here
#
# FSTD     = -fixed -extend_source 132 -cm -V -v -fno-alias
# FSTD     = -V -v -traceback -warn -noextend-source -nowarn nodeclarations
  FSTD     = -V -v -traceback -warn -noextend-source -nowarn 
#
# On intel 64 processors the following is the default and need not be set
# 
 Farch    = -mavx -qoverride_limits
#
#=>> choose the memory model level from ONLY ONE of these choices
#    "small" is the default, "medium" or "large" REQUIRE use of shared
#    libraries because static libraries are not allowed for these choices 
#
# Fmem      = -mcmodel=small -m64
  Fmem      = -shared-intel -mcmodel=medium 
# Fmem      = -shared-intel -mcmodel=large  -m64
#
#=>> choose the optimization level from ONLY ONE of these choices
#
# Fopt     = -O0
# Fopt     = -O1
# Fopt     = -O2
# Fopt     = -O2 -ftz -align all -fno-alias
  Fopt     = -O3 -ftz -align all -fno-alias
# 
#=>> inline and ipo
#
# Finline  = -ip-no-inlining -Winline
  Finline  = -inline-level=2
  Fipo     = -ip -ipo
#
#=>> profiling with gprof
#
  Fprof    = -p
#
#=>> bufferred IO may help performance
#
  Fbuff    = -assume buffered_io
#
#=>> optimization reports requires name of a file to follow
#
  Frpt     = -qopt-report=5 -qopt-report-phase=loop,openmp,vec -qopt-report-file=report.txt
#
#=>> these are for profile guided optimization [consult the documentation]
#    their use is not recommended since use is complex and yields little
#
# Fpgog    =  -prof-gen
# Fpgou    =  -prof-use
#
#
#=>> Arithmetic contstraints may reduce performance
#
  Fmath    = -mieee-fp -init=arrays,zero -double-size 128 -real-size 128 -integer-size 32
#
#
#=>> FOR mkl 
#
  Fmkl     = -mkl=sequential
#
#
#=>> For opennmp 
#
# Fomp     = -qopenmp -qopenmp-report1
  Fomp     = 
#
#=>> for MPI code using this in both compile and link flags
#
# Fmpi     = -I/opt/intel_compiler/impi/5.0.2.044/intel64/include
  Fmpi     = 
#
#
#=>> For Debug
#
#  DebugFLAGS = -g -O0 -check bounds -fpe-all=0 -fpe0
#
#
#=>> select preferences below by removing the '#' on only ONE of the F_FLAGS lines 
#
## optimization flags - basic
# F_FLAGS = $(FSTD) $(Farch) $(Fmem) $(Fopt) $(Fmath) $(Frpt) $(Fmkl)
#
## optimization + gprof profile flags
# F_FLAGS = $(FSTD) $(Farch) $(Fmem) $(Fopt) $(Finline) $(Fipo) $(Fmath) $(Frpt) $(Fmkl) $(Fprof)
#
## optimization + buffered io flags
# F_FLAGS = $(FSTD) $(Farch) $(Fmem) $(Fopt) $(Finline) $(Fipo) $(Fmath) $(Frpt) $(Fmkl) $(Fbuff)
#
## optimization + ipo + inline  flags
  F_FLAGS = $(FSTD) $(Farch) $(Fmem) $(Fopt) $(Finline) $(Fipo) $(Fmath) $(Frpt) $(Fmkl)
#
## pgo gen flags [consult the documentation]
# F_FLAGS = $(FSTD) $(Farch) $(Fmem) $(Fopt) $(Finline) $(Fipo) $(Fprof) $(Frpt) $(Fpgog) -prof-dir/home/explict-path-of-build-directory/pgos 
## pgo use flags [consult the documentation]
# F_FLAGS = $(FSTD) $(Farch) $(Fmem) $(Fopt) $(Finline) $(Fipo) $(Fprof) $(Frpt) $(Fpgou) -prof-dir/home/explict-path-of-build-directory/pgos
#
# after adding mpi need to set following correctly 
 C_FLAGS    = 
#
# For link flags remove -static to avoid problems with mpi libs
# for some reason the link to svml is required with ifort when sse in enabled
#
# If math constraints include IEEE arithmetic flags then they must also appear here
  LINK_FLAGS =  $(Farch) $(Fmem) $(Fopt) $(Finline) $(Fipo) $(Fmath) -lsvml 

#> Set location of libraries/include files
#
 LIBRARIES = 

#> Set of *.o files

 OBJECTS =\
 	MAIN.o CLEB.o ACOEF.o SIN.o 

.SUFFIXES: .f 

quadge: $(OBJECTS)
	$(FC) $(LINK_FLAGS) $(OBJECTS) $(LIBRARIES) -o $@

.f.o:
	$(FC) $(SRCS) $(F_FLAGS) -c $<

clean:
	rm -f quadge *.o *.mod *.f90
