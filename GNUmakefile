# Top level directory
top = .

# Flags
DEBUG	= FALSE
DIM	= 3
COMP    = gcc
USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_HIP   = FALSE

# Source directory
src_dirs := $(top)/src

INCLUDE_LOCATIONS += $(src_dirs)
VPATH_LOCATIONS   += $(src_dirs)
src_packages := $(foreach dir, $(src_dirs), $(dir)/Make.package)
include $(src_packages)

# AMReX home directory
AMREX_HOME = ../../../amrex

# Other AMReX directories
amrex_dirs 	:= Base Boundary LinearSolvers/MLMG AmrCore Amr
amrex_packages	+= $(foreach dir, $(amrex_dirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
include $(amrex_packages)
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules