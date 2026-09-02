# Toolchain selector. Defaults to GCC + OpenMPI; override via either
#   make COMPILER=intel ...
# or by setting COMPILER=intel in ../config.env (which run.sh exports).
COMPILER ?= gnu

ifeq ($(COMPILER),gnu)
  include makefile.inc.gnu
else ifeq ($(COMPILER),intel)
  include makefile.inc.intel
else
  $(error Unknown COMPILER='$(COMPILER)'. Supported values: gnu, intel.)
endif
