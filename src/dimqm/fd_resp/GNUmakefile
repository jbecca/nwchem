   OBJ = dimqm_FD_cycle.o\
         dimqm_iterative_FD.o\
         dimqm_pol_FD.o

   LIBRARY = libdimqm.a

   LIB_DEFINES = -DDBLEPREC_RYS

   LIB_INCLUDES = -I../../nwdft/include \
                  -I../../ddscf \
                  -I../../include \
                  -I../../util \
                  -I../../nwdft/grid \
                  -I../../NWints/hondo

include ../../config/makefile.h
include ../../config/makelib.h

#erf.o : erf.f90
#	$(FC) -c -fpp erf.f90
#FDEBUG += -O0
