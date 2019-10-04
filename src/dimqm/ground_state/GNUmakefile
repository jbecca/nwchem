   OBJ = dimqm_gs.o \
         dim_grid_gs.o\
         dim_eval_fnl_gs.o\
         dimqm_EqmE_gs.o\
         dimqm_EqmN_gs.o


   LIBRARY = libdimqm.a

   LIB_DEFINES = -DDBLEPREC_RYS

   LIB_INCLUDES = -I../../nwdft/include \
                  -I../../ddscf \
                  -I../../util \
                  -I../../nwdft/grid \
                  -I../../NWints/hondo

include ../../config/makefile.h
include ../../config/makelib.h

#erf.o : erf.f90
#	$(FC) -c -fpp erf.f90
#FDEBUG += -O0
