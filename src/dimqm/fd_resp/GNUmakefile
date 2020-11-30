   OBJ = dimqm_FD_cycle.o\
         dimqm_iterative_FD.o\
         dimqm_pol_FD.o\
         dim_grid_FD.o\
         dim_grid_FDa.o\
         dim_grid_FDb.o\
         dim_tabcd_FD.o\
         dim_eval_fnl_FD.o\
         scat_mat_cmplx.o\
         dimqm_EqmE_FD.o

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
