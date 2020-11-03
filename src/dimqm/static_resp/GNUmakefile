   OBJ = dim_grid_static.o\
         dim_eval_fnl_stat.o\
         dim_tabcd_stat.o\
         dimqm_indDipoles_stat.o\
         dimqm_pol_static.o

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
