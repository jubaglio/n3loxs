########################################################################
#                          -*- Makefile -*-                            #
########################################################################

# Compiler and flags C++:
CC = c++ -std=c++11 -pthread -O3
# only for debug:
#CC += -Wall -Wextra

########################################################################
## LHAPDF library for hadronic collision
## please adapt to your local installation
LHAPDFCONFIG = lhapdf-config
LIBSLHAPDF   = -Wl,-rpath,$(shell $(LHAPDFCONFIG) --libdir) $(shell $(LHAPDFCONFIG) --ldflags)

CC += $(shell $(GSL_CONFIG) --cflags) $(shell $(LHAPDFCONFIG) --cflags)
CC2 = -Wl,-rpath,$(shell $(GSL_CONFIG) --prefix)/lib $(shell $(GSL_CONFIG) --libs) $(LIBSLHAPDF)

## generate directory obj, if not yet existing
$(shell mkdir -p build)

## generate directory subprogs, if not yet existing
$(shell mkdir -p subprogs)

## working dir
WORKINGDIR = $(shell pwd)

# gsl library
GSL_CONFIG = $(shell pwd)/gsl-2.6/gsldir/bin/gsl-config

vpath %.cpp $(WORKINGDIR)/src/ $(WORKINGDIR)/src/dy_w $(WORKINGDIR)/src/ncdy
vpath %.cpp $(WORKINGDIR)/src/wh $(WORKINGDIR)/src/zh $(WORKINGDIR)/src/bbh $(WORKINGDIR)/src/ggh
vpath %.cpp $(WORKINGDIR)/src/dy_w_bins $(WORKINGDIR)/src/ncdy_bins $(WORKINGDIR)/src/wh_dynamical $(WORKINGDIR)/src/zh_dynamical

HEADERS = $(wildcard *.hpp $(WORKINGDIR)/include/*.hpp)
HEADERS += $(wildcard *.h $(WORKINGDIR)/include/*.h)

CC += -I$(WORKINGDIR)/include/

## create object files

./build/%.o : %.cpp $(HEADERS)
	$(CC) -c -o $@ $<

OBJ_DY_NC = alphaS.o pdffunctions.o softvirtual_ncdy.o softvirtual_ncdy_axial.o regularterms_qqb_ncdy.o \
	regularterms_qqb_ncdy_axial.o regularterms_gq_ncdy.o regularterms_gq_ncdy_axial.o regularterms_gg_ncdy.o \
	regularterms_gg_ncdy_axial.o regularterms_qq_ncdy.o regularterms_qq_ncdy_axial.o regularterms_q1q2_ncdy.o \
	regularterms_q1q2_ncdy_axial.o regularterms_q1q2b_ncdy.o regularterms_q1q2b_ncdy_axial.o \
	all_functions_ncdy.o born_normalization_ncdy.o main_ncdy_n3lo.o

OBJ_DY_NC_BINS = alphaS_simplified.o pdffunctions.o regularterms_qqb_ncdy.o regularterms_qqb_ncdy_axial.o \
	regularterms_gq_ncdy.o regularterms_gq_ncdy_axial.o regularterms_gg_ncdy.o regularterms_gg_ncdy_axial.o \
	regularterms_qq_ncdy.o regularterms_qq_ncdy_axial.o regularterms_q1q2_ncdy.o regularterms_q1q2_ncdy_axial.o \
	regularterms_q1q2b_ncdy.o regularterms_q1q2b_ncdy_axial.o all_functions_ncdy_bins.o born_normalization_ncdy.o main_ncdy_bins_n3lo.o

OBJ_W = alphaS.o pdffunctions_w.o softvirtual_w.o regularterms_w_all.o regularterms_w.o regularterms_gubar.o regularterms_gg_w.o \
	regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_dy_w_n3lo.o

OBJ_W_BINS = alphaS_simplified.o pdffunctions_w.o softvirtual_dy_w_bins.o regularterms_dy_w_bins_all.o regularterms_w.o \
	regularterms_gubar.o regularterms_gg_w.o regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_dy_w_bins_n3lo.o

OBJ_WH = alphaS.o pdffunctions_w.o all_functions_wh.o regularterms_w.o regularterms_gubar.o regularterms_gg_w.o \
	regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_wh_n3lo.o

OBJ_WH_DYN = alphaS_simplified.o pdffunctions_w.o all_functions_wh_dynamical.o regularterms_w.o \
	regularterms_gubar.o regularterms_gg_w.o regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_wh_dynamical_n3lo.o

OBJ_ZH = alphaS.o pdffunctions.o softvirtual_ncdy.o softvirtual_ncdy_axial.o regularterms_qqb_ncdy.o \
	regularterms_qqb_ncdy_axial.o regularterms_gq_ncdy.o regularterms_gq_ncdy_axial.o regularterms_gg_ncdy.o \
	regularterms_gg_ncdy_axial.o regularterms_qq_ncdy.o regularterms_qq_ncdy_axial.o regularterms_q1q2_ncdy.o \
	regularterms_q1q2_ncdy_axial.o regularterms_q1q2b_ncdy.o regularterms_q1q2b_ncdy_axial.o \
	all_functions_zh.o main_zh_n3lo.o

OBJ_ZH_DYN = alphaS_simplified.o pdffunctions.o regularterms_qqb_ncdy.o regularterms_qqb_ncdy_axial.o \
	regularterms_gq_ncdy.o  regularterms_gq_ncdy_axial.o regularterms_gg_ncdy.o regularterms_gg_ncdy_axial.o \
	regularterms_qq_ncdy.o regularterms_qq_ncdy_axial.o regularterms_q1q2_ncdy.o regularterms_q1q2_ncdy_axial.o \
	regularterms_q1q2b_ncdy.o regularterms_q1q2b_ncdy_axial.o all_functions_zh_dynamical.o main_zh_dynamical_n3lo.o

OBJ_BBH = alphaS.o pdffunctions_bbh.o softvirtual_bbh.o regularterms_bbh.o regularterms_bg_bbh.o regularterms_gg_bbh.o \
	regularterms_bq_bbh.o regularterms_bqbar_bbh.o regularterms_qqbar_bbh.o regularterms_bb_bbh.o \
	regularterms_qg_bbh.o main_bbh_n3lo.o

OBJ_GGH = alphaS.o pdffunctions_ggh.o exact_lo_ggh.o softvirtual_ggh.o regularterms_ggh.o regularterms_gq_ggh.o regularterms_qqbar_ggh.o \
	regularterms_qq_ggh.o regularterms_q1q2_ggh.o main_ggh_n3lo.o

all: run_dy_w run_dy_nc run_wh run_zh run_bbh run_ggh run_dy_w_bins run_dy_nc_bins run_wh_dyn run_zh_dyn

run_dy_nc: $(addprefix ./build/, $(OBJ_DY_NC))
	$(CC) $(patsubst %,./build/%,$(OBJ_DY_NC)) $(CC2) -o $@
	mv run_dy_nc ./subprogs/

run_dy_nc_bins: $(addprefix ./build/, $(OBJ_DY_NC_BINS))
	$(CC) $(patsubst %,./build/%,$(OBJ_DY_NC_BINS)) $(CC2) -o $@
	mv run_dy_nc_bins ./subprogs/

run_dy_w: $(addprefix ./build/, $(OBJ_W))
	$(CC) $(patsubst %,./build/%,$(OBJ_W)) $(CC2) -o $@
	mv run_dy_w ./subprogs/

run_dy_w_bins: $(addprefix ./build/, $(OBJ_W_BINS))
	$(CC) $(patsubst %,./build/%,$(OBJ_W_BINS)) $(CC2) -o $@
	mv run_dy_w_bins ./subprogs/

run_wh: $(addprefix ./build/, $(OBJ_WH))
	$(CC) $(patsubst %,./build/%,$(OBJ_WH)) $(CC2) -o $@
	mv run_wh ./subprogs/

run_wh_dyn: $(addprefix ./build/, $(OBJ_WH_DYN))
	$(CC) $(patsubst %,./build/%,$(OBJ_WH_DYN)) $(CC2) -o $@
	mv run_wh_dyn ./subprogs/

run_zh: $(addprefix ./build/, $(OBJ_ZH))
	$(CC) $(patsubst %,./build/%,$(OBJ_ZH)) $(CC2) -o $@
	mv run_zh ./subprogs/

run_zh_dyn: $(addprefix ./build/, $(OBJ_ZH_DYN))
	$(CC) $(patsubst %,./build/%,$(OBJ_ZH_DYN)) $(CC2) -o $@
	mv run_zh_dyn ./subprogs/

run_bbh: $(addprefix ./build/, $(OBJ_BBH))
	$(CC) $(patsubst %,./build/%,$(OBJ_BBH)) $(CC2) -o $@
	mv run_bbh ./subprogs/

run_ggh: $(addprefix ./build/, $(OBJ_GGH))
	$(CC) $(patsubst %,./build/%,$(OBJ_GGH)) $(CC2) -o $@
	mv run_ggh ./subprogs/

clean:
	rm -f ./build/*.o

distclean:
	rm -f ./build/*.o ./subprogs/*
