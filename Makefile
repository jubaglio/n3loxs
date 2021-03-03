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
#GSL_CONFIG = $(shell pwd)/gsl-2.6/gsldir/bin/gsl-config
GSL_CONFIG = gsl-config # test purposes, erase this line before commit

vpath %.cpp $(WORKINGDIR)/src/ $(WORKINGDIR)/src/dy_gamma $(WORKINGDIR)/src/dy_w
vpath %.cpp $(WORKINGDIR)/src/wh $(WORKINGDIR)/src/bbh $(WORKINGDIR)/src/ggh

HEADERS = $(wildcard *.hpp $(WORKINGDIR)/include/*.hpp)
HEADERS += $(wildcard *.h $(WORKINGDIR)/include/*.h)

CC += -I$(WORKINGDIR)/include/

## create object files

./build/%.o : %.cpp $(HEADERS)
	$(CC) -c -o $@ $<

OBJ_PHOT = alphaS.o pdffunctions.o softvirtual.o regularterms.o regularterms_gq.o regularterms_gg.o \
	regularterms_qq.o regularterms_q1q2.o regularterms_q1q2b.o main_dy_n3lo.o

OBJ_W = alphaS.o pdffunctions_w.o softvirtual_w.o regularterms_w_all.o regularterms_w.o regularterms_gubar.o regularterms_gg_w.o \
	regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_dy_w_n3lo.o

OBJ_WH = alphaS.o pdffunctions_w.o softvirtual_wh.o regularterms_wh_all.o regularterms_w.o regularterms_gubar.o regularterms_gg_w.o \
	regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_wh_n3lo.o

OBJ_BBH = alphaS.o pdffunctions_bbh.o softvirtual_bbh.o regularterms_bbh.o regularterms_bg_bbh.o regularterms_gg_bbh.o \
	regularterms_bq_bbh.o regularterms_bqbar_bbh.o regularterms_qqbar_bbh.o regularterms_bb_bbh.o \
	regularterms_qg_bbh.o main_bbh_n3lo.o

OBJ_GGH = alphaS.o pdffunctions_ggh.o exact_lo_ggh.o softvirtual_ggh.o regularterms_ggh.o regularterms_gq_ggh.o regularterms_qqbar_ggh.o \
	regularterms_qq_ggh.o regularterms_q1q2_ggh.o main_ggh_n3lo.o

all: run_dy_w run_dy_phot run_wh run_bbh run_ggh

run_dy_phot: $(addprefix ./build/, $(OBJ_PHOT))
	$(CC) $(patsubst %,./build/%,$(OBJ_PHOT)) $(CC2) -o $@
	mv run_dy_phot ./subprogs/

run_dy_w: $(addprefix ./build/, $(OBJ_W))
	$(CC) $(patsubst %,./build/%,$(OBJ_W)) $(CC2) -o $@
	mv run_dy_w ./subprogs/

run_wh: $(addprefix ./build/, $(OBJ_WH))
	$(CC) $(patsubst %,./build/%,$(OBJ_WH)) $(CC2) -o $@
	mv run_wh ./subprogs/

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
