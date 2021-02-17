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

## working dir
WORKINGDIR = $(shell pwd)

# gsl library
GSL_CONFIG = $(shell pwd)/gsl-2.6/gsldir/bin/gsl-config

vpath %.cpp $(WORKINGDIR)/src/dy_gamma $(WORKINGDIR)/src/dy_w $(WORKINGDIR)/src/wh

HEADERS = $(wildcard *.hpp $(WORKINGDIR)/include/*.hpp)
HEADERS += $(wildcard *.h $(WORKINGDIR)/include/*.h)

CC += -I$(WORKINGDIR)/include/

## create object files

./build/%.o : %.cpp $(HEADERS)
	$(CC) -c -o $@ $<

OBJ_PHOT = pdffunctions.o softvirtual.o regularterms.o regularterms_gq.o regularterms_gg.o \
	regularterms_qq.o regularterms_q1q2.o regularterms_q1q2b.o main_dy_n3lo.o

OBJ_W = pdffunctions_w.o softvirtual_w.o regularterms_w.o regularterms_gubar.o regularterms_gg_w.o \
	regularterms_gu.o regularterms_cubar.o regularterms_qqbar.o regularterms_qq_w.o \
	regularterms_qqprime.o regularterms_qbarqprimebar.o regularterms_ds.o regularterms_ubarcbar.o main_dy_w_n3lo.o

OBJ_WH = pdffunctions_w.o softvirtual_wh.o regularterms_wh.o regularterms_gubar_wh.o regularterms_gg_wh.o \
	regularterms_gu_wh.o regularterms_cubar_wh.o regularterms_qqbar_wh.o regularterms_qq_wh.o \
	regularterms_qqprime_wh.o regularterms_qbarqprimebar_wh.o regularterms_ds_wh.o regularterms_ubarcbar_wh.o main_wh_n3lo.o

all: run_dy_w run_dy_phot run_wh

run_dy_phot: $(addprefix ./build/, $(OBJ_PHOT))
	$(CC) $(patsubst %,./build/%,$(OBJ_PHOT)) $(CC2) -o $@
	mv run_dy_phot ./subprogs/

run_dy_w: $(addprefix ./build/, $(OBJ_W))
	$(CC) $(patsubst %,./build/%,$(OBJ_W)) $(CC2) -o $@
	mv run_dy_w ./subprogs/

run_wh: $(addprefix ./build/, $(OBJ_WH))
	$(CC) $(patsubst %,./build/%,$(OBJ_WH)) $(CC2) -o $@
	mv run_wh ./subprogs/

clean:
	rm -f ./build/*.o

distclean:
	rm -f ./build/*.o ./subprogs/*
