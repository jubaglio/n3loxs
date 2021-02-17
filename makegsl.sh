#!/bin/sh

########################################################################
#              -*- Compile gsl-2.6 for usage in dyn3lo -*-             #
########################################################################


WORKINGDIR=$(pwd)

cd ${WORKINGDIR}/gsl-2.6

./configure --prefix=${WORKINGDIR}/gsl-2.6/gsldir

make && make install
