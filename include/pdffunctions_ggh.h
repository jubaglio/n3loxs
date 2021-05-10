/**
 *
 * \file    pdffunctions_ggh.h
 * \author  Julien Baglio
 * \date    February 2021
 *
 */

/**
 *
 * \brief   The header for all pdf luminosities called by the program.
 *
 */

#include "LHAPDF/LHAPDF.h"

// g-g luminosity
double dlumgg(const double, const double, const double, LHAPDF::PDF const* const);

// g-q luminosities
double dlumgq(const double, const double, const double, LHAPDF::PDF const* const);

// q-qbar luminosity
double dlumqqb(const double, const double, const double, LHAPDF::PDF const* const);

// q-q luminosity
double dlumqq(const double, const double, const double, LHAPDF::PDF const* const);


// q1-q2 luminosity
double dlumq1q2(const double, const double, const double, LHAPDF::PDF const* const);
