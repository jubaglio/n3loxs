/**
 *
 * \file    pdffunctions_bbh.h
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

// b-bbar luminosity
double dlumbbb(const double, const double, const double, LHAPDF::PDF const* const);

// g-q luminosities
double dlumgb(const double, const double, const double, LHAPDF::PDF const* const);
double dlumgq(const double, const double, const double, LHAPDF::PDF const* const);

// g-g luminosity
double dlumgg(const double, const double, const double, LHAPDF::PDF const* const);

// q-qbar luminosity
double dlumqqb(const double, const double, const double, LHAPDF::PDF const* const);

// b-q luminosities
double dlumbb(const double, const double, const double, LHAPDF::PDF const* const);
double dlumbq(const double, const double, const double, LHAPDF::PDF const* const);
double dlumbqb(const double, const double, const double, LHAPDF::PDF const* const);
