/**
 *
 * \file    pdffunctions.h
 * \author  Julien Baglio
 * \date    September 2020
 *
 */

/**
 *
 * \brief   The header for all pdf luminosities called by the program.
 *
 */

#include "LHAPDF/LHAPDF.h"

// q-qbar luminosities
double dlumqqb(const double, const double, const double, LHAPDF::PDF const* const);
double dlumuub(const double, const double, const double, LHAPDF::PDF const* const);
double dlumddb(const double, const double, const double, LHAPDF::PDF const* const);

// g-q luminosities
double dlumgq(const double, const double, const double, LHAPDF::PDF const* const);
double dlumgu(const double, const double, const double, LHAPDF::PDF const* const);
double dlumgd(const double, const double, const double, LHAPDF::PDF const* const);

// g-g luminosities
double dlumgg(const double, const double, const double, LHAPDF::PDF const* const);

// q-q luminosities
double dlumqq(const double, const double, const double, LHAPDF::PDF const* const);

// non-singlet q-Q and u-d luminosities
double dlumqQq(const double, const double, const double, LHAPDF::PDF const* const);
double dlumud(const double, const double, const double, LHAPDF::PDF const* const);

// non-singlet q-Qbar and u-dbar luminosities
double dlumqQqb(const double, const double, const double, LHAPDF::PDF const* const);
double dlumudb(const double, const double, const double, LHAPDF::PDF const* const);
