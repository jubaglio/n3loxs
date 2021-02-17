/**
 *
 * \file    pdffunctions.h
 * \author  Julien Baglio
 * \date    September 2020
 *
 */

/**
 *
 * \brief   The header for all pdf luminosities called by the program. Interfaced to the Luminosity class where the interface with LHAPDF is defined
 *
 */

#include "LHAPDF/LHAPDF.h"

// d-ubar luminosity
double dlumdub(const double, const double, const double, LHAPDF::PDF const* const);
double dlumdub2(const double, const double, const double, LHAPDF::PDF const* const);

// g-ubar luminosities
double dlumgub(const double, const double, const double, LHAPDF::PDF const* const);
double dlumgub2(const double, const double, const double, LHAPDF::PDF const* const);

// g-u luminosities
double dlumgu(const double, const double, const double, LHAPDF::PDF const* const);
double dlumgu2(const double, const double, const double, LHAPDF::PDF const* const);

// g-g luminosities
double dlumgg(const double, const double, const double, LHAPDF::PDF const* const);

// c-ubar luminosities
double dlumcub(const double, const double, const double, LHAPDF::PDF const* const);
double dlumcub2(const double, const double, const double, LHAPDF::PDF const* const);

// q-qbar luminosities
double dlumqqb(const double, const double, const double, LHAPDF::PDF const* const);
double dlumqqb2(const double, const double, const double, LHAPDF::PDF const* const);

// q-q luminosities
double dlumqq(const double, const double, const double, LHAPDF::PDF const* const);

// q-qprime luminosities
double dlumqqprime(const double, const double, const double, LHAPDF::PDF const* const);
double dlumqqprime2(const double, const double, const double, LHAPDF::PDF const* const);

// qbar-qprimebar luminosities
double dlumqbqprimeb(const double, const double, const double, LHAPDF::PDF const* const);
double dlumqbqprimeb2(const double, const double, const double, LHAPDF::PDF const* const);

// d-s and ubar-cbar luminosities
double dlumds(const double, const double, const double, LHAPDF::PDF const* const);
double dlumubcb(const double, const double, const double, LHAPDF::PDF const* const);
