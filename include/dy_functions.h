/**
 *
 * \file    dy_functions.h
 * \author  Julien Baglio
 * \date    September 2020
 *
 */

/**
 *
 * \brief   The header for all neutral Drell-Yan coefficient functions
 *
 */

double intpow(const double& ,int);

double delta(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double qqb_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NLO terms
double qqb_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO qqbar terms
double qqb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO qqbar terms

double gq_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq NLO terms
double gq_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq NNLO terms
double gq_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms

double gg_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms
double gg_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms

double qq_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms
double qq_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms

double qQq_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ NNLO terms
double ud_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud NNLO terms
double qQq_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ N3LO terms
double ud_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud N3LO terms

double qQqb_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar NNLO terms
double udb_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar NNLO terms
double qQqb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar N3LO terms
double udb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar N3LO terms
