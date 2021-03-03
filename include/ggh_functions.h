/**
 *
 * \file    ggh_function.h
 * \author  Julien Baglio
 * \date    February 2021
 *
 */

/**
 *
 * \brief   The header for all ggH coefficient functions
 *
 */

double intpow(const double& ,int);

double delta_ggh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst_ggh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1_ggh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2_ggh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double gg_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular g-g NLO terms
double gg_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular g-g NNLO terms
double gg_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular g-g N3LO terms

double gq_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq NLO terms
double gq_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq NNLO terms
double gq_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms

double qqb_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q-qbar NLO terms
double qqb_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q-qbar NNLO terms
double qqb_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q-qbar N3LO terms

double qq_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q-q NNLO terms
double qq_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q-q N3LO terms

double q1q2_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q1-q2 NNLO terms
double q1q2_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular q1-q2 N3LO terms
