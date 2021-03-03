/**
 *
 * \file    bbh_function.h
 * \author  Julien Baglio
 * \date    February 2021
 *
 */

/**
 *
 * \brief   The header for all bbH coefficient functions
 *
 */

double intpow(const double& ,int);

double delta_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double bbb_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar NLO terms
double bbb_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar N3LO terms
double bbb_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar N3LO terms

double bg_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NLO terms
double bg_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NNLO terms
double bg_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NNLO terms

double gg_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms
double gg_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms

double qg_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms

double bb_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) NNLO terms
double bb_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) N3LO terms

double qqb_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms
double qqb_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms

double bq_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq NNLO terms
double bq_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq N3LO terms

double bqb_regular_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar NNLO terms
double bqb_regular_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar N3LO terms
