/**
 *
 * \file    zh_functions.h
 * \author  Julien Baglio
 * \date    Oct 2021
 *
 */

/**
 *
 * \brief   The header for all Z Higgs-strahlung coefficient functions
 *
 */

double delta_zh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst_zh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1_zh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2_zh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double qqb_regular_zh_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NLO terms
double qqb_regular_zh_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular NNLO qqbar terms
double qqb_regular_zh_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular N3LO qqbar terms

double gq_regular_zh_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq NLO terms
double gq_regular_zh_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq NNLO terms
double gq_regular_zh_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms

double gg_regular_zh_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms
double gg_regular_zh_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms

double qq_regular_zh_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms
double qq_regular_zh_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms

double qQq_regular_zh_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qQ + ud NNLO terms
double qQq_regular_zh_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qQ + ud N3LO terms

double qQqb_regular_zh_nnlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qQbar + udbar NNLO terms
double qQqb_regular_zh_n3lo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qQbar + udbar N3LO terms
