/**
 *
 * \file    dy_functions_w.h
 * \author  Julien Baglio
 * \date    October 2020
 *
 */

/**
 *
 * \brief   The header for all charged Drell-Yan coefficient functions
 *
 */

double delta(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double dub_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular dubar NLO terms

double dub_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO dubar terms
double dub_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO dubar terms

double gub_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar NLO terms
double gub_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar NNLO terms
double gub_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar N3LO terms

double gg_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms
double gg_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms

double gdb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gdbar N3LO terms

double cub_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar NNLO terms
double cub_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar N3LO terms

double qqb_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms
double qqb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms

double qq_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms
double qq_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms

double qqprime_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime NNLO terms
double qqprime_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime N3LO terms

double qbqprimeb_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbqprimeb NNLO terms
double qbqprimeb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbqprimeb N3LO terms

double ds_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms
double ds_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms

double ubcb_regular_nnlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ubcb NNLO terms
double ubcb_regular_n3lo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ubcb N3LO terms
