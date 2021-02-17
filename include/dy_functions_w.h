/**
 *
 * \file    dy_functions.h
 * \author  Julien Baglio
 * \date    October 2020
 *
 */

/**
 *
 * \brief   The header for all Drell-Yann coefficient functions
 *
 */

double delta(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double dub_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular dubar NLO terms
double dub_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO dubar terms, low x1-region
double dub_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO dubar terms, intermediate x1-region
double dub_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO dubar terms, large x1-region
double dub_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO dubar terms, low x1-region
double dub_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO dubar terms, intermediate x1-region
double dub_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO dubar terms, large x1-region

double gub_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar NLO terms
double gub_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar NNLO terms, low x1-region
double gub_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar NNLO terms, intermediate x1-region
double gub_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar NNLO terms, large x1-region
double gub_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar N3LO terms, low x1-region
double gub_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar N3LO terms, intermediate x1-region
double gub_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar N3LO terms, large x1-region

double gg_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, low x1-region
double gg_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, intermediate x1-region
double gg_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, large x1-region
double gg_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, low x1-region
double gg_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, intermediate x1-region
double gg_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, large x1-region

double gdb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gdbar N3LO terms, low x1-region
double gdb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gdbar N3LO terms, intermediate x1-region
double gdb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gdbar N3LO terms, large x1-region

double cub_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar NNLO terms, low x1-region
double cub_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar NNLO terms, intermediate x1-region
double cub_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar NNLO terms, large x1-region
double cub_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar N3LO terms, low x1-region
double cub_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar N3LO terms, intermediate x1-region
double cub_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar N3LO terms, large x1-region

double qqb_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms, low x1-region
double qqb_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms, intermediate x1-region
double qqb_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms, large x1-region
double qqb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms, low x1-region
double qqb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms, intermediate x1-region
double qqb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms, large x1-region

double qq_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms, low x1-region
double qq_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms, intermediate x1-region
double qq_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms, large x1-region
double qq_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms, low x1-region
double qq_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms, intermediate x1-region
double qq_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms, large x1-region

double qqprime_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime NNLO terms, low x1-region
double qqprime_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime NNLO terms, intermediate x1-region
double qqprime_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime NNLO terms, large x1-region
double qqprime_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime N3LO terms, low x1-region
double qqprime_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime N3LO terms, intermediate x1-region
double qqprime_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime N3LO terms, large x1-region

double qbqprimeb_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar NNLO terms, low x1-region
double qbqprimeb_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar NNLO terms, intermediate x1-region
double qbqprimeb_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar NNLO terms, large x1-region
double qbqprimeb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar N3LO terms, low x1-region
double qbqprimeb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar N3LO terms, intermediate x1-region
double qbqprimeb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar N3LO terms, large x1-region

double ds_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms, low x1-region
double ds_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms, intermediate x1-region
double ds_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms, large x1-region
double ds_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms, low x1-region
double ds_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms, intermediate x1-region
double ds_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms, large x1-region

double ubcb_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms, low x1-region
double ubcb_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms, intermediate x1-region
double ubcb_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds NNLO terms, large x1-region
double ubcb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms, low x1-region
double ubcb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms, intermediate x1-region
double ubcb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds N3LO terms, large x1-region
