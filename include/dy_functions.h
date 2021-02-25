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

double delta(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2(const double*, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double qqb_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NLO terms
double qqb_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO qqbar terms, low x1-region
double qqb_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO qqbar terms, intermediate x1-region
double qqb_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular NNLO qqbar terms, large x1-region
double qqb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO qqbar terms, low x1-region
double qqb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO qqbar terms, intermediate x1-region
double qqb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular N3LO qqbar terms, large x1-region


double gq_regular_nlo(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq NLO terms
double gq_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq NNLO terms, low x1-region
double gq_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq NNLO terms, intermediate x1-region
double gq_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq NNLO terms, large x1-region
double gq_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms, low x1-region
double gq_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms, intermediate x1-region
double gq_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms, large x1-region


double gg_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, low x1-region
double gg_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, intermediate x1-region
double gg_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, large x1-region
double gg_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, low x1-region
double gg_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, intermediate x1-region
double gg_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, large x1-region


double qq_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms, low x1-region
double qq_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms, intermediate x1-region
double qq_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq NNLO terms, large x1-region
double qq_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms, low x1-region
double qq_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms, intermediate x1-region
double qq_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq N3LO terms, large x1-region


double qQq_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ NNLO terms, low x1-region
double ud_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud NNLO terms, low x1-region
double qQq_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ NNLO terms, intermediate x1-region
double ud_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud NNLO terms, intermediate x1-region
double qQq_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ NNLO terms, large x1-region
double ud_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud NNLO terms, large x1-region
double qQq_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ N3LO terms, low x1-region
double ud_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud NNLO terms, low x1-region
double qQq_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ N3LO terms, intermediate x1-region
double ud_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud N3LO terms, intermediate x1-region
double qQq_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQ N3LO terms, large x1-region
double ud_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular ud N3LO terms, large x1-region


double qQqb_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar NNLO terms, low x1-region
double udb_regular_nnlo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar NNLO terms, low x1-region
double qQqb_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar NNLO terms, intermediate x1-region
double udb_regular_nnlo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar NNLO terms, intermediate x1-region
double qQqb_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar NNLO terms, large x1-region
double udb_regular_nnlo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar NNLO terms, large x1-region
double qQqb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar N3LO terms, low x1-region
double udb_regular_n3lo_z(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar NNLO terms, low x1-region
double qQqb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar N3LO terms, intermediate x1-region
double udb_regular_n3lo_w(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar N3LO terms, intermediate x1-region
double qQqb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular qQbar N3LO terms, large x1-region
double udb_regular_n3lo_zb(const double*, const double, const double, const double, LHAPDF::PDF const* const); // Regular udbar N3LO terms, large x1-region

