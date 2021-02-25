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

double delta_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2_bbh(const double*, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double bbb_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar NLO terms
double bbb_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar NNLO terms, low x1-region
double bbb_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar NNLO terms, intermediate x1-region
double bbb_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar NNLO terms, large x1-region
double bbb_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar N3LO terms, low x1-region
double bbb_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar N3LO terms, intermediate x1-region
double bbb_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-bbar N3LO terms, large x1-region

double bg_regular_nlo(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NLO terms
double bg_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NNLO terms, low x1-region
double bg_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NNLO terms, intermediate x1-region
double bg_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb NNLO terms, large x1-region
double bg_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb N3LO terms, low x1-region
double bg_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb N3LO terms, intermediate x1-region
double bg_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gb N3LO terms, large x1-region

double gg_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, low x1-region
double gg_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, intermediate x1-region
double gg_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg NNLO terms, large x1-region
double gg_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, low x1-region
double gg_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, intermediate x1-region
double gg_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gg N3LO terms, large x1-region

double qg_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms, low x1-region
double qg_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms, intermediate x1-region
double qg_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular gq N3LO terms, large x1-region

double bb_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) NNLO terms, low x1-region
double bb_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) NNLO terms, intermediate x1-region
double bb_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) NNLO terms, large x1-region
double bb_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) N3LO terms, low x1-region
double bb_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) N3LO terms, intermediate x1-region
double bb_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular b-b (and bbar-bbar) N3LO terms, large x1-region

double qqb_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms, low x1-region
double qqb_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms, intermediate x1-region
double qqb_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar NNLO terms, large x1-region
double qqb_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms, low x1-region
double qqb_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms, intermediate x1-region
double qqb_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular qqbar N3LO terms, large x1-region

double bq_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq NNLO terms, low x1-region
double bq_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq NNLO terms, intermediate x1-region
double bq_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq NNLO terms, large x1-region
double bq_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq N3LO terms, low x1-region
double bq_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq N3LO terms, intermediate x1-region
double bq_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bq N3LO terms, large x1-region

double bqb_regular_nnlo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar NNLO terms, low x1-region
double bqb_regular_nnlo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar NNLO terms, intermediate x1-region
double bqb_regular_nnlo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar NNLO terms, large x1-region
double bqb_regular_n3lo_z(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar N3LO terms, low x1-region
double bqb_regular_n3lo_w(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar N3LO terms, intermediate x1-region
double bqb_regular_n3lo_zb(const double*, const double, const double, LHAPDF::PDF const* const); // Regular bqbar N3LO terms, large x1-region
