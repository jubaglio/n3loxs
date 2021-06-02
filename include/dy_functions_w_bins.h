/**
 *
 * \file    dy_functions_w_bins.h
 * \author  Julien Baglio
 * \date    May 2021
 *
 */

/**
 *
 * \brief   The header for all charged Drell-Yan coefficient functions, including binning
 *
 */

double delta(const double*, const double, const double, const double, const double,
	     const double, const double, const double, const double, const int, LHAPDF::PDF const* const); // delta-terms
double PlusConst(const double*, const double, const double, const double, const double,
		 const double, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt1(const double*, const double, const double, const double, const double,
		const double, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double PlusInt2(const double*, const double, const double, const double, const double,
		const double, const double, const double, const double, const int, LHAPDF::PDF const* const); // PlusDistributions terms
double dub_regular_nlo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular dubar full NLO terms
double dub_regular_nnlo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular full NNLO dubar terms
double dub_regular_n3lo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular full N3LO dubar terms

double gub_regular_nlo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar full NLO terms
double gub_regular_nnlo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar full NNLO terms
double gub_regular_n3lo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular gubar full N3LO terms

double gg_regular_nnlo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg full NNLO terms
double gg_regular_n3lo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular gg full N3LO terms

double gdb_regular_n3lo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular gdbar full N3LO terms

double cub_regular_nnlo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar full NNLO terms
double cub_regular_n3lo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular cubar full N3LO terms

double qqb_regular_nnlo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar full NNLO terms
double qqb_regular_n3lo(const double*, const double, const double, const double, const double,
			const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqbar full N3LO terms

double qq_regular_nnlo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq full NNLO terms
double qq_regular_n3lo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qq full N3LO terms

double qqprime_regular_nnlo(const double*, const double, const double, const double, const double,
			    const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime full NNLO terms
double qqprime_regular_n3lo(const double*, const double, const double, const double, const double,
			    const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qqprime full N3LO terms

double qbqprimeb_regular_nnlo(const double*, const double, const double, const double, const double,
			      const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar full NNLO terms
double qbqprimeb_regular_n3lo(const double*, const double, const double, const double, const double,
			      const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular qbarqprimebar full N3LO terms

double ds_regular_nnlo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds full NNLO terms
double ds_regular_n3lo(const double*, const double, const double, const double, const double,
		       const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular ds full N3LO terms

double ubcb_regular_nnlo(const double*, const double, const double, const double, const double,
			 const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular ubar-cbar full NNLO terms
double ubcb_regular_n3lo(const double*, const double, const double, const double, const double,
			 const double, const double, const double, const double, LHAPDF::PDF const* const); // Regular ubar-cbar full N3LO terms
