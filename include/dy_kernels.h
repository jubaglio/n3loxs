/**
 *
 * \file    dy_kernels.h
 * \author  Julien Baglio
 * \date    May 2021
 *
 */

/**
 *
 * \brief   The header for all Drell-Yan coefficient functions kernels
 *
 */

#include <tuple>

double intpow(const double& ,int);

double qqb_regular_kernel_nlo(const double, const double);
std::pair<double, double> qqb_regular_kernel_nnlo(const double, const double);
std::tuple<double, double, double> qqb_regular_kernel_n3lo(const double, const double);

double gq_regular_kernel_nlo(const double, const double);
double gq_regular_kernel_nnlo(const double, const double);
std::tuple<double, double, double> gq_regular_kernel_n3lo(const double, const double);

double gg_regular_kernel_nnlo(const double, const double);
std::pair<double, double> gg_regular_kernel_n3lo(const double, const double);

double qq_regular_kernel_nnlo(const double, const double);
double qq_regular_kernel_n3lo(const double, const double);

double qQq_regular_kernel_nnlo(const double, const double);
double ud_regular_kernel_nnlo(const double, const double);
double qQq_regular_kernel_n3lo(const double, const double);
double ud_regular_kernel_n3lo(const double, const double);

double qQqb_regular_kernel_nnlo(const double, const double);
double udb_regular_kernel_nnlo(const double, const double);
double qQqb_regular_kernel_n3lo(const double, const double);
double udb_regular_kernel_n3lo(const double, const double);
