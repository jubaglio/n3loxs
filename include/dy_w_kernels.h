/**
 *
 * \file    dy_functions_w_kernels.h
 * \author  Julien Baglio
 * \date    October 2020
 *
 */

/**
 *
 * \brief   The header for all charged Drell-Yan coefficient functions kernels
 *
 */

#include <tuple>

double intpow(const double& ,int);

double dub_regular_kernel_nlo(const double, const double);
std::pair<double, double> dub_regular_kernel_nnlo(const double, const double);
std::pair<double, double> dub_regular_kernel_n3lo(const double, const double);

double gub_regular_kernel_nlo(const double, const double);
double gub_regular_kernel_nnlo(const double, const double);
std::pair<double, double> gub_regular_kernel_n3lo(const double, const double);

double gg_regular_kernel_nnlo(const double, const double);
double gg_regular_kernel_n3lo(const double, const double);

std::pair<double, double> gdb_regular_kernel_n3lo(const double, const double);

double cub_regular_kernel_nnlo(const double, const double);
double cub_regular_kernel_n3lo(const double, const double);

std::pair<double, double> qqb_regular_kernel_nnlo(const double, const double);
std::pair<double, double> qqb_regular_kernel_n3lo(const double, const double);

double qq_regular_kernel_nnlo(const double, const double);
double qq_regular_kernel_n3lo(const double, const double);

std::pair<double, double> qqprime_regular_kernel_nnlo(const double, const double);
std::pair<double, double> qqprime_regular_kernel_n3lo(const double, const double);

std::pair<double, double> qbqprimeb_regular_kernel_nnlo(const double, const double);
std::pair<double, double> qbqprimeb_regular_kernel_n3lo(const double, const double);

double ds_regular_kernel_nnlo(const double, const double);
double ds_regular_kernel_n3lo(const double, const double);

double ubcb_regular_kernel_nnlo(const double, const double);
double ubcb_regular_kernel_n3lo(const double, const double);
