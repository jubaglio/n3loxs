/**
 *
 * \file    ncdy_kernels.h
 * \author  Julien Baglio
 * \date    Sep 2021
 *
 */

/**
 *
 * \brief   The header for all neutral Drell-Yan coefficient functions kernels (vector+axial parts)
 *
 */

#include <tuple>

double intpow(const double& ,int);


/**
 *    Vector part:
 */

std::tuple<double, double> delta_kernel(const double, const int);
double PlusConst_kernel(const double, const double, const int);
double PlusInt_kernel(const double, const double, const int);

double qqb_regular_kernel_nlo(const double, const double);
std::tuple<double, double> qqb_regular_kernel_nnlo(const double, const double);
std::tuple<double, double, double> qqb_regular_kernel_n3lo(const double, const double);

double gq_regular_kernel_nlo(const double, const double);
double gq_regular_kernel_nnlo(const double, const double);
std::tuple<double, double, double> gq_regular_kernel_n3lo(const double, const double);

double gg_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> gg_regular_kernel_n3lo(const double, const double);

double qq_regular_kernel_nnlo(const double, const double);
double qq_regular_kernel_n3lo(const double, const double);

double qQq_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> ud_regular_kernel_nnlo(const double, const double);
double qQq_regular_kernel_n3lo(const double, const double);
std::tuple<double, double> ud_regular_kernel_n3lo(const double, const double);

double qQqb_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> udb_regular_kernel_nnlo(const double, const double);
double qQqb_regular_kernel_n3lo(const double, const double);
std::tuple<double, double> udb_regular_kernel_n3lo(const double, const double);

/**
 *    Axial part:
 */

std::tuple<double, double> delta_axial_kernel(const double, const int);
std::tuple<double, double> PlusConst_axial_kernel(const double, const double, const int);
std::tuple<double, double> PlusInt_axial_kernel(const double, const double, const int);

double qqb_axial_regular_kernel_nlo(const double, const double);
std::tuple<double, double, double> qqb_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double, double, double> qqb_axial_regular_kernel_n3lo(const double, const double);

double gq_axial_regular_kernel_nlo(const double, const double);
std::tuple<double, double> gq_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double, double, double> gq_axial_regular_kernel_n3lo(const double, const double);

double gg_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> gg_axial_regular_kernel_n3lo(const double, const double);

double qq_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> qq_axial_regular_kernel_n3lo(const double, const double);

double qQq_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> ud_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> qQq_axial_regular_kernel_n3lo(const double, const double);
std::tuple<double, double, double> ud_axial_regular_kernel_n3lo(const double, const double);

double qQqb_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> udb_axial_regular_kernel_nnlo(const double, const double);
std::tuple<double, double> qQqb_axial_regular_kernel_n3lo(const double, const double);
std::tuple<double, double, double> udb_axial_regular_kernel_n3lo(const double, const double);
