#include <iostream>
#include <fstream>
//#include <cmath> // sin, cos, exp
#include "qmc.hpp"

#include <stdlib.h>
#include <cstring>

#include "LHAPDF/LHAPDF.h"

// Auxiliary functions for the integration
#include "dy_functions_wh.h"

// Global constants and QCD parameters
#include "constants.h"

// Header for the routines alphaS(muR) and mb_msbar(muR)
#include "alphaS.h"

struct {
  double s;
  double xmuf;
  const LHAPDF::PDF* pdf;
} global_param;

#include "pdfpar_w.h"
struct parampdf_struc parampdf;

const int scale0 = constants::MH + constants::MW;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Functors for the integration
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// q-qbar partonic channel
struct functor_delta_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = delta(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_delta;


struct functor_PlusConst_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = PlusConst(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_PlusConst;

struct functor_PlusInt1_t  {
  unsigned long long int number_of_integration_variables = 3;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand = PlusInt1(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_PlusInt1;

struct functor_PlusInt2_t  {
  unsigned long long int number_of_integration_variables = 3;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand = PlusInt2(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_PlusInt2;

struct functor_RegNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand = dub_regular_nlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegNLO;

struct functor_RegNNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      dub_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      dub_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      dub_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegNNLO;

struct functor_RegN3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      dub_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      dub_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      dub_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegN3LO;


// g-ubar partonic channel
struct functor_gubar_NLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand = gub_regular_nlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gubar_NLO;

struct functor_gubar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      gub_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      gub_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      gub_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gubar_NNLO;

struct functor_gubar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      gub_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      gub_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      gub_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gubar_N3LO;


// g-g partonic channel
struct functor_gg_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      gg_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      gg_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      gg_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gg_NNLO;

struct functor_gg_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      gg_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      gg_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      gg_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gg_N3LO;


// g-dbar partonic channel
struct functor_gdbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      gdb_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      gdb_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      gdb_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gdbar_N3LO;


// c-ubar partonic channel
struct functor_cubar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      cub_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      cub_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      cub_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_cubar_NNLO;

struct functor_cubar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      cub_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      cub_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      cub_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_cubar_N3LO;


// q-qbar partonic channel
struct functor_qqbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qqb_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      qqb_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      qqb_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqbar_NNLO;

struct functor_qqbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qqb_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      qqb_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      qqb_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqbar_N3LO;


// q-q partonic channel
struct functor_qq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qq_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      qq_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      qq_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qq_NNLO;

struct functor_qq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qq_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      qq_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      qq_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qq_N3LO;


// q-qprime partonic channel
struct functor_qqprime_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qqprime_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      qqprime_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      qqprime_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqprime_NNLO;

struct functor_qqprime_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qqprime_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      qqprime_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      qqprime_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqprime_N3LO;


// qbar-qprimebar partonic channel
struct functor_qbarqprimebar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qbqprimeb_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      qbqprimeb_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      qbqprimeb_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qbarqprimebar_NNLO;

struct functor_qbarqprimebar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      qbqprimeb_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      qbqprimeb_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      qbqprimeb_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qbarqprimebar_N3LO;


// ubar-cbar + d-s partonic channel
struct functor_ubarcbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      ubcb_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      ds_regular_nnlo_z(x, global_param.s, muf, global_param.pdf) +
      ubcb_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      ds_regular_nnlo_w(x, global_param.s, muf, global_param.pdf) +
      ubcb_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf) +
      ds_regular_nnlo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_ubarcbar_NNLO;

struct functor_ubarcbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
    muf = global_param.xmuf*scale0;

    integrand =
      ubcb_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      ds_regular_n3lo_z(x, global_param.s, muf, global_param.pdf) +
      ubcb_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      ds_regular_n3lo_w(x, global_param.s, muf, global_param.pdf) +
      ubcb_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf) +
      ds_regular_n3lo_zb(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_ubarcbar_N3LO;


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  if(argc==1)
    {
      printf("\nNo argument given, program will stop!!\n");
      exit(1);
    }
  if(argc==2 && (std::strcmp(argv[1],"--help") == 0 || std::strcmp(argv[1],"-h") == 0))
    {
      std::cout << "Usage:  " << argv[0] << " a b c d e f g h i (j) with:" << std::endl;
      std::cout << "a:  Lattice size (integer)" << std::endl;
      std::cout << "b:  Seed (integer)" << std::endl;
      std::cout << "c:  QCD order (integer, between O for LO and 3 for N3LO)" << std::endl;
      std::cout << "d:  p-p (0) or p-pbar (1) collider" << std::endl;
      std::cout << "e:  Hadronic energy in TeV (double)" << std::endl;
      std::cout << "f:  W+ (1) or W- (-1) production (integer)" << std::endl;
      std::cout << "g:  x_muf so that mu_F = x_muf*(MH+MW) (double)" << std::endl;
      std::cout << "h:  PDF set (string)" << std::endl;
      std::cout << "i:  PDF member (integer)" << std::endl;
      std::cout << "j:  --scale: optional flag to calculate various mu_R predictions. If absent, mu_R = Q" << std::endl;
      return 0;
    }
  if(argc < 10)
    {
      printf("\nNot enough arguments, program will stop!!\n");
      exit(1);
    }
  else
    {	
      int lattice  = atoi(argv[1]);
      int seed     = atoi(argv[2]);
      int qcdorder = atoi(argv[3]);

      const unsigned int MAXVAR = 3;

      int collider = atoi(argv[4]);
      if(collider==0)
	{
	  parampdf.collidertype = 1;
	}
      else
	{
	  parampdf.collidertype = -1;
	}

      double energy = atof(argv[5]); // energy in TeV
      double s;
      s = energy*energy*1.e6;
      int wchoice = atoi(argv[6]);
      wchoice = -wchoice; // internal conversion: all the contributions are written for W- by default.
      parampdf.wchoice = wchoice;
      double xmuf = atof(argv[7]);

      // init PDF set
      const std::string setname = argv[8];
      const int setimem = atoi(argv[9]);
      const LHAPDF::PDF* basepdf = LHAPDF::mkPDF( setname, setimem);
      LHAPDF::setVerbosity(0); // default is 1;
      
      // init parameters for all functors
      
      global_param.s       = s;
      global_param.xmuf    = xmuf;
      global_param.pdf     = basepdf;

      //integrators::Qmc<double,double,MAXVAR,integrators::transforms::Korobov<2>::type,integrators::fitfunctions::PolySingular::type> real_integrator;
      integrators::Qmc<double,double,MAXVAR,integrators::transforms::Korobov<2>::type> real_integrator;

      real_integrator.minn = lattice; // (optional) lattice size
      real_integrator.maxeval = 1;
      //      real_integrator.maxeval = 4;
      
      //      real_integrator.verbosity = 3; // default is zero
      real_integrator.minm = 12; // number of random shifts of the lattice, default is 32

      std::cout.precision(18);

      real_integrator.devices = {-1}; // only CPUs
      real_integrator.cputhreads = 1; // only one CPU
      real_integrator.randomgenerator.seed(seed);

      // define variables to store the results of QMC integrations

      integrators::result<double> resultdelta;
      double dubar_lo_result,dubar_lo_error;

      integrators::result<double> resultPlusConst;
      integrators::result<double> resultPlusInt1;
      integrators::result<double> resultPlusInt2;
      integrators::result<double> resultRegNLO;
      integrators::result<double> result_gubar_NLO;
      double dubar_nlo_result,dubar_nlo_error;

      integrators::result<double> resultRegNNLO;
      integrators::result<double> result_gubar_NNLO;
      integrators::result<double> result_gg_NNLO;
      integrators::result<double> result_cubar_NNLO;
      integrators::result<double> result_qqbar_NNLO;
      integrators::result<double> result_qq_NNLO;
      integrators::result<double> result_qqprime_NNLO;
      integrators::result<double> result_qbarqprimebar_NNLO;
      integrators::result<double> result_ubarcbar_NNLO;
      double dubar_nnlo_result,dubar_nnlo_error;

      integrators::result<double> resultRegN3LO;
      integrators::result<double> result_gubar_N3LO;
      integrators::result<double> result_gg_N3LO;
      integrators::result<double> result_gdbar_N3LO;
      integrators::result<double> result_cubar_N3LO;
      integrators::result<double> result_qqbar_N3LO;
      integrators::result<double> result_qq_N3LO;
      integrators::result<double> result_qqprime_N3LO;
      integrators::result<double> result_qbarqprimebar_N3LO;
      integrators::result<double> result_ubarcbar_N3LO;
      double dubar_n3lo_result,dubar_n3lo_error;


      // Perform all the required QMC integrations

      if(qcdorder>=0)
	{
	  real_integrator.minn = lattice/10;
	  functor_delta.k = 0;
	  resultdelta = real_integrator.integrate(functor_delta);
	  dubar_lo_result = resultdelta.integral;
	  dubar_lo_error  = resultdelta.error;

	  if(qcdorder>=1)
	    {
	      functor_delta.k     = 1;
	      functor_PlusConst.k = 1;
	      functor_PlusInt1.k  = 1;
	      functor_PlusInt2.k  = 1;
	      resultdelta = real_integrator.integrate(functor_delta);
	      resultPlusConst = real_integrator.integrate(functor_PlusConst);

	      real_integrator.minn = lattice;
	      resultPlusInt1 = real_integrator.integrate(functor_PlusInt1);
	      resultPlusInt2 = real_integrator.integrate(functor_PlusInt2);
	      resultRegNLO = real_integrator.integrate(functor_RegNLO);
	      dubar_nlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNLO.integral;
	      dubar_nlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
	      			    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNLO.error,2));
	      result_gubar_NLO = real_integrator.integrate(functor_gubar_NLO);

	      if(qcdorder>=2)
		{
		  functor_delta.k     = 2;
		  functor_PlusConst.k = 2;
		  functor_PlusInt1.k  = 2;
		  functor_PlusInt2.k  = 2;
		  real_integrator.minn = lattice/10;
		  resultdelta = real_integrator.integrate(functor_delta);
		  resultPlusConst = real_integrator.integrate(functor_PlusConst);
		  real_integrator.minn = lattice;
		  resultPlusInt1 = real_integrator.integrate(functor_PlusInt1);
		  resultPlusInt2 = real_integrator.integrate(functor_PlusInt2);
		  resultRegNNLO    = real_integrator.integrate(functor_RegNNLO);
		  dubar_nnlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNNLO.integral;
		  dubar_nnlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		  			 pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNNLO.error,2));
		  result_gubar_NNLO   = real_integrator.integrate(functor_gubar_NNLO);
		  result_gg_NNLO   = real_integrator.integrate(functor_gg_NNLO);
		  result_cubar_NNLO   = real_integrator.integrate(functor_cubar_NNLO);
		  result_qqbar_NNLO   = real_integrator.integrate(functor_qqbar_NNLO);
		  result_qq_NNLO   = real_integrator.integrate(functor_qq_NNLO);
		  result_qqprime_NNLO   = real_integrator.integrate(functor_qqprime_NNLO);
		  result_qbarqprimebar_NNLO   = real_integrator.integrate(functor_qbarqprimebar_NNLO);
		  result_ubarcbar_NNLO   = real_integrator.integrate(functor_ubarcbar_NNLO);

		  if(qcdorder==3)
		    {
		      functor_delta.k     = 3;
		      functor_PlusConst.k = 3;
		      functor_PlusInt1.k  = 3;
		      functor_PlusInt2.k  = 3;
		      real_integrator.minn = lattice/10;
		      resultdelta = real_integrator.integrate(functor_delta);
		      resultPlusConst = real_integrator.integrate(functor_PlusConst);
		      real_integrator.minn = lattice;
		      real_integrator.minn = 10*lattice;
		      resultPlusInt1 = real_integrator.integrate(functor_PlusInt1);
		      real_integrator.minn = lattice;
		      resultPlusInt2 = real_integrator.integrate(functor_PlusInt2);
		      real_integrator.minn = 10*lattice;
		      resultRegN3LO    = real_integrator.integrate(functor_RegN3LO);
		      real_integrator.minn = 10*lattice;
		      dubar_n3lo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegN3LO.integral;
		      dubar_n3lo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		      			    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegN3LO.error,2));
		      result_gubar_N3LO   = real_integrator.integrate(functor_gubar_N3LO);
		      real_integrator.minn = lattice;
		      result_gg_N3LO   = real_integrator.integrate(functor_gg_N3LO);
		      result_cubar_N3LO   = real_integrator.integrate(functor_cubar_N3LO);
		      result_qqbar_N3LO   = real_integrator.integrate(functor_qqbar_N3LO);
		      result_qq_N3LO   = real_integrator.integrate(functor_qq_N3LO);
		      result_qqprime_N3LO   = real_integrator.integrate(functor_qqprime_N3LO);
		      result_qbarqprimebar_N3LO   = real_integrator.integrate(functor_qbarqprimebar_N3LO);
		      result_ubarcbar_N3LO   = real_integrator.integrate(functor_ubarcbar_N3LO);
		      result_gdbar_N3LO   = real_integrator.integrate(functor_gdbar_N3LO);
		    }
		}
	    }
	}

      // Building the cross section: adding back alphaS, evolve to muR, add DY normalization
      // Result of this code: xs(p p / p pbar -> W+/W- + H + X) in pb

      double BornDY = constants::gevtopb*constants::MW*constants::MW*constants::MW*constants::MW
      	/(48*constants::Pi*constants::Nc*constants::vev*constants::vev*constants::vev*constants::vev);

      double muf  = xmuf*scale0;
      double muf2 = muf*muf;
      double xmur;
      double mur;
      double mur2;

      double asopi;
      double asopi2;
      double asopi3;
      double asopi4;
      double asopi6;

      double logmu1;
      double logmu2;

      double xslo_result;
      double xslo_error;
      double xsnlo_result;
      double xsnlo_error;
      double xsnnlo_result;
      double xsnnlo_error;
      double xsn3lo_result;
      double xsn3lo_error;

      std::string finalfile;
      std::stringstream filename;
      std::string header;

      int imax;
      double dxmur;

      double asopimz = (basepdf->alphasQ(constants::MZ))/constants::Pi;

      if(argc>=11 && std::strcmp(argv[10],"--scale") == 0)
	{
	  imax = 16;
	  dxmur = 1.5/(imax-1);
	}
      else
	{
	  imax = 1;
	}

      if(collider==0)
	{
	  if(wchoice==1)
	    {
	      filename << "WminusH_xs_pp_" << energy << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p p -> W- H)";
	    }
	  else
	    {
	      filename << "WplusH_xs_pp_" << energy << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p p -> W+ H)";
	    }
	}
      else
	{
	  if(wchoice==1)
	    {
	      filename << "WminusH_xs_ppbar_" << energy << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p pbar -> W- H)";
	    }
	  else
	    {
	      filename << "WplusH_xs_ppbar_" << energy << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p pbar -> W+ H)";
	    }
	}
      
      
      filename >> finalfile;
      std::ofstream fa(finalfile);

      switch(qcdorder)
	{
	case 0:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central scale mu0 = MH+MW\n# mu_R/mu0\t" << "mu_F/mu0\t" << "xs_LO (pb)" << std::endl;
	  break;
	case 1:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central scale mu0 = MH+MW\n# mu_R/mu0\t" << "mu_F/mu0\t" << "xs_LO (pb)\t" << "xs_NLO (pb)" << std::endl;
	  break;
	case 2:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central scale mu0 = MH+MW\n# mu_R/mu0\t" << "mu_F/mu0\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t"
	     << "xs_NNLO (pb)" << std::endl;
	  break;
	case 3:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central scale mu0 = MH+MW\n# mu_R/mu0\t" << "mu_F/mu0\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t"
	     << "xs_NNLO (pb)\t" << "xs_N3LO (pb)" << std::endl;
	  break;
	}
      
      for(int i = 0; i<imax; i++)
	{
	  if(imax==1)
	    {
	      xmur = xmuf;
	    }
	  else
	    {
	      xmur = 0.5 + i*dxmur;
	    }
	  mur  = xmur*scale0;
	  mur2 = mur*mur;

      if(qcdorder>=0)
      	{
      	  xslo_result = BornDY*dubar_lo_result;
      	  xslo_error  = BornDY*dubar_lo_error;
	  if(qcdorder==0)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/scale0) << "\t" << (muf/scale0) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << std::scientific << xslo_error << std::endl;
	    }
      	}
      if(qcdorder>=1)
      	{
	  // alphaS(mur) at NLO
	  asopi = as_n3loxs(mur, 1, asopimz);
      	  asopi2 = asopi*asopi;

      	  xsnlo_result = BornDY*(dubar_lo_result + asopi*(dubar_nlo_result + result_gubar_NLO.integral));
      	  xsnlo_error  = BornDY*sqrt(pow(dubar_lo_error,2) + asopi2*(pow(dubar_nlo_error,2) + pow(result_gubar_NLO.error,2)));

	  if(qcdorder==1)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/scale0) << "\t" << (muf/scale0) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << std::endl;
	    }
      	}
      if(qcdorder>=2)
      	{
	  // alphaS(mur) at NNLO
	  asopi = as_n3loxs(mur, 2, asopimz);
	  asopi2 = asopi*asopi;
      	  asopi4 = asopi2*asopi2;
      	  logmu1 = log(mur2/muf2);

      	  xsnnlo_result = BornDY*
      	    (dubar_lo_result +
      	     asopi*(dubar_nlo_result + result_gubar_NLO.integral) +
      	     asopi2*(dubar_nnlo_result + result_gubar_NNLO.integral + result_gg_NNLO.integral +
      		     result_cubar_NNLO.integral + result_qqbar_NNLO.integral + result_qq_NNLO.integral + 
      		     result_qqprime_NNLO.integral + result_qbarqprimebar_NNLO.integral + result_ubarcbar_NNLO.integral + 
		     constants::b0*logmu1*(dubar_nlo_result + result_gubar_NLO.integral))
      	     );
      	  xsnnlo_error  = BornDY*
      	    sqrt(pow(dubar_lo_error,2) +
      		 asopi2*(pow(dubar_nlo_error,2) + pow(result_gubar_NLO.error,2)) +
      		 asopi4*(pow(dubar_nnlo_error,2) + pow(result_gubar_NNLO.error,2) +
      			 pow(result_gg_NNLO.error,2) + pow(result_cubar_NNLO.error,2) +
      			 pow(result_qqbar_NNLO.error,2) + pow(result_qq_NNLO.error,2) +
			 pow(result_qqprime_NNLO.error,2) + pow(result_qbarqprimebar_NNLO.error,2) + pow(result_ubarcbar_NNLO.error,2) + 
      			 pow(constants::b0*logmu1,2)*(pow(dubar_nlo_error,2) + pow(result_gubar_NLO.error,2)))
      		 );

	  if(qcdorder==2)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/scale0) << "\t" << (muf/scale0) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << std::endl;
	    }
      	}
      if(qcdorder==3)
      	{
	  // alphaS(mur) at N3LO
	  asopi = as_n3loxs(mur, 3, asopimz);
	  asopi2 = asopi*asopi;
	  asopi3 = asopi*asopi2;
      	  asopi4 = asopi2*asopi2;
      	  asopi6 = asopi3*asopi3;

      	  logmu2 = logmu1*logmu1;

      	  xsn3lo_result = BornDY*
      	    (dubar_lo_result +
      	     asopi*(dubar_nlo_result + result_gubar_NLO.integral) +
      	     asopi2*(dubar_nnlo_result + result_gubar_NNLO.integral + result_gg_NNLO.integral +
      		     result_cubar_NNLO.integral + result_qqbar_NNLO.integral + result_qq_NNLO.integral +
		     result_qqprime_NNLO.integral + result_qbarqprimebar_NNLO.integral + result_ubarcbar_NNLO.integral + 
      		     constants::b0*logmu1*(dubar_nlo_result + result_gubar_NLO.integral)) +
      	     asopi3*(dubar_n3lo_result + result_gubar_N3LO.integral + result_gg_N3LO.integral +
      		     result_cubar_N3LO.integral + result_qqbar_N3LO.integral + result_qq_N3LO.integral +
		     result_qqprime_N3LO.integral + result_qbarqprimebar_N3LO.integral +
		     result_ubarcbar_N3LO.integral + result_gdbar_N3LO.integral +
      		     (pow(constants::b0,2)*logmu2 + constants::b1*logmu1)*(dubar_nlo_result + result_gubar_NLO.integral) +
      		     2*constants::b0*logmu1*
      		     (dubar_nnlo_result + result_gubar_NNLO.integral + result_gg_NNLO.integral +
      		      result_cubar_NNLO.integral + result_qqbar_NNLO.integral +
		      result_qq_NNLO.integral + result_qqprime_NNLO.integral +
		      result_qbarqprimebar_NNLO.integral + result_ubarcbar_NNLO.integral))
      	     );
      	  xsn3lo_error = BornDY*
      	    sqrt(pow(dubar_lo_error,2) +
      		 asopi2*(pow(dubar_nlo_error,2) + pow(result_gubar_NLO.error,2)) +
      		 asopi4*(pow(dubar_nnlo_error,2) + pow(result_gubar_NNLO.error,2) + pow(result_gg_NNLO.error,2) +
      			 pow(result_cubar_NNLO.error,2) + pow(result_qqbar_NNLO.error,2) + pow(result_qq_NNLO.error,2) +
			 pow(result_qqprime_NNLO.error,2) + pow(result_qbarqprimebar_NNLO.error,2) + pow(result_ubarcbar_NNLO.error,2) + 
      			 pow(constants::b0*logmu1,2)*(pow(dubar_nlo_error,2) + pow(result_gubar_NLO.error,2))) +
      		 asopi6*(pow(dubar_n3lo_error,2) + pow(result_gubar_N3LO.error,2) + pow(result_gg_N3LO.error,2) +
      			 pow(result_cubar_N3LO.error,2) + pow(result_qqbar_N3LO.error,2) + pow(result_qq_N3LO.error,2) +
			 pow(result_qqprime_N3LO.error,2) + pow(result_qbarqprimebar_N3LO.error,2) +
			 pow(result_ubarcbar_N3LO.error,2) + pow(result_gdbar_N3LO.error,2) +
      			 pow((pow(constants::b0,2)*logmu2 + constants::b1*logmu1),2)*(pow(dubar_nlo_error,2) + pow(result_gubar_NLO.error,2)) +
      			 pow(2*constants::b0*logmu1,2)*
      			 (pow(dubar_nnlo_error,2) + pow(result_gubar_NNLO.error,2) + pow(result_gg_NNLO.error,2) +
      			  pow(result_cubar_NNLO.error,2) + pow(result_qqbar_NNLO.error,2) +
			  pow(result_qq_NNLO.error,2) + pow(result_qqprime_NNLO.error,2) +
			  pow(result_qbarqprimebar_NNLO.error,2) + pow(result_ubarcbar_NNLO.error,2)))
		 );

	  fa << std::fixed << std::setprecision(3) << (mur/scale0) << "\t" << (muf/scale0) << "\t"
	     << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result << "\t" << xsn3lo_result
	     << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << "\t" << xsn3lo_error << std::endl;
      	}
	}
      delete basepdf;
      fa.close();
      return 0;
    }
}
