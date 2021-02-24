#include <iostream>
#include <fstream>
//#include <cmath>
#include "qmc.hpp"

#include <stdlib.h>
#include <cstring>

#include "LHAPDF/LHAPDF.h"

// Auxiliary functions for the integration
#include "dy_functions.h"

// Global constants and QCD parameters
#include "constants.h"

// Header for the routines alphaS(muR) and mb_msbar(muR)
#include "alphaS.h"

struct {
  double s;
  double Q;
  double xmuf;
  const LHAPDF::PDF* pdf;
} global_param_gamma;

#include "pdfpar.h"
struct parampdf_struc parampdf;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Functors for the integration
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// q-qbar partonic channel
struct functor_delta_t  {
  unsigned long long int number_of_integration_variables = 1;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand = delta(x, global_param_gamma.s, Q2, muf, k, global_param_gamma.pdf);

    return integrand;

  }
} functor_delta;


struct functor_PlusConst_t  {
  unsigned long long int number_of_integration_variables = 1;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand = PlusConst(x, global_param_gamma.s, Q2, muf, k, global_param_gamma.pdf);

    return integrand;

  }
} functor_PlusConst;

struct functor_PlusInt1_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand = PlusInt1(x, global_param_gamma.s, Q2, muf, k, global_param_gamma.pdf);

    return integrand;

  }
} functor_PlusInt1;

struct functor_PlusInt2_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand = PlusInt2(x, global_param_gamma.s, Q2, muf, k, global_param_gamma.pdf);

    return integrand;

  }
} functor_PlusInt2;

struct functor_RegNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand = qqb_regular_nlo(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_RegNLO;

struct functor_RegNNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qqb_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qqb_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qqb_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_RegNNLO;

struct functor_RegN3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qqb_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qqb_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qqb_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_RegN3LO;


// g-q partonic channel
struct functor_gq_NLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand = gq_regular_nlo(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_gq_NLO;

struct functor_gq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      gq_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gq_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gq_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_gq_NNLO;

struct functor_gq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      gq_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gq_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gq_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_gq_N3LO;


// g-g partonic channel
struct functor_gg_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      gg_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gg_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gg_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_gg_NNLO;

struct functor_gg_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      gg_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gg_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      gg_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_gg_N3LO;


// q-q partonic channel
struct functor_qq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qq_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qq_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qq_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_qq_NNLO;

struct functor_qq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qq_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qq_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qq_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_qq_N3LO;


// q-Q partonic channel
struct functor_qQq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qQq_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      ud_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQq_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      ud_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQq_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      ud_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_qQq_NNLO;

struct functor_qQq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qQq_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      ud_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQq_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      ud_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQq_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      ud_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_qQq_N3LO;


// q-Qbar partonic channel
struct functor_qQqb_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qQqb_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      udb_regular_nnlo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQqb_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      udb_regular_nnlo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQqb_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      udb_regular_nnlo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_qQqb_NNLO;

struct functor_qQqb_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double Q2;
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    Q2 = global_param_gamma.Q*global_param_gamma.Q;
    muf = global_param_gamma.xmuf*global_param_gamma.Q;

    integrand =
      qQqb_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      udb_regular_n3lo_z(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQqb_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      udb_regular_n3lo_w(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      qQqb_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf) +
      udb_regular_n3lo_zb(x, global_param_gamma.s, Q2, muf, global_param_gamma.pdf);

    return integrand;

  }
} functor_qQqb_N3LO;


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
      std::cout << "f:  Invariant lepton-mass Q in GeV (double)" << std::endl;
      std::cout << "g:  x_muf so that mu_F = x_muf*Q (double)" << std::endl;
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

      const unsigned int MAXVAR = 2;

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
      double Q    = atof(argv[6]);
      double xmuf = atof(argv[7]);

      // init PDF set
      const std::string setname = argv[8];
      const int setimem = atoi(argv[9]);
      const LHAPDF::PDF* basepdf = LHAPDF::mkPDF( setname, setimem);
      LHAPDF::setVerbosity(0); // default is 1;

      // init parameters for all functors      
      global_param_gamma.s       = s;
      global_param_gamma.Q       = Q;
      global_param_gamma.xmuf    = xmuf;
      global_param_gamma.pdf     = basepdf;

      //      integrators::Qmc<double,double,MAXVAR,integrators::transforms::Korobov<2>::type,integrators::fitfunctions::PolySingular::type> real_integrator;
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
      double qqb_lo_result,qqb_lo_error;

      integrators::result<double> resultPlusConst;
      integrators::result<double> resultPlusInt1;
      integrators::result<double> resultPlusInt2;
      integrators::result<double> resultRegNLO;
      integrators::result<double> result_gq_NLO;
      double qqb_nlo_result,qqb_nlo_error;

      integrators::result<double> resultRegNNLO;
      integrators::result<double> result_gq_NNLO;
      integrators::result<double> result_gg_NNLO;
      integrators::result<double> result_qq_NNLO;
      integrators::result<double> result_qQq_NNLO;
      integrators::result<double> result_qQqb_NNLO;
      double qqb_nnlo_result,qqb_nnlo_error;

      integrators::result<double> resultRegN3LO;
      integrators::result<double> result_gq_N3LO;
      integrators::result<double> result_gg_N3LO;
      integrators::result<double> result_qq_N3LO;
      integrators::result<double> result_qQq_N3LO;
      integrators::result<double> result_qQqb_N3LO;
      double qqb_n3lo_result,qqb_n3lo_error;


      // Perform all the required QMC integrations

      if(qcdorder>=0)
	{
	  real_integrator.minn = lattice/10;
	  functor_delta.k = 0;
	  resultdelta = real_integrator.integrate(functor_delta);
	  qqb_lo_result = resultdelta.integral;
	  qqb_lo_error  = resultdelta.error;

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
	      qqb_nlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNLO.integral;
	      qqb_nlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
				    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNLO.error,2));
	      result_gq_NLO = real_integrator.integrate(functor_gq_NLO);

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
		  qqb_nnlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNNLO.integral;
		  qqb_nnlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
					pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNNLO.error,2));
		  result_gq_NNLO   = real_integrator.integrate(functor_gq_NNLO);
		  result_gg_NNLO   = real_integrator.integrate(functor_gg_NNLO);
		  result_qq_NNLO   = real_integrator.integrate(functor_qq_NNLO);
		  result_qQq_NNLO  = real_integrator.integrate(functor_qQq_NNLO);
		  result_qQqb_NNLO = real_integrator.integrate(functor_qQqb_NNLO);

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
		      qqb_n3lo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegN3LO.integral;
		      qqb_n3lo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
					    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegN3LO.error,2));
		      result_gq_N3LO   = real_integrator.integrate(functor_gq_N3LO);
		      real_integrator.minn = lattice;
		      result_gg_N3LO   = real_integrator.integrate(functor_gg_N3LO);
		      result_qq_N3LO   = real_integrator.integrate(functor_qq_N3LO);
		      result_qQq_N3LO  = real_integrator.integrate(functor_qQq_N3LO);
		      result_qQqb_N3LO = real_integrator.integrate(functor_qQqb_N3LO);
		    }
		}
	    }
	}

      // Building the cross section: adding back alphaS, evolve to muR, add DY normalization
      // Result of this code: xs(p p / p pbar -> gamma* + X) in pb for an off-shell photon
      // To obtain the differential cross section Q^2*dxs/dQ^2 including gamma* -> l+ l-: [result of this code] * e^2/(12*Pi^2)

      double BornDY;
      BornDY = constants::gevtopb*constants::Pi*constants::ee2/(constants::Nc*Q*Q);
      
      double muf  = xmuf*Q;
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
	  filename << "dy_xs_gamma_pp_" << energy << "tev_q" << int(Q) << "_pdf" << setimem << "_muf" << xmuf << ".txt";
	  header = "# Drell-Yan cross section xs(p p -> gamma*) at a given photon virtuality Q = ";
	}
      else
	{
	  filename << "dy_xs_gamma_ppbar_" << energy << "tev_q" << int(Q) << "_pdf" << setimem << "_muf" << xmuf << ".txt";
	  header = "# Drell-Yan cross section xs(p pbar -> gamma*) at a given photon virtuality Q = ";
	}

      filename >> finalfile;
      std::ofstream fa(finalfile);

      switch(qcdorder)
	{
	case 0:
	  fa << header << std::fixed << std::setprecision(2) << Q << " GeV, sqrt(S) = " << energy
	     << " TeV\n# mu_R/Q\t" << "mu_F/Q\t" << "xs_LO (pb)\t" << "num error (respective xs)" << std::endl;
	  break;
	case 1:
	  fa << header << std::fixed << std::setprecision(2) << Q << " GeV, sqrt(S) = " << energy
	     << " TeV\n# mu_R/Q\t" << "mu_F/Q\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t" << "num error (respective xs)" << std::endl;
	  break;
	case 2:
	  fa << header << std::fixed << std::setprecision(2) << Q << " GeV, sqrt(S) = " << energy
	     << " TeV\n# mu_R/Q\t" << "mu_F/Q\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t"
	     << "dxs_NNLO/dQ (pb/GeV)\t" << "num error (respective xs)" << std::endl;
	  break;
	case 3:
	  fa << header << std::fixed << std::setprecision(2) << Q << " GeV, sqrt(S) = " << energy
	     << " TeV\n# mu_R/Q\t" << "mu_F/Q\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t"
	     << "xs_NNLO (pb)\t" << "xs_N3LO (pb)\t" << "num error (respective xs)" << std::endl;
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
	  mur  = xmur*Q;
	  mur2 = mur*mur;

      if(qcdorder>=0)
      	{
      	  xslo_result = BornDY*qqb_lo_result;
      	  xslo_error  = BornDY*qqb_lo_error;
	  if(qcdorder==0)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/Q) << "\t" << (muf/Q) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << std::scientific << xslo_error << std::endl;
	    }
      	}
      if(qcdorder>=1)
      	{
	  // alphaS(mur) and a mb(mur) at NLO
	  asopi   = as_n3loxs(mur, 1, asopimz);
      	  asopi2 = asopi*asopi;

      	  xsnlo_result = BornDY*(qqb_lo_result + asopi*(qqb_nlo_result + result_gq_NLO.integral));
      	  xsnlo_error  = BornDY*sqrt(pow(qqb_lo_error,2) + asopi2*(pow(qqb_nlo_error,2) + pow(result_gq_NLO.error,2)));

	  if(qcdorder==1)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/Q) << "\t" << (muf/Q) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << std::endl;
	    }
      	}
      if(qcdorder>=2)
      	{
	  // alphaS(mur) and a mb(mur) at NNLO
	  asopi   = as_n3loxs(mur, 2, asopimz);
	  asopi2 = asopi*asopi;
      	  asopi4 = asopi2*asopi2;
      	  logmu1 = log(mur2/muf2);

      	  xsnnlo_result = BornDY*
      	    (qqb_lo_result +
      	     asopi*(qqb_nlo_result + result_gq_NLO.integral) +
      	     asopi2*(qqb_nnlo_result + result_gq_NNLO.integral + result_gg_NNLO.integral +
      		     result_qq_NNLO.integral + result_qQq_NNLO.integral + result_qQqb_NNLO.integral +
      		     constants::b0*logmu1*(qqb_nlo_result + result_gq_NLO.integral))
      	     );
      	  xsnnlo_error  = BornDY*
      	    sqrt(pow(qqb_lo_error,2) +
      		 asopi2*(pow(qqb_nlo_error,2) + pow(result_gq_NLO.error,2)) +
      		 asopi4*(pow(qqb_nnlo_error,2) + pow(result_gq_NNLO.error,2) +
      			 pow(result_gg_NNLO.error,2) + pow(result_qq_NNLO.error,2) +
      			 pow(result_qQq_NNLO.error,2) + pow(result_qQqb_NNLO.error,2) +
      			 pow(constants::b0*logmu1,2)*(pow(qqb_nlo_error,2) + pow(result_gq_NLO.error,2)))
      		 );

	  if(qcdorder==2)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/Q) << "\t" << (muf/Q) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << std::endl;
	    }
      	}
      if(qcdorder==3)
      	{
	  // alphaS(mur) and a mb(mur) at NNLO
	  asopi   = as_n3loxs(mur, 3, asopimz);
	  asopi2 = asopi*asopi;
      	  asopi3 = asopi*asopi2;
	  asopi4 = asopi2*asopi2;
      	  asopi6 = asopi3*asopi3;

      	  logmu2 = logmu1*logmu1;

      	  xsn3lo_result = BornDY*
      	    (qqb_lo_result +
      	     asopi*(qqb_nlo_result + result_gq_NLO.integral) +
      	     asopi2*(qqb_nnlo_result + result_gq_NNLO.integral + result_gg_NNLO.integral +
      		     result_qq_NNLO.integral + result_qQq_NNLO.integral + result_qQqb_NNLO.integral +
      		     constants::b0*logmu1*(qqb_nlo_result + result_gq_NLO.integral)) +
      	     asopi3*(qqb_n3lo_result + result_gq_N3LO.integral + result_gg_N3LO.integral +
      		     result_qq_N3LO.integral + result_qQq_N3LO.integral + result_qQqb_N3LO.integral +
      		     (pow(constants::b0,2)*logmu2 + constants::b1*logmu1)*(qqb_nlo_result + result_gq_NLO.integral) +
      		     2*constants::b0*logmu1*
      		     (qqb_nnlo_result + result_gq_NNLO.integral + result_gg_NNLO.integral +
      		      result_qq_NNLO.integral + result_qQq_NNLO.integral + result_qQqb_NNLO.integral))
      	     );
      	  xsn3lo_error = BornDY*
      	    sqrt(pow(qqb_lo_error,2) +
      		 asopi2*(pow(qqb_nlo_error,2) + pow(result_gq_NLO.error,2)) +
      		 asopi4*(pow(qqb_nnlo_error,2) + pow(result_gq_NNLO.error,2) + pow(result_gg_NNLO.error,2) +
      			 pow(result_qq_NNLO.error,2) + pow(result_qQq_NNLO.error,2) + pow(result_qQqb_NNLO.error,2) +
      			 pow(constants::b0*logmu1,2)*(pow(qqb_nlo_error,2) + pow(result_gq_NLO.error,2))) +
      		 asopi6*(pow(qqb_n3lo_error,2) + pow(result_gq_N3LO.error,2) + pow(result_gg_N3LO.error,2) +
      			 pow(result_qq_N3LO.error,2) + pow(result_qQq_N3LO.error,2) + pow(result_qQqb_N3LO.error,2) +
      			 pow((pow(constants::b0,2)*logmu2 + constants::b1*logmu1),2)*(pow(qqb_nlo_error,2) + pow(result_gq_NLO.error,2)) +
      			 pow(2*constants::b0*logmu1,2)*
      			 (pow(qqb_nnlo_error,2) + pow(result_gq_NNLO.error,2) + pow(result_gg_NNLO.error,2) +
      			  pow(result_qq_NNLO.error,2) + pow(result_qQq_NNLO.error,2) + pow(result_qQqb_NNLO.error,2)))
      	     );

	  fa << std::fixed << std::setprecision(3) << (mur/Q) << "\t" << (muf/Q) << "\t"
	     << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result << "\t" << xsn3lo_result
	     << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << "\t" << xsn3lo_error << std::endl;
      	}
	}
      delete basepdf;
      fa.close();
      return 0;
    }
}
