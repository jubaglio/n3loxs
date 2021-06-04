#include <iostream>
#include <fstream>
//#include <cmath> // sin, cos, exp
#include "qmc.hpp"

#include <stdlib.h>
#include <cstring>

#include "LHAPDF/LHAPDF.h"

// Auxiliary functions for the integration
#include "dy_functions_wh_dyn.h"

// Global constants and QCD parameters
#include "constants.h"

struct {
  double s;
  double xmuf;
  double xmur;
  double asopimz;
  const LHAPDF::PDF* pdf;
} global_param;

#include "pdfpar_w.h"
struct parampdf_struc parampdf;

double constants::MW;
double constants::MZ;
double constants::MH;
double constants::Mb;
double constants::Mt;
double constants::Mbmb;

double constants::vev;
double constants::alphainv;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Functors for the integration
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// d-ubar partonic channel
struct functor_delta_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];

    integrand = delta(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, k, global_param.pdf);

    return integrand;

  }
} functor_delta;


struct functor_PlusConst_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];

    integrand = PlusConst(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, k, global_param.pdf);

    return integrand;

  }
} functor_PlusConst;

struct functor_PlusInt1_t  {
  unsigned long long int number_of_integration_variables = 3;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = PlusInt1(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, k, global_param.pdf);

    return integrand;

  }
} functor_PlusInt1;

struct functor_PlusInt2_t  {
  unsigned long long int number_of_integration_variables = 3;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = PlusInt2(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, k, global_param.pdf);

    return integrand;

  }
} functor_PlusInt2;

struct functor_RegNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = dub_regular_nlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_RegNLO;

struct functor_RegNNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = dub_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_RegNNLO;

struct functor_RegN3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = dub_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_RegN3LO;


// g-ubar partonic channel
struct functor_gubar_NLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = gub_regular_nlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_gubar_NLO;

struct functor_gubar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = gub_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_gubar_NNLO;

struct functor_gubar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = gub_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_gubar_N3LO;


// g-g partonic channel
struct functor_gg_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = gg_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_gg_NNLO;

struct functor_gg_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = gg_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_gg_N3LO;


// g-dbar partonic channel
struct functor_gdbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = gdb_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_gdbar_N3LO;


// c-ubar partonic channel
struct functor_cubar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = cub_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_cubar_NNLO;

struct functor_cubar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = cub_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_cubar_N3LO;


// q-qbar partonic channel
struct functor_qqbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qqb_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qqbar_NNLO;

struct functor_qqbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qqb_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qqbar_N3LO;


// q-q partonic channel
struct functor_qq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qq_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qq_NNLO;

struct functor_qq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qq_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qq_N3LO;


// q-qprime partonic channel
struct functor_qqprime_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qqprime_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qqprime_NNLO;

struct functor_qqprime_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qqprime_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qqprime_N3LO;


// qbar-qprimebar partonic channel
struct functor_qbarqprimebar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qbqprimeb_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qbarqprimebar_NNLO;

struct functor_qbarqprimebar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand = qbqprimeb_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_qbarqprimebar_N3LO;


// ubar-cbar + d-s partonic channel
struct functor_ubarcbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand =
      ubcb_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf) +
      ds_regular_nnlo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_ubarcbar_NNLO;

struct functor_ubarcbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 3;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];

    integrand =
      ubcb_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf) +
      ds_regular_n3lo(x, global_param.s, global_param.xmuf, global_param.xmur, global_param.asopimz, global_param.pdf);

    return integrand;

  }
} functor_ubarcbar_N3LO;


static void removeTrailingCharacters(std::string &str, const char charToRemove) {
  str.erase (str.find_last_not_of(charToRemove) + 1, std::string::npos );
  double numb = std::stod(str);
  if(int(numb)/numb==1)
    {
      str.erase (str.find_last_not_of('.') + 1, std::string::npos );
    }
}

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
      std::cout << "Usage:  " << argv[0] << " a b c d e f g h i j k l m n with:" << std::endl;
      std::cout << "a:  Lattice size (integer)" << std::endl;
      std::cout << "b:  Seed (integer)" << std::endl;
      std::cout << "c:  QCD order (integer, between O for LO and 3 for N3LO)" << std::endl;
      std::cout << "d:  p-p (0) or p-pbar (1) collider" << std::endl;
      std::cout << "e:  Hadronic energy in TeV (double)" << std::endl;
      std::cout << "f:  W+ (1) or W- (-1) production (integer)" << std::endl;
      std::cout << "g:  x_muf so that mu_F = x_muf*M_(HW) (double)" << std::endl;
      std::cout << "h:  x_mur so that mu_R = x_mur*M_(HW) (double)" << std::endl;
      std::cout << "i:  PDF set (string)" << std::endl;
      std::cout << "j:  PDF member (integer)" << std::endl;
      std::cout << "k:  W mass in GeV (double)" << std::endl;
      std::cout << "l:  Z mass in GeV (double)" << std::endl;
      std::cout << "m:  Higgs mass in GeV (double)" << std::endl;
      std::cout << "n:  Vacuum expectation value in GeV (double)" << std::endl;
      return 0;
    }
  if(argc < 14)
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

      std::string energyheader = argv[5];
      removeTrailingCharacters(energyheader, '0');
      double energy = std::stod(energyheader); // energy in TeV
      double s;
      s = energy*energy*1.e6;
      int wchoice = atoi(argv[6]);
      wchoice = -wchoice; // internal conversion: all the contributions are written for W- by default.
      parampdf.wchoice = wchoice;
      double xmuf = atof(argv[7]);
      double xmur = atof(argv[8]);

      // init PDF set
      const std::string setname = argv[9];
      const int setimem = atoi(argv[10]);
      const LHAPDF::PDF* basepdf = LHAPDF::mkPDF( setname, setimem);
      LHAPDF::setVerbosity(0); // default is 1;

      constants::MW     = atof(argv[11]);
      constants::MZ     = atof(argv[12]);
      constants::MH     = atof(argv[13]);
      constants::vev    = atof(argv[14]);

      // init parameters for all functors
      global_param.s       = s;
      global_param.xmuf    = xmuf;
      global_param.xmur    = xmur;
      global_param.pdf     = basepdf;
      double asopimz = (basepdf->alphasQ(constants::MZ))/constants::Pi;
      global_param.asopimz    = asopimz;

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
	      dubar_nlo_error  = pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNLO.error,2);
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
		  real_integrator.minn = 10*lattice;
		  resultPlusInt1 = real_integrator.integrate(functor_PlusInt1);
		  real_integrator.minn = lattice;
		  resultPlusInt2 = real_integrator.integrate(functor_PlusInt2);
		  real_integrator.minn = 10*lattice;
		  resultRegNNLO    = real_integrator.integrate(functor_RegNNLO);
		  dubar_nnlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNNLO.integral;
		  dubar_nnlo_error  = pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNNLO.error,2);
		  result_gubar_NNLO   = real_integrator.integrate(functor_gubar_NNLO);
		  real_integrator.minn = lattice;
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
		      real_integrator.minn = 10*lattice;
		      resultPlusInt1 = real_integrator.integrate(functor_PlusInt1);
		      real_integrator.minn = lattice;
		      resultPlusInt2 = real_integrator.integrate(functor_PlusInt2);
		      real_integrator.minn = 10*lattice;
		      resultRegN3LO    = real_integrator.integrate(functor_RegN3LO);
		      dubar_n3lo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegN3LO.integral;
		      dubar_n3lo_error  = pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
			pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegN3LO.error,2);
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

      // Building the cross section: adding back alphaS, evolve to muR, add Born normalization
      // Result of this code: xs(p p / p pbar -> W+/W- + H + X) in pb

      double BornDY = constants::gevtopb*constants::MW*constants::MW*constants::MW*constants::MW
      	/(48*constants::Pi*constants::Nc*constants::vev*constants::vev*constants::vev*constants::vev);

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
      std::string muf0header;
      std::string mur0header;

      if(collider==0)
	{
	  if(wchoice==1)
	    {
	      filename << "WminusH_xs_pp_dyn_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << "_mur" << xmur << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p p -> W- H), sqrt(S) = ";
	    }
	  else
	    {
	      filename << "WplusH_xs_pp_dyn_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << "_mur" << xmur << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p p -> W+ H), sqrt(S) = ";
	    }
	}
      else
	{
	  if(wchoice==1)
	    {
	      filename << "WminusH_xs_ppbar_dyn_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << "_mur" << xmur << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p pbar -> W- H), sqrt(S) = ";
	    }
	  else
	    {
	      filename << "WplusH_xs_ppbar_dyn_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << "_mur" << xmur << ".txt";
	      header = "# Standard Model Higgs-strahlung cross section xs(p pbar -> W+ H), sqrt(S) = ";
	    }
	}
      
      
      filename >> finalfile;
      std::ofstream fa(finalfile);

      energyheader = std::to_string(energy);
      removeTrailingCharacters(energyheader, '0');

      header += energyheader +
	" TeV, dynamical central factorization and renormalization scales mu_F0 = mu_R0 = M_(HW)\n# mu_R/mu_R0\t" +
	"mu_F/mu_F0\t";

      switch(qcdorder)
	{
	case 0:
	  fa << header << "xs_LO (pb)\t" << "num error (respective xs)" << std::endl;
	  break;
	case 1:
	  fa << header << "xs_LO (pb)\t" << "xs_NLO (pb)\t" << "num error (respective xs)" << std::endl;
	  break;
	case 2:
	  fa << header << "xs_LO (pb)\t" << "xs_NLO (pb)\t" << "xs_NNLO (pb)\t" << "num error (respective xs)" << std::endl;
	  break;
	case 3:
	  fa << header << "xs_LO (pb)\t" << "xs_NLO (pb)\t" << "xs_NNLO (pb)\t" << "xs_N3LO (pb)\t" << "num error (respective xs)" << std::endl;
	  break;
	}
      
      if(qcdorder>=0)
      	{
      	  xslo_result = BornDY*dubar_lo_result;
      	  xslo_error  = BornDY*dubar_lo_error;
	  if(qcdorder==0)
	    {
	      fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << std::scientific << xslo_error << std::endl;
	    }
      	}
      if(qcdorder>=1)
      	{
      	  xsnlo_result = BornDY*(dubar_nlo_result + result_gubar_NLO.integral);
      	  xsnlo_error  = BornDY*sqrt(dubar_nlo_error + pow(result_gubar_NLO.error,2));

	  if(qcdorder==1)
	    {
	      fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << std::endl;
	    }
      	}
      if(qcdorder>=2)
      	{
      	  xsnnlo_result = BornDY*
      	    (dubar_nnlo_result + result_gubar_NNLO.integral + result_gg_NNLO.integral +
	     result_cubar_NNLO.integral + result_qqbar_NNLO.integral + result_qq_NNLO.integral + 
	     result_qqprime_NNLO.integral + result_qbarqprimebar_NNLO.integral + result_ubarcbar_NNLO.integral
      	     );
      	  xsnnlo_error  = BornDY*
      	    sqrt(
		 dubar_nnlo_error + pow(result_gubar_NNLO.error,2) +
		 pow(result_gg_NNLO.error,2) + pow(result_cubar_NNLO.error,2) +
		 pow(result_qqbar_NNLO.error,2) + pow(result_qq_NNLO.error,2) +
		 pow(result_qqprime_NNLO.error,2) + pow(result_qbarqprimebar_NNLO.error,2) + pow(result_ubarcbar_NNLO.error,2)
      		 );

	  if(qcdorder==2)
	    {
	      fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << std::endl;
	    }
      	}
      if(qcdorder==3)
      	{
      	  xsn3lo_result = BornDY*
      	    (dubar_n3lo_result + result_gubar_N3LO.integral + result_gg_N3LO.integral +
	     result_cubar_N3LO.integral + result_qqbar_N3LO.integral + result_qq_N3LO.integral +
	     result_qqprime_N3LO.integral + result_qbarqprimebar_N3LO.integral +
	     result_ubarcbar_N3LO.integral + result_gdbar_N3LO.integral
      	     );
      	  xsn3lo_error = BornDY*
      	    sqrt(
      		 dubar_n3lo_error + pow(result_gubar_N3LO.error,2) + pow(result_gg_N3LO.error,2) +
		 pow(result_cubar_N3LO.error,2) + pow(result_qqbar_N3LO.error,2) + pow(result_qq_N3LO.error,2) +
		 pow(result_qqprime_N3LO.error,2) + pow(result_qbarqprimebar_N3LO.error,2) +
		 pow(result_ubarcbar_N3LO.error,2) + pow(result_gdbar_N3LO.error,2)
		 );

	  fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
	     << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result << "\t" << xsn3lo_result
	     << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << "\t" << xsn3lo_error << std::endl;
      	}
      delete basepdf;
      fa.close();
      return 0;
    }
}
