#include <iostream>
#include <fstream>
//#include <cmath> // sin, cos, exp
#include "qmc.hpp"

#include <stdlib.h>
#include <cstring>

#include "LHAPDF/LHAPDF.h"

// Auxiliary functions for the integration
#include "bbh_functions.h"

// Global constants and QCD parameters
#include "constants.h"

// Header for the routines alphaS(muR) and mb_msbar(muR)
#include "alphaS.h"

struct {
  double s;
  double xmuf;
  double scalemuF0;
  const LHAPDF::PDF* pdf;
} global_param;

#include "pdfpar.h"
struct parampdf_struc parampdf;

//const double scalemuF0 = (constants::MH + 2*constants::Mbpole)/4.0;
//const double scalemuR0 = constants::MH;

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

// b-bbar partonic channel
struct functor_delta_t  {
  unsigned long long int number_of_integration_variables = 1;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = delta_bbh(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_delta;


struct functor_PlusConst_t  {
  unsigned long long int number_of_integration_variables = 1;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = PlusConst_bbh(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_PlusConst;

struct functor_PlusInt1_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = PlusInt1_bbh(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_PlusInt1;

struct functor_PlusInt2_t  {
  unsigned long long int number_of_integration_variables = 2;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = PlusInt2_bbh(x, global_param.s, muf, k, global_param.pdf);

    return integrand;

  }
} functor_PlusInt2;

struct functor_RegNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bbb_regular_nlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegNLO;

struct functor_RegNNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bbb_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegNNLO;

struct functor_RegN3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bbb_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegN3LO;


// g-b partonic channel
struct functor_bg_NLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bg_regular_nlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bg_NLO;

struct functor_bg_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bg_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bg_NNLO;

struct functor_bg_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bg_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bg_N3LO;


// g-g partonic channel
struct functor_gg_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = gg_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gg_NNLO;

struct functor_gg_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = gg_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gg_N3LO;


// q-g (+ qbar-g) partonic channel
struct functor_qg_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = qg_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qg_N3LO;


// b-b (+ bbar-bbar) partonic channel
struct functor_bb_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bb_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bb_NNLO;

struct functor_bb_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bb_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bb_N3LO;


// b-q (+ bbar-qbar) partonic channel
struct functor_bq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bq_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bq_NNLO;

struct functor_bq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bq_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bq_N3LO;


// b-qbar (+ bbar-q) partonic channel
struct functor_bqbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bqb_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bqbar_NNLO;

struct functor_bqbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = bqb_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_bqbar_N3LO;


// q-qbar partonic channel
struct functor_qqbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = qqb_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqbar_NNLO;

struct functor_qqbar_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*global_param.scalemuF0;

    integrand = qqb_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqbar_N3LO;


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
      std::cout << "Usage:  " << argv[0] << " a b c d e f g h i j k l m n o p (q) with:" << std::endl;
      std::cout << "a:  Lattice size (integer)" << std::endl;
      std::cout << "b:  Seed (integer)" << std::endl;
      std::cout << "c:  QCD order (integer, between O for LO and 3 for N3LO)" << std::endl;
      std::cout << "d:  p-p (0) or p-pbar (1) collider" << std::endl;
      std::cout << "e:  Hadronic energy in TeV (double)" << std::endl;
      std::cout << "f:  x_muf so that mu_F = x_muf*mu_F0 (double)" << std::endl;
      std::cout << "g:  mu_F0: central factorization scale (double); if set to -1, default value is (MH + Mb)/4" << std::endl;
      std::cout << "h:  x_mur so that mu_R = x_mur*mu_R0 (double)" << std::endl;
      std::cout << "i:  mu_R0: central renormalization scale (double); if set to -1, default value is MH" << std::endl;
      std::cout << "j:  PDF set (string)" << std::endl;
      std::cout << "k:  PDF member (integer)" << std::endl;
      std::cout << "l:  Z mass in GeV (double)" << std::endl;
      std::cout << "m:  Higgs mass in GeV (double)" << std::endl;
      std::cout << "n:  Bottom-quark pole mass in GeV (double)" << std::endl;
      std::cout << "o:  Bottom-quark MSbar mass mb(mb) in GeV (double)" << std::endl;
      std::cout << "p:  Vacuum expectation value in GeV (double)" << std::endl;
      std::cout << "q:  --scale: optional flag to calculate various mu_R predictions. If absent, mu_R = x_mur*mu_R0" << std::endl;
      return 0;
    }
  if(argc < 17)
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

      std::string energyheader = argv[5];
      removeTrailingCharacters(energyheader, '0');
      double energy = std::stod(energyheader); // energy in TeV
      double s;
      s = energy*energy*1.e6;
      double xmuf = atof(argv[6]);

      double muf0 = atof(argv[7]);
      double xmur = atof(argv[8]);
      double mur0 = atof(argv[9]);
      int muf_flag;
      int mur_flag;
      double scalemuF0;
      double scalemuR0;

      // init PDF set
      const std::string setname = argv[10];
      const int setimem = atoi(argv[11]);
      const LHAPDF::PDF* basepdf = LHAPDF::mkPDF( setname, setimem);
      LHAPDF::setVerbosity(0); // default is 1;

      constants::MZ     = atof(argv[12]);
      constants::MH     = atof(argv[13]);
      constants::Mb     = atof(argv[14]);
      constants::Mbmb   = atof(argv[15]);
      constants::vev    = atof(argv[16]);
      if(muf0 == -1)
	{
	  scalemuF0     = (constants::MH + 2*constants::Mb)/4.0;
	  muf_flag      = 0;
	}
      else
	{
	  scalemuF0     = muf0;
	  muf_flag      = 1;
	}
      if(mur0 == -1)
	{
	  scalemuR0     = constants::MH;
	  mur_flag      = 0;
	}
      else
	{
	  scalemuR0     = mur0;
	  mur_flag      = 1;
	}
      global_param.scalemuF0  = scalemuF0;

      // init parameters for all functors
      global_param.s       = s;
      global_param.xmuf    = xmuf;
      global_param.pdf     = basepdf;

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
      double bbbar_lo_result,bbbar_lo_error;

      integrators::result<double> resultPlusConst;
      integrators::result<double> resultPlusInt1;
      integrators::result<double> resultPlusInt2;
      integrators::result<double> resultRegNLO;
      integrators::result<double> result_bg_NLO;
      double bbbar_nlo_result,bbbar_nlo_error;

      integrators::result<double> resultRegNNLO;
      integrators::result<double> result_bg_NNLO;
      integrators::result<double> result_gg_NNLO;
      integrators::result<double> result_bb_NNLO;
      integrators::result<double> result_bq_NNLO;
      integrators::result<double> result_bqbar_NNLO;
      integrators::result<double> result_qqbar_NNLO;
      double bbbar_nnlo_result,bbbar_nnlo_error;

      integrators::result<double> resultRegN3LO;
      integrators::result<double> result_bg_N3LO;
      integrators::result<double> result_gg_N3LO;
      integrators::result<double> result_bb_N3LO;
      integrators::result<double> result_bq_N3LO;
      integrators::result<double> result_bqbar_N3LO;
      integrators::result<double> result_qqbar_N3LO;
      integrators::result<double> result_qg_N3LO;
      double bbbar_n3lo_result,bbbar_n3lo_error;


      // Perform all the required QMC integrations

      if(qcdorder>=0)
	{
	  real_integrator.minn = lattice/10;
	  functor_delta.k = 0;
	  resultdelta = real_integrator.integrate(functor_delta);
	  bbbar_lo_result = resultdelta.integral;
	  bbbar_lo_error  = resultdelta.error;

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
	      bbbar_nlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNLO.integral;
	      bbbar_nlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
	      			    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNLO.error,2));
	      result_bg_NLO = real_integrator.integrate(functor_bg_NLO);

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
		  bbbar_nnlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNNLO.integral;
		  bbbar_nnlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		  			 pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNNLO.error,2));
		  result_bg_NNLO   = real_integrator.integrate(functor_bg_NNLO);
		  result_gg_NNLO   = real_integrator.integrate(functor_gg_NNLO);
		  result_bb_NNLO   = real_integrator.integrate(functor_bb_NNLO);
		  result_bq_NNLO   = real_integrator.integrate(functor_bq_NNLO);
		  result_bqbar_NNLO   = real_integrator.integrate(functor_bqbar_NNLO);
		  result_qqbar_NNLO   = real_integrator.integrate(functor_qqbar_NNLO);

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
		      bbbar_n3lo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegN3LO.integral;
		      bbbar_n3lo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		      			    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegN3LO.error,2));
		      result_bg_N3LO   = real_integrator.integrate(functor_bg_N3LO);
		      real_integrator.minn = lattice;
		      result_gg_N3LO   = real_integrator.integrate(functor_gg_N3LO);
		      result_bb_N3LO   = real_integrator.integrate(functor_bb_N3LO);
		      result_bq_N3LO   = real_integrator.integrate(functor_bq_N3LO);
		      result_bqbar_N3LO   = real_integrator.integrate(functor_bqbar_N3LO);
		      result_qqbar_N3LO   = real_integrator.integrate(functor_qqbar_N3LO);
		      result_qg_N3LO   = real_integrator.integrate(functor_qg_N3LO);
		    }
		}
	    }
	}

      // Building the cross section: adding back alphaS, evolve to muR, add Born normalization
      // Result of this code: xs(b bbar -> H + X) in pb

      double BornbbH;
      double mbatmur;

      double muf  = xmuf*scalemuF0;
      double muf2 = muf*muf;
      double mur;
      double mur2;

      double asmbopi;
      double asopi;
      double asopi2;
      double asopi3;
      double asopi4;
      double asopi6;

      double logmu1;
      double logmu2;
      double logmu3;

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
      

      int imax;
      double dxmur;

      double asopimz = (basepdf->alphasQ(constants::MZ))/constants::Pi;

      if(argc>=18 && std::strcmp(argv[17],"--scale") == 0)
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
	  if(imax==1)
	    {
	      filename << "bbH_xs_pp_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << "_mur" << xmur << ".txt";
	    }
	  else
	    {
	      filename << "bbH_xs_pp_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	    }
	  header = "# Inclusive cross section for SM Higgs production in bottom-quark fusion xs(b bbbar -> H) in the 5FS, p-p collider, sqrt(S) = ";
	}
      else
	{
	  if(imax==1)
	    {
	      filename << "bbH_xs_ppbar_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << "_mur" << xmur << ".txt";
	    }
	  else
	    {
	      filename << "bbH_xs_ppbar_" << energyheader << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	    }
	  header = "# Inclusive cross section for SM Higgs production in bottom-quark fusion xs(b bbbar -> H) in the 5FS, p-pbar collider, sqrt(S) = ";
	}

      
      filename >> finalfile;
      std::ofstream fa(finalfile);

      if(muf_flag == 0)
	{
	  if(mur_flag == 0)
	    {
	      header += energyheader +
		" TeV, central factorization scale mu_F0 = (MH+2*Mb)/4, renormalization scale mu_R0 = MH\n# mu_R/mu_R0\t" +
		"mu_F/mu_F0\t";
	    }
	  else
	    {
	      mur0header = std::to_string(scalemuR0);
	      removeTrailingCharacters(mur0header, '0');
	      header += energyheader +
		" TeV, central factorization scale mu_F0 = (MH+2*Mb)/4, renormalization scale mu_R0 = " +
		mur0header +
		" GeV\n# mu_R/mu_R0\t" +
		"mu_F/mu_F0\t";
	    }
	}
      else
	{
	  muf0header = std::to_string(scalemuF0);
	  removeTrailingCharacters(muf0header, '0');
	  if(mur_flag == 0)
	    {
	      header += energyheader +
		" TeV, central factorization scale mu_F0 = " +
		muf0header +
		" GeV, renormalization scale mu_R0 = MH\n# mu_R/mu_R0\t" +
		"mu_F/mu_F0\t";
	    }
	  else
	    {
	      mur0header = std::to_string(scalemuR0);
	      removeTrailingCharacters(mur0header, '0');
	      header += energyheader +
		" TeV, central factorization scale mu_F0 = " +
		muf0header +
		" GeV, renormalization scale mu_R0 = " +
		mur0header +
		" GeV\n# mu_R/mu_R0\t" + "mu_F/mu_F0\t";
	    }
	}


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

      for(int i = 0; i<imax; i++)
	{
	  if(imax!=1)
	    {
	      xmur = 0.5 + i*dxmur;
	    }
	  mur  = xmur*scalemuR0;
	  mur2 = mur*mur;


      if(qcdorder>=0)
      	{

	  // alphaS(mur) and mb(mur) at LO
	  asopi   = as_n3loxs(mur, 0, asopimz);
	  asmbopi = as_n3loxs(constants::Mbmb, 0, asopimz);
	  mbatmur = mb_n3loxs(mur, 0, constants::Mbmb, constants::Mbmb, asmbopi);
	  //std::cout << "\n\t The program uses mb_LO(muR)   = " << mbatmur << " GeV" << std::endl;

	  BornbbH = constants::gevtopb*constants::Pi*mbatmur*mbatmur
	    /(2*constants::Nc*constants::vev*constants::vev*constants::MH*constants::MH);

      	  xslo_result = BornbbH*bbbar_lo_result;
      	  xslo_error  = BornbbH*bbbar_lo_error;
	  if(qcdorder==0)
	    {
	      fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << std::scientific << xslo_error << std::endl;
	    }
      	}
      if(qcdorder>=1)
      	{
	  // alphaS(mur) and mb(mur) at NLO
	  asopi   = as_n3loxs(mur, 1, asopimz);
	  asmbopi = as_n3loxs(constants::Mbmb, 1, asopimz);
	  mbatmur = mb_n3loxs(mur, 1, constants::Mbmb, constants::Mbmb, asmbopi);
	  //std::cout << "\t The program uses mb_NLO(muR)  = " << mbatmur << " GeV" << std::endl;

	  BornbbH = constants::gevtopb*constants::Pi*mbatmur*mbatmur
	    /(2*constants::Nc*constants::vev*constants::vev*constants::MH*constants::MH);

	  asopi2 = asopi*asopi;
      	  logmu1 = log(mur2/muf2);

      	  xsnlo_result = BornbbH*
	    (bbbar_lo_result +
	     asopi*(bbbar_nlo_result + result_bg_NLO.integral -2*constants::byuk0*logmu1*bbbar_lo_result)
	     );
      	  xsnlo_error  = BornbbH*
	    sqrt(pow(bbbar_lo_error,2) +
		 asopi2*(pow(bbbar_nlo_error,2) + pow(result_bg_NLO.error,2) + 4*pow(constants::byuk0*logmu1*bbbar_lo_error,2))
		 );

	  if(qcdorder==1)
	    {
	      fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << std::endl;
	    }
      	}
      if(qcdorder>=2)
      	{
	  // alphaS(mur) and mb(mur) at NNLO
	  asopi   = as_n3loxs(mur, 2, asopimz);
	  asmbopi = as_n3loxs(constants::Mbmb, 2, asopimz);
	  mbatmur = mb_n3loxs(mur, 2, constants::Mbmb, constants::Mbmb, asmbopi);
	  //std::cout << "\t The program uses mb_NNLO(muR) = " << mbatmur << " GeV" << std::endl;

	  BornbbH = constants::gevtopb*constants::Pi*mbatmur*mbatmur
	    /(2*constants::Nc*constants::vev*constants::vev*constants::MH*constants::MH);

	  asopi2 = asopi*asopi;
      	  asopi4 = asopi2*asopi2;
	  logmu2 = logmu1*logmu1;

      	  xsnnlo_result = BornbbH*
      	    (bbbar_lo_result +
      	     asopi*(bbbar_nlo_result + result_bg_NLO.integral -2*constants::byuk0*logmu1*bbbar_lo_result) +
      	     asopi2*(bbbar_nnlo_result + result_bg_NNLO.integral + result_gg_NNLO.integral + result_bb_NNLO.integral + 
      		     result_bq_NNLO.integral + result_bqbar_NNLO.integral + result_qqbar_NNLO.integral -
		     ((constants::b0-2*constants::byuk0)*constants::byuk0*logmu2 + 2*constants::byuk1*logmu1)*bbbar_lo_result +
		     (constants::b0-2*constants::byuk0)*logmu1*(bbbar_nlo_result + result_bg_NLO.integral)
		     )
      	     );
      	  xsnnlo_error  = BornbbH*
      	    sqrt(pow(bbbar_lo_error,2) +
      		 asopi2*(pow(bbbar_nlo_error,2) + pow(result_bg_NLO.error,2) + 4*pow(constants::byuk0*logmu1*bbbar_lo_error,2)) +
      		 asopi4*(pow(bbbar_nnlo_error,2) + pow(result_bg_NNLO.error,2) + pow(result_gg_NNLO.error,2) + pow(result_bb_NNLO.error,2) +
			 pow(result_bq_NNLO.error,2) + pow(result_bqbar_NNLO.error,2) + pow(result_qqbar_NNLO.error,2) +
			 pow(((constants::b0-2*constants::byuk0)*constants::byuk0*logmu2 + 2*constants::byuk1*logmu1)*bbbar_lo_error,2) +
			 pow((constants::b0-2*constants::byuk0)*logmu1,2)*(pow(bbbar_nlo_error,2) + pow(result_bg_NLO.error,2))
			 )
      		 );

	  if(qcdorder==2)
	    {
	      fa << std::fixed << std::setprecision(3) << xmur<< "\t" << xmuf << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result
		 << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << std::endl;
	    }
      	}
      if(qcdorder==3)
      	{
	  // alphaS(mur) and mb(mur) at N3LO
	  asopi   = as_n3loxs(mur, 3, asopimz);
	  asmbopi = as_n3loxs(constants::Mbmb, 3, asopimz);
	  mbatmur = mb_n3loxs(mur, 3, constants::Mbmb, constants::Mbmb, asmbopi);
	  //std::cout << "\t The program uses mb_N3LO(muR) = " << mbatmur << " GeV" << std::endl;

	  BornbbH = constants::gevtopb*constants::Pi*mbatmur*mbatmur
	    /(2*constants::Nc*constants::vev*constants::vev*constants::MH*constants::MH);

	  asopi2 = asopi*asopi;
      	  asopi3 = asopi*asopi2;
	  asopi4 = asopi2*asopi2;
      	  asopi6 = asopi3*asopi3;
	  logmu3 = logmu1*logmu2;

      	  xsn3lo_result = BornbbH*
      	    (bbbar_lo_result +
      	     asopi*(bbbar_nlo_result + result_bg_NLO.integral -2*constants::byuk0*logmu1*bbbar_lo_result) +
      	     asopi2*(bbbar_nnlo_result + result_bg_NNLO.integral + result_gg_NNLO.integral + result_bb_NNLO.integral +
		     result_bq_NNLO.integral + result_bqbar_NNLO.integral + result_qqbar_NNLO.integral -
		     ((constants::b0-2*constants::byuk0)*constants::byuk0*logmu2 + 2*constants::byuk1*logmu1)*bbbar_lo_result +
		     (constants::b0-2*constants::byuk0)*logmu1*(bbbar_nlo_result + result_bg_NLO.integral)
		     ) +
      	     asopi3*(bbbar_n3lo_result + result_bg_N3LO.integral + result_gg_N3LO.integral +
      	     	     result_bb_N3LO.integral + result_bq_N3LO.integral + result_bqbar_N3LO.integral +
	     	     result_qqbar_N3LO.integral + result_qg_N3LO.integral -
	     	     (2*constants::byuk2*logmu1 + constants::byuk0*(constants::b1-2*constants::byuk1)*logmu2 +
	     	      2*(constants::b0-constants::byuk0)*(constants::byuk1*logmu2+constants::byuk0*logmu3*(constants::b0-2*constants::byuk0)/3.0))
		     *bbbar_lo_result +
	     	     ((constants::b1-2*constants::byuk1)*logmu1 + (constants::b0-constants::byuk0)*(constants::b0-2*constants::byuk0)*logmu2)*
	     	     (bbbar_nlo_result + result_bg_NLO.integral) +
	     	     2*(constants::b0-constants::byuk0)*logmu1*
	     	     (bbbar_nnlo_result + result_bg_NNLO.integral + result_gg_NNLO.integral + result_bb_NNLO.integral +
	     	      result_bq_NNLO.integral + result_bqbar_NNLO.integral + result_qqbar_NNLO.integral)
	     	     )
      	     );
      	  xsn3lo_error = BornbbH*
      	    sqrt(pow(bbbar_lo_error,2) +
      		 asopi2*(pow(bbbar_nlo_error,2) + pow(result_bg_NLO.error,2) + 4*pow(constants::byuk0*logmu1*bbbar_lo_error,2)) +
      		 asopi4*(pow(bbbar_nnlo_error,2) + pow(result_bg_NNLO.error,2) + pow(result_gg_NNLO.error,2) +
      			 pow(result_bb_NNLO.error,2) + pow(result_bq_NNLO.error,2) + pow(result_bqbar_NNLO.error,2) +
			 pow(result_qqbar_NNLO.error,2) +
			 pow(((constants::b0-2*constants::byuk0)*constants::byuk0*logmu2 + 2*constants::byuk1*logmu1)*bbbar_lo_error,2) +
			 pow((constants::b0-2*constants::byuk0)*logmu1,2)*(pow(bbbar_nlo_error,2) + pow(result_bg_NLO.error,2))
			 ) +
      		 asopi6*(pow(bbbar_n3lo_error,2) + pow(result_bg_N3LO.error,2) + pow(result_gg_N3LO.error,2) +
      		 	 pow(result_bb_N3LO.error,2) + pow(result_bq_N3LO.error,2) + pow(result_bqbar_N3LO.error,2) +
		 	 pow(result_qqbar_N3LO.error,2) + pow(result_qg_N3LO.error,2) +
			 pow((2*constants::byuk2*logmu1 + constants::byuk0*(constants::b1-2*constants::byuk1)*logmu2 +
			      2*(constants::b0-constants::byuk0)*(constants::byuk1*logmu2+constants::byuk0*logmu3*(constants::b0-2*constants::byuk0)/3.0))
			     *bbbar_lo_error,2) +
			 pow(((constants::b1-2*constants::byuk1)*logmu1 + (constants::b0-constants::byuk0)*(constants::b0-2*constants::byuk0)*logmu2),2)*
			 (pow(bbbar_nlo_error,2) + pow(result_bg_NLO.error,2)) + 4*pow((constants::b0-constants::byuk0)*logmu1,2)*
			 (pow(bbbar_nnlo_error,2) + pow(result_bg_NNLO.error,2) + pow(result_gg_NNLO.error,2) + pow(result_bb_NNLO.error,2) +
			  pow(result_bq_NNLO.error,2) + pow(result_bqbar_NNLO.error,2) + pow(result_qqbar_NNLO.error,2))
			 )
		 );

	  fa << std::fixed << std::setprecision(3) << xmur << "\t" << xmuf << "\t"
	     << std::setprecision(18) << xslo_result << "\t" << xsnlo_result << "\t" << xsnnlo_result << "\t" << xsn3lo_result
	     << "\t" << std::scientific << xslo_error << "\t" << xsnlo_error << "\t" << xsnnlo_error << "\t" << xsn3lo_error << std::endl;
      	}
	}
      delete basepdf;
      fa.close();
      return 0;
    }
}
