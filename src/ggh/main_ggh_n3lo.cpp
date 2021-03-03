#include <iostream>
#include <fstream>
//#include <cmath> // sin, cos, exp
#include "qmc.hpp"

#include <stdlib.h>
#include <cstring>

#include "LHAPDF/LHAPDF.h"

// Auxiliary functions for the integration
#include "ggh_functions.h"

// Global constants and QCD parameters
#include "constants.h"

// Header for the routines alphaS(muR) and mb_msbar(muR)
#include "alphaS.h"

// Header for exact LO form factor: the code gives the Born-improved cross sections
#include "exact_lo_ggh.h"

struct {
  double s;
  double xmuf;
  const LHAPDF::PDF* pdf;
} global_param;

#include "pdfpar.h"
struct parampdf_struc parampdf;

const double scale0 = constants::MH/2.0;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Functors for the integration
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// gg partonic channel
struct functor_delta_t  {
  unsigned long long int number_of_integration_variables = 1;
  int k;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    muf = global_param.xmuf*scale0;

    integrand = delta_ggh(x, global_param.s, muf, k, global_param.pdf);

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
    muf = global_param.xmuf*scale0;

    integrand = PlusConst_ggh(x, global_param.s, muf, k, global_param.pdf);

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
    muf = global_param.xmuf*scale0;

    integrand = PlusInt1_ggh(x, global_param.s, muf, k, global_param.pdf);

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
    muf = global_param.xmuf*scale0;

    integrand = PlusInt2_ggh(x, global_param.s, muf, k, global_param.pdf);

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
    muf = global_param.xmuf*scale0;

    integrand = gg_regular_nlo(x, global_param.s, muf, global_param.pdf);

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
    muf = global_param.xmuf*scale0;

    integrand = gg_regular_nnlo(x, global_param.s, muf, global_param.pdf);

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
    muf = global_param.xmuf*scale0;

    integrand = gg_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_RegN3LO;


// g-q partonic channel
struct functor_bg_NLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = gq_regular_nlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gq_NLO;

struct functor_gq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = gq_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gq_NNLO;

struct functor_gq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = gq_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_gq_N3LO;


// q-qbar partonic channel
struct functor_qqbar_NLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = qqb_regular_nlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqbar_NLO;

struct functor_qqbar_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

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
    muf = global_param.xmuf*scale0;

    integrand = qqb_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qqbar_N3LO;


// q-q partonic channel
struct functor_qq_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = qq_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qq_NNLO;

struct functor_qq_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = qq_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_qq_N3LO;


// q1-q2 partonic channel
struct functor_q1q2_NNLO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = q1q2_regular_nnlo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_q1q2_NNLO;

struct functor_q1q2_N3LO_t  {
  unsigned long long int number_of_integration_variables = 2;

  double operator()(double* y) const {
    double x[number_of_integration_variables];
    double muf;
    double integrand;

    x[0] = y[0];
    x[1] = y[1];
    muf = global_param.xmuf*scale0;

    integrand = q1q2_regular_n3lo(x, global_param.s, muf, global_param.pdf);

    return integrand;

  }
} functor_q1q2_N3LO;


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
      std::cout << "Usage:  " << argv[0] << " a b c d e f g h (i) with:" << std::endl;
      std::cout << "a:  Lattice size (integer)" << std::endl;
      std::cout << "b:  Seed (integer)" << std::endl;
      std::cout << "c:  QCD order (integer, between O for LO and 3 for N3LO)" << std::endl;
      std::cout << "d:  p-p (0) or p-pbar (1) collider" << std::endl;
      std::cout << "e:  Hadronic energy in TeV (double)" << std::endl;
      std::cout << "f:  x_muf so that mu_F = x_muf*MH/2 (double)" << std::endl;
      std::cout << "g:  PDF set (string)" << std::endl;
      std::cout << "h:  PDF member (integer)" << std::endl;
      std::cout << "i:  --scale: optional flag to calculate various mu_R predictions. If absent, mu_R = mu_F" << std::endl;
      return 0;
    }
  if(argc < 9)
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
      double xmuf = atof(argv[6]);

      // init PDF set
      const std::string setname = argv[7];
      const int setimem = atoi(argv[8]);
      const LHAPDF::PDF* basepdf = LHAPDF::mkPDF( setname, setimem);
      LHAPDF::setVerbosity(0); // default is 1;

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
      double gg_lo_result,gg_lo_error;

      integrators::result<double> resultPlusConst;
      integrators::result<double> resultPlusInt1;
      integrators::result<double> resultPlusInt2;
      integrators::result<double> resultRegNLO;
      integrators::result<double> result_gq_NLO;
      integrators::result<double> result_qqbar_NLO;
      double gg_nlo_result,gg_nlo_error;

      integrators::result<double> resultRegNNLO;
      integrators::result<double> result_gq_NNLO;
      integrators::result<double> result_qqbar_NNLO;
      integrators::result<double> result_qq_NNLO;
      integrators::result<double> result_q1q2_NNLO;
      double gg_nnlo_result,gg_nnlo_error;

      integrators::result<double> resultRegN3LO;
      integrators::result<double> result_gq_N3LO;
      integrators::result<double> result_qqbar_N3LO;
      integrators::result<double> result_qq_N3LO;
      integrators::result<double> result_q1q2_N3LO;
      double gg_n3lo_result,gg_n3lo_error;


      // Perform all the required QMC integrations

      if(qcdorder>=0)
	{
	  real_integrator.minn = lattice/10;
	  functor_delta.k = 0;
	  resultdelta = real_integrator.integrate(functor_delta);
	  gg_lo_result = resultdelta.integral;
	  gg_lo_error  = resultdelta.error;

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
	      gg_nlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNLO.integral;
	      gg_nlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
	      			    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNLO.error,2));
	      result_gq_NLO = real_integrator.integrate(functor_gq_NLO);
	      result_qqbar_NLO = real_integrator.integrate(functor_qqbar_NLO);

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
		  gg_nnlo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegNNLO.integral;
		  gg_nnlo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		  			 pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegNNLO.error,2));
		  result_gq_NNLO    = real_integrator.integrate(functor_gq_NNLO);
		  result_qqbar_NNLO = real_integrator.integrate(functor_qqbar_NNLO);
		  result_qq_NNLO    = real_integrator.integrate(functor_qq_NNLO);
		  result_q1q2_NNLO  = real_integrator.integrate(functor_q1q2_NNLO);

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
		      real_integrator.minn = 10*lattice;
		      gg_n3lo_result = resultdelta.integral + resultPlusConst.integral + resultPlusInt1.integral + resultPlusInt2.integral + resultRegN3LO.integral;
		      gg_n3lo_error  = sqrt(pow(resultdelta.error,2) + pow(resultPlusConst.error,2) +
		      			    pow(resultPlusInt1.error,2) + pow(resultPlusInt2.error,2) + pow(resultRegN3LO.error,2));
		      result_gq_N3LO   = real_integrator.integrate(functor_gq_N3LO);
		      real_integrator.minn = lattice;
		      result_qqbar_N3LO = real_integrator.integrate(functor_qqbar_N3LO);
		      result_qq_N3LO    = real_integrator.integrate(functor_qq_N3LO);
		      result_q1q2_N3LO  = real_integrator.integrate(functor_q1q2_N3LO);
		    }
		}
	    }
	}

      // Building the cross section: adding back alphaS, evolve to muR, add LO normalization
      // Result of this code: xs(g g -> H + X) in pb

      double BornggH0, BornggH;

      BornggH0 = constants::gevtopb*constants::Pi/
	(72.0*(constants::Nc*constants::Nc-1)*constants::vev*constants::vev);

      BornggH0 = BornggH0*oneloopfac(constants::Mt); // rescale by exact LO

      double muf  = xmuf*scale0;
      double muf2 = muf*muf;
      double xmur;
      double mur;
      double mur2;

      // wilson = 1 + as*wilson1 + as^2*wilson2 + as^3*wilson3
      // from arXiv:1607.05548 both for OS and MSbar schemes
      double logt  = log(constants::Mt*constants::Mt/muf2);
      double logt2 = logt*logt;

      const double wilson1 = 11.0/4.0;
      double wilson2 = (2777 - 201*constants::nf - 6*logt*(57 + 16*constants::nf))/288.0;
      // MSbar:
      // double wilson3 = (-432*logt2*(-33 + 2*constants::nf)*(57 + 16*constants::nf) - 
      // 			144*logt*(10398 + 11*constants::nf*(160 + 7*constants::nf)) - 
      // 			2*constants::nf*(-241746 + 27460*constants::nf + 997011*constants::Zeta3) + 
      // 			3*(-5785318 + 8081487*constants::Zeta3))/248832.0;
      // OS:
      double wilson3 = (-432*logt2*(-33 + 2*constants::nf)*(57 + 16*constants::nf) - 
			144*logt*(14502 + 7*constants::nf*(416 + 11*constants::nf)) - 
			2*constants::nf*(-352338 + 27460*constants::nf + 997011*constants::Zeta3) + 
			3*(-5522662 + 8081487*constants::Zeta3))/248832.0;

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

      double eta1, eta2, eta3;
      double deta12, deta22, deta32;

      std::string finalfile;
      std::stringstream filename;
      std::string header;

      int imax;
      double dxmur;

      double asopimz = (basepdf->alphasQ(constants::MZ))/constants::Pi;

      if(argc>=10 && std::strcmp(argv[9],"--scale") == 0)
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
	  filename << "ggH_xs_pp_" << energy << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	  header = "# Inclusive cross section for SM Higgs production in gluon fusion xs(g g -> H) in the heavy-top limit, p-p collider";
	}
      else
	{
	  filename << "ggH_xs_ppbar_" << energy << "tev_pdf" << setimem << "_muf" << xmuf << ".txt";
	  header = "# Inclusive cross section for SM Higgs production in gluon fusion xs(g g -> H) in the heavy-top limit, p-pbar collider";
	}

      
      filename >> finalfile;
      std::ofstream fa(finalfile);

      switch(qcdorder)
	{
	case 0:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central factorization and renormalization scales muF0 = muR0 = MH/2\n# mu_R/muR0\t" << "mu_F/muF0\t" << "xs_LO (pb)" << std::endl;
	  break;
	case 1:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central factorization and renormalization scales muF0 = muR0 = MH/2\n# mu_R/muR0\t" << "mu_F/muF0\t" << "xs_LO (pb)\t" << "xs_NLO (pb)" << std::endl;
	  break;
	case 2:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central factorization and renormalization scales muF0 = muR0 = MH/2\n# mu_R/muR0\t" << "mu_F/muF0\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t"
	     << "xs_NNLO (pb)" << std::endl;
	  break;
	case 3:
	  fa << header << std::fixed << std::setprecision(2) << ", sqrt(S) = " << energy
	     << " TeV, central factorization and renormalization scales muF0 = muR0 = MH/2\n# mu_R/muR0\t" << "mu_F/muF0\t" << "xs_LO (pb)\t" << "xs_NLO (pb)\t"
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

	  // alphaS(mur) at LO
	  asopi   = as_n3loxs(mur, 0, asopimz);
	  asopi2  = asopi*asopi;

	  BornggH = BornggH0*asopi2;

      	  xslo_result = BornggH*gg_lo_result;
      	  xslo_error  = BornggH*gg_lo_error;
	  if(qcdorder==0)
	    {
	      fa << std::fixed << std::setprecision(3) << (mur/scale0) << "\t" << (muf/scale0) << "\t"
		 << std::setprecision(18) << xslo_result << "\t" << std::scientific << xslo_error << std::endl;
	    }
      	}
      if(qcdorder>=1)
      	{
	  // alphaS(mur) at NLO
	  asopi   = as_n3loxs(mur, 1, asopimz);
	  asopi2  = asopi*asopi;

	  BornggH = BornggH0*asopi2;

      	  logmu1 = log(mur2/muf2);

      	  xsnlo_result = BornggH*
	    (gg_lo_result*(1 + 2*asopi*wilson1) +
	     asopi*(gg_nlo_result + result_gq_NLO.integral + result_qqbar_NLO.integral +
		    2*constants::b0*logmu1*gg_lo_result)
	     );
      	  xsnlo_error  = BornggH*
	    sqrt(pow((1 + 2*asopi*wilson1)*gg_lo_error,2) +
		 asopi2*(pow(gg_nlo_error,2) + pow(result_gq_NLO.error,2) + pow(result_qqbar_NLO.error,2) +
			 4*pow(constants::b0*logmu1*gg_lo_error,2))
		 );

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
	  asopi   = as_n3loxs(mur, 2, asopimz);
	  asopi2  = asopi*asopi;

	  BornggH = BornggH0*asopi2;

      	  asopi4 = asopi2*asopi2;
	  logmu2 = logmu1*logmu1;

      	  xsnnlo_result = BornggH*
      	    (gg_lo_result*(1 + 2*asopi*wilson1 + asopi2*(2*wilson2 + wilson1*wilson1)) +
      	     asopi*((gg_nlo_result + result_gq_NLO.integral + result_qqbar_NLO.integral)*(1 + 2*asopi*wilson1) +
		    2*constants::b0*logmu1*gg_lo_result) +
      	     asopi2*(gg_nnlo_result + result_gq_NNLO.integral + result_qqbar_NNLO.integral +
		     result_qq_NNLO.integral + result_q1q2_NNLO.integral + 
		     (3*constants::b0*constants::b0*logmu2 + 2*constants::b1*logmu1)*gg_lo_result +
		     3*constants::b0*logmu1*(gg_nlo_result + result_gq_NLO.integral + result_qqbar_NLO.integral + 2*wilson1*gg_lo_result)
		     )
      	     );
      	  xsnnlo_error  = BornggH*
      	    sqrt(pow((1 + 2*asopi*wilson1 + asopi2*(2*wilson2 + wilson1*wilson1))*gg_lo_error,2) +
      		 asopi2*(pow(1+2*asopi*wilson1,2)*(pow(gg_nlo_error,2) + pow(result_gq_NLO.error,2) + pow(result_qqbar_NLO.error,2)) +
			 4*pow(constants::b0*logmu1*gg_lo_error,2)) +
      		 asopi4*(pow(gg_nnlo_error,2) + pow(result_gq_NNLO.error,2) + pow(result_qqbar_NNLO.error,2) +
			 pow(result_qq_NNLO.error,2) + pow(result_q1q2_NNLO.error,2) + 
			 pow((3*constants::b0*constants::b0*logmu2 + 2*constants::b1*logmu1)*gg_lo_error,2) +
			 9*constants::b0*constants::b0*logmu2*
			 (pow(gg_nlo_error,2) + pow(result_gq_NLO.error,2) +
			  pow(result_qqbar_NLO.error,2) + 4*wilson1*wilson1*pow(gg_lo_error,2))
			 )
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
	  asopi   = as_n3loxs(mur, 3, asopimz);
	  asopi2 = asopi*asopi;

	  BornggH = BornggH0*asopi2;

      	  asopi3 = asopi*asopi2;
	  asopi4 = asopi2*asopi2;
      	  asopi6 = asopi3*asopi3;
	  logmu3 = logmu1*logmu2;

	  eta1 = gg_nlo_result + result_gq_NLO.integral + result_qqbar_NLO.integral + 2*wilson1*gg_lo_result;
	  deta12 = pow(gg_nlo_error, 2) + pow(result_gq_NLO.error, 2) + pow(result_qqbar_NLO.error, 2) +
	    4*pow(wilson1*gg_lo_error, 2);
	  eta2 = gg_nnlo_result + result_gq_NNLO.integral + result_qqbar_NNLO.integral +
	    result_qq_NNLO.integral + result_q1q2_NNLO.integral + 2*wilson1*
	    (gg_nlo_result + result_gq_NLO.integral + result_qqbar_NLO.integral) +
	    (wilson1*wilson1 + 2*wilson2)*gg_lo_result;
	  deta22 = pow(gg_nnlo_error, 2) + pow(result_gq_NNLO.error, 2) + pow(result_qqbar_NNLO.error, 2) + 
	    pow(result_qq_NNLO.error, 2) + pow(result_q1q2_NNLO.error, 2) +
	    4*wilson1*wilson1*(pow(gg_nlo_error, 2) + pow(result_gq_NLO.error, 2) +
			       pow(result_qqbar_NLO.error, 2)) +
	    pow((wilson1*wilson1 + 2*wilson2)*gg_lo_error, 2);
	  eta3 = gg_n3lo_result + result_gq_N3LO.integral + result_qqbar_N3LO.integral +
	    result_qq_N3LO.integral + result_q1q2_N3LO.integral +
	    2*(wilson1*wilson2 + wilson3)*gg_lo_result +
	    (wilson1*wilson1 + 2*wilson2)*
	    (gg_nlo_result + result_gq_NLO.integral + result_qqbar_NLO.integral) +
	    2*wilson1*(gg_nnlo_result + result_gq_NNLO.integral + result_qqbar_NNLO.integral +
		       result_qq_NNLO.integral + result_q1q2_NNLO.integral);
	  deta32 = pow(gg_n3lo_error, 2) + pow(result_gq_N3LO.error, 2) + pow(result_qqbar_N3LO.error, 2) + 
	    pow(result_qq_N3LO.error, 2) + pow(result_q1q2_N3LO.error, 2) +
	    4*pow((wilson1*wilson2 + wilson3)*gg_lo_error, 2) +
	    pow(wilson1*wilson1 + 2*wilson2, 2)*
	    (pow(gg_nlo_error, 2) + pow(result_gq_NLO.error, 2) + pow(result_qqbar_NLO.error, 2)) +
	    4*wilson1*wilson1*
	    (pow(gg_nnlo_error, 2) + pow(result_gq_NNLO.error, 2) + pow(result_qqbar_NNLO.error, 2) + 
	     pow(result_qq_NNLO.error, 2) + pow(result_q1q2_NNLO.error, 2));

	  xsn3lo_result = BornggH*
	    (gg_lo_result +
	     asopi*(eta1 + 2*constants::b0*logmu1*gg_lo_result) +
	     asopi2*(eta2 + 3*constants::b0*logmu1*eta1 +
	  	     (3*constants::b0*constants::b0*logmu2 + 2*constants::b1*logmu1)*gg_lo_result) +
	     asopi3*(eta3 + 4*constants::b0*logmu1*eta2 +
	  	     3*(constants::b1*logmu1 + 2*constants::b0*constants::b0*logmu2)*eta1 +
	  	     (2*constants::b2*logmu1 + 7*constants::b0*constants::b1*logmu2 +
	  	      4*constants::b0*constants::b0*constants::b0*logmu3)*gg_lo_result)
	     );
	  xsn3lo_error = BornggH*
	    sqrt(pow(gg_lo_error, 2) +
	  	 asopi2*(deta12 + 4*logmu2*pow(constants::b0*gg_lo_error, 2)) +
	  	 asopi4*(deta22 + 9*logmu2*constants::b0*constants::b0*deta12 +
	  		 pow((3*constants::b0*constants::b0*logmu2 + 2*constants::b1*logmu1)*gg_lo_error,2)) +
	  	 asopi6*(deta32 + 16*constants::b0*constants::b0*logmu2*deta22 +
	  		 9*pow(constants::b1*logmu1 + 2*constants::b0*constants::b0*logmu2,2)*deta12 +
	  		 pow((2*constants::b2*logmu1 + 7*constants::b0*constants::b1*logmu2 +
	  		      4*constants::b0*constants::b0*constants::b0*logmu3)*gg_lo_error,2))
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
