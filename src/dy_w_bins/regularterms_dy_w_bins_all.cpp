/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 07/05/2021
Regular hard terms for all DY W subprocesses up to N3LO QCD including binning
*********************************************************************
********************************************************************* */

#include "dy_w_kernels.h"

// pdf functions
#include "pdffunctions_w.h"

#include "dy_functions_w_bins.h"

#include "constants.h"

// Header for the routines alphaS(muR)
#include "alphaS.h"

static const double eps = 1.e-10;

static const double GammaW = 2.085;


//////////////////////////////////////////////////
///////////////// d-ubar channel /////////////////
//////////////////////////////////////////////////

// NLO d-ubar regular term
double dub_regular_nlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 1, asopimz);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi*dub_regular_kernel_nlo(x1, log1);

  res = fac*res;
  res = res*dlumdub(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}


// NNLO d-ubar regular term
double dub_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2, res0;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res0 = (asopi + asopi2*constants::b0*logmu1)*dub_regular_kernel_nlo(x1, log1);
  
  auto result0 = dub_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res1 = asopi2*res1;
  res2 = std::get<1>(result0);
  res2 = asopi2*res2;
  
  res = fac*((res0 + res1)*dlumdub(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumdub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}

// N3LO d-ubar regular term
double dub_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1,logmu2;
  double res1, res2, res0, res01, res02;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res0 = (asopi +
	  asopi2*constants::b0*logmu1 +
	  asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*dub_regular_kernel_nlo(x1, log1);
  
  auto result0 = dub_regular_kernel_nnlo(x1, log1);
  res01 = std::get<0>(result0);
  res01 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res01;
  res02 = std::get<1>(result0);
  res02 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res02;

  auto result1 = dub_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result1);
  res1 = asopi3*res1;
  res2 = std::get<1>(result1);
  res2 = asopi3*res2;

  res = fac*((res0 + res01 + res1)*dlumdub(x2,tau/x1/x2,muf2,pdf) +
	     (res02 + res2)*dlumdub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
///////////////// g-ubar channel /////////////////
//////////////////////////////////////////////////

// NLO g-ub + d-g regular term
double gub_regular_nlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 1, asopimz);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi*gub_regular_kernel_nlo(x1, log1);

  res = fac*res;
  res = res*dlumgub(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}

// NNLO g-ub + d-g regular term
double gub_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double res1, res2;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi + asopi2*constants::b0*logmu1)*gub_regular_kernel_nlo(x1, log1);
  res2 = asopi2*gub_regular_kernel_nnlo(x1, log1);

  res = fac*(res1 + res2)*dlumgub(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;
  
  return res;
}

// N3LO g-ub + d-g regular term
double gub_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1,logmu2;
  double res1, res2, res00, res01;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res00 = (asopi +
	   asopi2*constants::b0*logmu1 +
	   asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*gub_regular_kernel_nlo(x1, log1);
  res01 = (asopi2 + 2*asopi3*constants::b0*logmu1)*gub_regular_kernel_nnlo(x1, log1);

  auto result0 = gub_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res1 = asopi3*res1;
  res2 = std::get<1>(result0);
  res2 = asopi3*res2;

  res = fac*((res00 + res01 + res1)*dlumgub(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumgub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
////////////////// g-g channel ///////////////////
//////////////////////////////////////////////////

// NNLO g-g regular term
double gg_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*gg_regular_kernel_nnlo(x1, log1);
  
  res = fac*res*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}

// N3LO g-g regular term
double gg_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = (asopi2 + 2*asopi3*constants::b0*logmu1)*gg_regular_kernel_nnlo(x1, log1) +
    asopi3*gg_regular_kernel_n3lo(x1, log1);
  
  res = fac*res*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
//////////////// g-dbar channel //////////////////
//////////////////////////////////////////////////

// N3LO g-dbar + dbar-g regular term
double gdb_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1, res2;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = gdb_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res1 = asopi3*res1;
  res2 = std::get<1>(result0);
  res2 = asopi3*res2;
  
  res = fac*(res1*dlumgu(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumgu2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
//////////////// c-ubar channel //////////////////
//////////////////////////////////////////////////

// NNLO c-ubar + u-cbar regular term
double cub_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*cub_regular_kernel_nnlo(x1, log1);

  res = fac*res*(dlumcub(x2,tau/x1/x2,muf2,pdf)+dlumcub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}

// N3LO c-ubar + u-cbar regular term
double cub_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = (asopi2 + 2*asopi3*constants::b0*logmu1)*cub_regular_kernel_nnlo(x1, log1) +
    asopi3*cub_regular_kernel_n3lo(x1, log1);
  
  res = fac*res*(dlumcub(x2,tau/x1/x2,muf2,pdf)+dlumcub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
//////////////// q-qbar channel //////////////////
//////////////////////////////////////////////////

// NNLO q-qbar regular term
double qqb_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1, res2;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res1 = asopi2*res1;
  res2 = std::get<1>(result0);
  res2 = asopi2*res2;

  res = fac*(res1*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqqb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}

// N3LO q-qbar regular term
double qqb_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double res1, res2, res01, res02;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  res01 = std::get<0>(result0);
  res01 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res01;
  res02 = std::get<1>(result0);
  res02 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res02;

  auto result1 = qqb_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result1);
  res1 = asopi3*res1;
  res2 = std::get<1>(result1);
  res2 = asopi3*res2;

  res = fac*((res01 + res1)*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	     (res02 + res2)*dlumqqb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
//////////// q-q / qbar-qbar channel /////////////
//////////////////////////////////////////////////

// NNLO q-q + qbar-qbar regular term
double qq_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*qq_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}

// N3LO q-q regular term
double qq_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = (asopi2 + 2*asopi3*constants::b0*logmu1)*qq_regular_kernel_nnlo(x1, log1) +
    asopi3*qq_regular_kernel_n3lo(x1, log1);

  res = fac*res*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
////////////////// u-d channel ///////////////////
//////////////////////////////////////////////////

// NNLO u-d regular term
double qqprime_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1, res2;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqprime_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res1 = asopi2*res1;
  res2 = std::get<1>(result0);
  res2 = asopi2*res2;

  res = fac*(res1*dlumqqprime(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqqprime2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}

// N3LO u-d regular term
double qqprime_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double res1, res2, res01, res02;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqprime_regular_kernel_nnlo(x1, log1);
  res01 = std::get<0>(result0);
  res01 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res01;
  res02 = std::get<1>(result0);
  res02 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res02;

  auto result1 = qqprime_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result1);
  res1 = asopi3*res1;
  res2 = std::get<1>(result1);
  res2 = asopi3*res2;

  res = fac*((res01 + res1)*dlumqqprime(x2,tau/x1/x2,muf2,pdf) +
	     (res02 + res2)*dlumqqprime2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
/////////////// ubar-dbar channel ////////////////
//////////////////////////////////////////////////

// NNLO ubar-dbar regular term
double qbqprimeb_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1, res2;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qbqprimeb_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res1 = asopi2*res1;
  res2 = std::get<1>(result0);
  res2 = asopi2*res2;

  res = fac*(res1*dlumqbqprimeb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqbqprimeb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}

// N3LO ubar-dbar regular term
double qbqprimeb_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double res1, res2, res01, res02;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qbqprimeb_regular_kernel_nnlo(x1, log1);
  res01 = std::get<0>(result0);
  res01 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res01;
  res02 = std::get<1>(result0);
  res02 = (asopi2 + 2*asopi3*constants::b0*logmu1)*res02;

  auto result1 = qbqprimeb_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result1);
  res1 = asopi3*res1;
  res2 = std::get<1>(result1);
  res2 = asopi3*res2;

  res = fac*((res01 + res1)*dlumqbqprimeb(x2,tau/x1/x2,muf2,pdf) +
	     (res02 + res2)*dlumqbqprimeb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
////////////////// d-s channel ///////////////////
//////////////////////////////////////////////////

// NNLO d-s regular term
double ds_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*ds_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumds(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}

// N3LO d-s regular term
double ds_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = (asopi2 + 2*asopi3*constants::b0*logmu1)*ds_regular_kernel_nnlo(x1, log1) +
    asopi3*ds_regular_kernel_n3lo(x1, log1);

  res = fac*res*dlumds(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}


//////////////////////////////////////////////////
/////////////// ubar-cbar channel ////////////////
//////////////////////////////////////////////////

// NNLO ub-cb regular term
double ubcb_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  asopi = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*ubcb_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumubcb(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}

// N3LO ub-cb regular term
double ubcb_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, BW;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  BW = (Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW) +
    constants::MW*constants::MW*GammaW*GammaW;
  BW = tau/BW;

  if(muf0==-1)
    {
      muf2 = Q2*xmuf*xmuf;
    }
  else
    {
      muf2 = muf0*muf0*xmuf*xmuf;
    }
  if(mur0==-1)
    {
      mur = sqrt(Q2)*xmur;
    }
  else
    {
      mur = mur0*xmur;
    }

  mur2 = mur*mur;

  asopi = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;
  logmu1 = log(mur2/muf2);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = (asopi2 + 2*asopi3*constants::b0*logmu1)*ubcb_regular_kernel_nnlo(x1, log1) +
    asopi3*ubcb_regular_kernel_n3lo(x1, log1);

  res = fac*res*dlumubcb(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*BW;

  return res;
}
