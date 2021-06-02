/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 11/05/2021
Regular hard terms for all DY subprocesses up to N3LO QCD including binning
*********************************************************************
********************************************************************* */

#include "dy_kernels.h"

// pdf functions
#include "pdffunctions.h"

#include "dy_functions_bins.h"

#include "constants.h"

// Header for the routines alphaS(muR)
#include "alphaS.h"

static const double eps = 1.e-10;


//////////////////////////////////////////////////
///////////////// q-qbar channel /////////////////
//////////////////////////////////////////////////

// NLO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_nlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi*qqb_regular_kernel_nlo(x1, log1);

  res = fac*res;
  res = res*dlumqqb(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;

  return res;
}


// NNLO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double csum2 = 11.0/9.0;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res0, resuec2_nnlo, rescsum2_nnlo;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res0 = (asopi + asopi2*constants::b0*logmu1)*qqb_regular_kernel_nlo(x1, log1);
  
  auto result0 = qqb_regular_kernel_nnlo(x1, log1);

  resuec2_nnlo = std::get<0>(result0);
  resuec2_nnlo = asopi2*resuec2_nnlo;
  rescsum2_nnlo = std::get<1>(result0);
  rescsum2_nnlo = asopi2*rescsum2_nnlo;
  
  res = fac*((res0 + resuec2_nnlo)*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	     rescsum2_nnlo*csum2*(9.0/4.0*dlumuub(x2,tau/x1/x2,muf2,pdf) +
			     9.0*dlumddb(x2,tau/x1/x2,muf2,pdf)))/x1/x2;
  
  res = res/Den;

  return res;
}

// N3LO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double csum2 = 11.0/9.0;
  const double csumQuoverQu2 = 0.5;
  const double csumQdoverQd2 = -1.0;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1,logmu2;
  double res0;
  double resuec2_nnlo, rescsum2_nnlo;
  double resuec2_n3lo, rescsum2_n3lo, resueccsum_n3lo;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res0 = (asopi +
	  asopi2*constants::b0*logmu1 +
	  asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*qqb_regular_kernel_nlo(x1, log1);
  
  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  resuec2_nnlo = std::get<0>(result0);
  resuec2_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*resuec2_nnlo;
  rescsum2_nnlo = std::get<1>(result0);
  rescsum2_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*rescsum2_nnlo;

  auto result1 = qqb_regular_kernel_n3lo(x1, log1);
  resuec2_n3lo = std::get<0>(result1);
  resuec2_n3lo = asopi3*resuec2_n3lo;
  rescsum2_n3lo = std::get<1>(result1);
  rescsum2_n3lo = asopi3*rescsum2_n3lo;
  resueccsum_n3lo = std::get<2>(result1);
  resueccsum_n3lo = asopi3*resueccsum_n3lo;

  res = fac*((res0 + resuec2_nnlo + resuec2_n3lo)*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	     (rescsum2_nnlo + rescsum2_n3lo)*csum2*(9.0/4.0*dlumuub(x2,tau/x1/x2,muf2,pdf) +
						    9.0*dlumddb(x2,tau/x1/x2,muf2,pdf)) +
	     resueccsum_n3lo*(csumQuoverQu2*dlumuub(x2,tau/x1/x2,muf2,pdf) +
			      csumQdoverQd2*dlumddb(x2,tau/x1/x2,muf2,pdf)))/x1/x2;

  res = res/Den;

  return res;
}


//////////////////////////////////////////////////
/////////////////// g-q channel //////////////////
//////////////////////////////////////////////////

// NLO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_nlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi*gq_regular_kernel_nlo(x1, log1);

  res = fac*res;
  res = res*dlumgq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;

  return res;
}

// NNLO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double res1, res_nnlo;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi + asopi2*constants::b0*logmu1)*gq_regular_kernel_nlo(x1, log1);
  res_nnlo = asopi2*gq_regular_kernel_nnlo(x1, log1);

  res = fac*(res1 + res_nnlo)*dlumgq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;
  
  return res;
}

// N3LO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double csum2 = 11.0/9.0;
  const double csumQuoverQu2 = 0.5;
  const double csumQdoverQd2 = -1.0;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1,logmu2;
  double res1, res2, resuec2_n3lo, rescsum2_n3lo, resueccsum_n3lo;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi +
	  asopi2*constants::b0*logmu1 +
	  asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*gq_regular_kernel_nlo(x1, log1);
  res2 = (asopi2 + 2*asopi3*constants::b0*logmu1)*gq_regular_kernel_nnlo(x1, log1);

  auto result0 = gq_regular_kernel_n3lo(x1, log1);
  resuec2_n3lo = std::get<0>(result0);
  resuec2_n3lo = asopi3*resuec2_n3lo;
  rescsum2_n3lo = std::get<1>(result0);
  rescsum2_n3lo = asopi3*rescsum2_n3lo;
  resueccsum_n3lo = std::get<2>(result0);
  resueccsum_n3lo = asopi3*resueccsum_n3lo;

  res = fac*((res1 + res2 + resuec2_n3lo)*dlumgq(x2,tau/x1/x2,muf2,pdf) +
	     rescsum2_n3lo*csum2*(9.0/4.0*dlumgu(x2,tau/x1/x2,muf2,pdf) +
			     9.0*dlumgd(x2,tau/x1/x2,muf2,pdf)) +
	     resueccsum_n3lo*(csumQuoverQu2*dlumgu(x2,tau/x1/x2,muf2,pdf) +
			 csumQdoverQd2*dlumgd(x2,tau/x1/x2,muf2,pdf)))/x1/x2;

  res = res/Den;

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
  double tau, Q2, Den;
  double x1, x2;
  double csum2 = 11.0/9.0;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*gg_regular_kernel_nnlo(x1, log1);
  
  res = fac*res*csum2*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;

  return res;
}

// N3LO g-g regular term
double gg_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double csum2 = 11.0/9.0;
  double csumsq = 1.0/9.0;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double rescsum2,rescsumsq;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = gg_regular_kernel_n3lo(x1, log1);
  rescsum2 = std::get<0>(result0);
  rescsum2 = asopi3*rescsum2;
  rescsumsq = std::get<1>(result0);
  rescsumsq = asopi3*rescsumsq;

  res = (asopi2 + 2*asopi3*constants::b0*logmu1)*csum2*gg_regular_kernel_nnlo(x1, log1) +
    csum2*rescsum2 + csumsq*rescsumsq;
  
  res = fac*res*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;

  return res;
}


//////////////////////////////////////////////////
///////////////// q-q channel ////////////////////
//////////////////////////////////////////////////

// NNLO q-q regular term, electric charge stripped out and included in dlumqq
double qq_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*qq_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;

  return res;
}

// N3LO q-q regular term, electric charge stripped out and included in dlumqq
double qq_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*qq_regular_kernel_nnlo(x1, log1);
  res2 = asopi3*qq_regular_kernel_n3lo(x1, log1);
  
  res = fac*(res1 + res2)*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res/Den;

  return res;
}


//////////////////////////////////////////////////
////////////////// q-Q channel ///////////////////
//////////////////////////////////////////////////

// NNLO q-Q + u-d regular term, electric charge stripped out and included in dlumqQq
double qQq_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double res1, res2;
  double muf2, mur;
  double log1;
  double asopi, asopi2;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = asopi2*qQq_regular_kernel_nnlo(x1, log1);
  res2 = asopi2*ud_regular_kernel_nnlo(x1, log1);

  res = fac*(res1*dlumqQq(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res/Den;

  return res;
}

// N3LO q-Q + u-d regular term, electric charge stripped out and included in dlumqQq
double qQq_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double res1, res2;
  double res3, res4;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*qQq_regular_kernel_nnlo(x1, log1);    
  res2 = asopi3*qQq_regular_kernel_n3lo(x1, log1);
  res3 = (asopi2 + 2*asopi3*constants::b0*logmu1)*ud_regular_kernel_nnlo(x1, log1);    
  res4 = asopi3*ud_regular_kernel_n3lo(x1, log1);

  res = fac*(
	     (res1 + res2)*dlumqQq(x2,tau/x1/x2,muf2,pdf) +
	     (res3 + res4)*dlumud(x2,tau/x1/x2,muf2,pdf)
	     )/x1/x2;

  res = res/Den;

  return res;
}


//////////////////////////////////////////////////
//////////////// q-Qbar channel //////////////////
//////////////////////////////////////////////////

// NNLO q-Qbar + u-dbar regular term, electric charge stripped out and included in dlumqQqb
double qQqb_regular_nnlo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
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

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = asopi2*qQqb_regular_kernel_nnlo(x1, log1);
  res2 = asopi2*udb_regular_kernel_nnlo(x1, log1);

  res = fac*(res1*dlumqQqb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res/Den;

  return res;
}

// N3LO q-Qbar + u-dbar regular term, electric charge stripped out and included in dlumqQqb
double qQqb_regular_n3lo(const double X[], const double s,
		       const double muf0, const double xmuf, const double mur0, const double xmur,
		       const double q2min, const double q2max, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, Q2, Den;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double res1, res2;
  double res3, res4;
  double logmu1;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[2];
  fac = (q2max*(1.0-eps)-q2min*(1.0+eps));
  tau = Q2/s;

  Den = s*Q2;

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
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*qQqb_regular_kernel_nnlo(x1, log1);    
  res2 = asopi3*qQqb_regular_kernel_n3lo(x1, log1);
  res3 = (asopi2 + 2*asopi3*constants::b0*logmu1)*udb_regular_kernel_nnlo(x1, log1);    
  res4 = asopi3*udb_regular_kernel_n3lo(x1, log1);

  res = fac*(
	     (res1 + res2)*dlumqQqb(x2,tau/x1/x2,muf2,pdf) +
	     (res3 + res4)*dlumudb(x2,tau/x1/x2,muf2,pdf)
	     )/x1/x2;

  res = res/Den;

  return res;
}
