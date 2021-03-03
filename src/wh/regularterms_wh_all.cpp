/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 07/02/2021
Regular hard terms for all W H subprocesses (DY-type) up to N3LO QCD
*********************************************************************
********************************************************************* */

#include "dy_w_kernels.h"

// pdf functions
#include "pdffunctions_w.h"

#include "dy_functions_wh.h"

#include "constants.h"

static const double eps = 1.e-8;

static const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);


//////////////////////////////////////////////////
///////////////// d-ubar channel /////////////////
//////////////////////////////////////////////////

// NLO d-ubar regular term
double dub_regular_nlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = dub_regular_kernel_nlo(x1, log1);

  res = fac*res;
  res = res*dlumdub(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}


// NNLO d-ubar regular term
double dub_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1, res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = dub_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);
  
  res = fac*(res1*dlumdub(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumdub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO d-ubar regular term
double dub_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1, res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = dub_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumdub(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumdub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
///////////////// g-ubar channel /////////////////
//////////////////////////////////////////////////

// NLO g-ub + d-g regular term
double gub_regular_nlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = gub_regular_kernel_nlo(x1, log1);

  res = fac*res;
  res = res*dlumgub(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}

// NNLO g-ub + d-g regular term
double gub_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = gub_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumgub(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO g-ub + d-g regular term
double gub_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1, res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = gub_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumgub(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumgub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
////////////////// g-g channel ///////////////////
//////////////////////////////////////////////////

// NNLO g-g regular term
double gg_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = gg_regular_kernel_nnlo(x1, log1);
  
  res = fac*res*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO g-g regular term
double gg_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = gg_regular_kernel_n3lo(x1, log1);
  
  res = fac*res*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
//////////////// g-dbar channel //////////////////
//////////////////////////////////////////////////

// N3LO g-dbar + dbar-g regular term
double gdb_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1, res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = gdb_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);
  
  res = fac*(res1*dlumgu(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumgu2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
//////////////// c-ubar channel //////////////////
//////////////////////////////////////////////////

// NNLO c-ubar + u-cbar regular term
double cub_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = cub_regular_kernel_nnlo(x1, log1);

  res = fac*res*(dlumcub(x2,tau/x1/x2,muf2,pdf)+dlumcub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO c-ubar + u-cbar regular term
double cub_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = cub_regular_kernel_n3lo(x1, log1);
  
  res = fac*res*(dlumcub(x2,tau/x1/x2,muf2,pdf)+dlumcub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
//////////////// q-qbar channel //////////////////
//////////////////////////////////////////////////

// NNLO q-qbar regular term
double qqb_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1,res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqqb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO q-qbar regular term
double qqb_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1,res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqqb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
//////////// q-q / qbar-qbar channel /////////////
//////////////////////////////////////////////////

// NNLO q-q + qbar-qbar regular term
double qq_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = qq_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO q-q regular term
double qq_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = qq_regular_kernel_n3lo(x1, log1);

  res = fac*res*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
////////////////// u-d channel ///////////////////
//////////////////////////////////////////////////

// NNLO u-d regular term
double qqprime_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1,res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqprime_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumqqprime(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqqprime2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO u-d regular term
double qqprime_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1,res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqprime_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumqqprime(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqqprime2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
/////////////// ubar-dbar channel ////////////////
//////////////////////////////////////////////////

// NNLO ubar-dbar regular term
double qbqprimeb_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1,res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qbqprimeb_regular_kernel_nnlo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumqbqprimeb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqbqprimeb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO ubar-dbar regular term
double qbqprimeb_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;
  double res1,res2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qbqprimeb_regular_kernel_n3lo(x1, log1);
  res1 = std::get<0>(result0);
  res2 = std::get<1>(result0);

  res = fac*(res1*dlumqbqprimeb(x2,tau/x1/x2,muf2,pdf) +
	     res2*dlumqbqprimeb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
////////////////// d-s channel ///////////////////
//////////////////////////////////////////////////

// NNLO d-s regular term
double ds_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = ds_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumds(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO d-s regular term
double ds_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = ds_regular_kernel_n3lo(x1, log1);

  res = fac*res*dlumds(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}


//////////////////////////////////////////////////
/////////////// ubar-cbar channel ////////////////
//////////////////////////////////////////////////

// NNLO ub-cb regular term
double ubcb_regular_nnlo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1,log2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = ubcb_regular_kernel_nnlo(x1, log1);

  res = fac*res*dlumubcb(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}

// N3LO ub-cb regular term
double ubcb_regular_n3lo(const double X[], const double s, const double muf, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, lambda, gammlo, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2;
  double log1,log2;

  muf2 = pow(muf,2);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = ubcb_regular_kernel_n3lo(x1, log1);

  res = fac*res*dlumubcb(x2,tau/x1/x2,muf2,pdf)/x1/x2;
  res = res*gammlo;

  return res;
}
