/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 07/10/2020
Soft+Virtual contributions for the DY process D Ubar -> W- -> l- ~nu_l up to N3LO QCD
*********************************************************************
********************************************************************* */

// pdf functions
#include "pdffunctions_w.h"

#include "dy_functions_wh.h"

#include "constants.h"

static const double eps = 1.e-8;

//static const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

// Virtual delta(z) contribution up to N3LO
double delta(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, lambda, gammlo, Q2;
  double delterms;
  double z;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  muf2 = muf*muf;

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[1])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  delterms = 0.0;
  
  switch(k)
    {
    case 0: // LO
      delterms = 1.0;
      break;
    case 1: // NLO
      log1 = log(Q2/muf2);
      delterms = 2.0/9.0*(-24.0 + 9*log1 + 2*constants::Pi2);
      break;
    case 2: // NNLO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      delterms = 1.0/6480.0*
	(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3));
      break;
    case 3: // N3LO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      delterms =
	(684831*constants::Pi4 - 187912*constants::Pi6 - 3240*log3*(33 + 32*constants::Pi2 - 4096*constants::Zeta3) + 
	 270*constants::Pi2*(47045 + 195068*constants::Zeta3) - 108*log2*
	 (36675 + 35200*constants::Pi2 + 1408*constants::Pi4 + 243000*constants::Zeta3) +
	 45*(-2117173 + 461868*constants::Zeta3 + 9038016*constants::Zeta3*constants::Zeta3 - 12599856*constants::Zeta5) - 
	 180*log1*(267*constants::Pi4 + 8*constants::Pi2*(-5773 + 47448*constants::Zeta3) - 
		   72*(2547 + 4373*constants::Zeta3 + 39654*constants::Zeta5)))/2.09952e6;
      break;
    }

  z   = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  fac = -(1.0-2.0*eps)*z*log(tau)*fac;
  
  res = fac*delterms;
  res = res*dlumdub(z,tau/z,muf2,pdf)/z;

  res = res*gammlo;

  return res;
}

////////////////////////////////////////////////////
// Soft PlusDistributions contributions up to N3LO
/*
     with f(x1) = Log[1-x1]^k/(1-x1),
     and  g(x1) = Integrate[pdfa(x2)/x2*pdfb(tau/x2/x1)/x1,{x2,tau/x1,1}:

     Plus = Integrate[(f(x1))_+*g(x1),{x1,tau,1}] 
          = Log[1-tau]^(k+1)/(k+1)*g(1)
            + Integrate[f(x1)*(g(x1)-g(1)),{x1,tau,1}]

          = PlusConst + PlusInt1 + PlusInt2

     with
     
     PlusConst = Log[1-tau]^(k+1)/(k+1)*g(1)
     PlusInt1  = Integrate[Log[1-x1]^k/(1-x1)*pdfa(x2)/x2*(pdfb(tau/x1/x2)/x1 - pdfb(tau/x2)),{x1,tau,1},{x2,tau/x1,1}]
     PlusInt2  = -Integrate[Log[1-x1]^k/(1-x1)*pdfa(x2)/x2*pdfb(tau/x2),{x1,tau,1},{x2,tau,tau/x1}]
*/


// PlusConst term
double PlusConst(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, lambda, gammlo, Q2;
  double plusterms[6];
  double x2;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  muf2 = muf*muf;

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[1])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;
  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);
  gammlo = (Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  x2   = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  fac = -(1.0-2.0*eps)*x2*log(tau)*fac;

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      log1 = log(Q2/muf2);
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]*log(1.0-tau) + plusterms[1]*pow(log(1.0-tau),2)/2.0;
      break;
    case 2: // NNLO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(41 + 25*constants::Pi2 - 3*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*pow(log(1.0-tau),2)/2.0 +
	plusterms[2]*pow(log(1.0-tau),3)/3.0 +
	plusterms[3]*pow(log(1.0-tau),4)/4.0;
      break;
    case 3: // N3LO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 1.0/87480.0*
	(8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	 540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	 120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	 18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	 5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5));
      plusterms[1] = -1.0/2430.0*
	(45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	 1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	 20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3));
      plusterms[2] = -2.0/81.0*
	(5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 + 756*log2 - 384*log3 - 19896*constants::Zeta3);
      plusterms[3] = -(4.0/81.0)*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2);
      plusterms[4] = 160.0/81.0*(-23 + 24*log1);
      plusterms[5] = 512.0/27.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*pow(log(1.0-tau),2)/2.0 +
	plusterms[2]*pow(log(1.0-tau),3)/3.0 +
	plusterms[3]*pow(log(1.0-tau),4)/4.0 +
	plusterms[4]*pow(log(1.0-tau),5)/5.0 +
	plusterms[5]*pow(log(1.0-tau),6)/6.0;
      break;
    }
  res = fac*res;
  res = res*dlumdub(x2,tau/x2,muf2,pdf)/x2;

  res = res*gammlo;

  return res;
}

// PlusInt1 term
double PlusInt1(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, lambda, gammlo, Q2;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  muf2 = muf*muf;

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

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      log1 = log(Q2/muf2);
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      //      std::cout << "DEBUG\t" << res << std::endl;
      break;
    case 2: // NNLO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(41 + 25*constants::Pi2 - 3*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*pow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*pow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 1.0/87480.0*
	(8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	 540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	 120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	 18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	 5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5));
      plusterms[1] = -1.0/2430.0*
	(45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	 1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	 20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3));
      plusterms[2] = -2.0/81.0*
	(5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 + 756*log2 - 384*log3 - 19896*constants::Zeta3);
      plusterms[3] = -(4.0/81.0)*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2);
      plusterms[4] = 160.0/81.0*(-23 + 24*log1);
      plusterms[5] = 512.0/27.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*pow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*pow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*pow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*pow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = fac*res;
  res = res*( dlumdub(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumdub(x2,tau/x2,muf2,pdf)/x2 );

  res = res*gammlo;

  return res;
}


// PlusInt2 term
double PlusInt2(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, lambda, gammlo, Q2;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  muf2 = muf*muf;

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
  x2 = tau*exp(-(eps+(1.0-2.0*eps)*X[1])*log(x1));
  fac = pow(1.0-2.0*eps,2)*x1*log(tau)*x2*log(x1)*fac;

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      log1 = log(Q2/muf2);
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(41 + 25*constants::Pi2 - 3*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*pow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*pow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 1.0/87480.0*
	(8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	 540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	 120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	 18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	 5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5));
      plusterms[1] = -1.0/2430.0*
	(45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	 1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	 20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3));
      plusterms[2] = -2.0/81.0*
	(5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 + 756*log2 - 384*log3 - 19896*constants::Zeta3);
      plusterms[3] = -(4.0/81.0)*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2);
      plusterms[4] = 160.0/81.0*(-23 + 24*log1);
      plusterms[5] = 512.0/27.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*pow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*pow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*pow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*pow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = -fac*res;
  res = res*dlumdub(x2,tau/x2,muf2,pdf)/x2;

  res = res*gammlo;

  return res;
}
