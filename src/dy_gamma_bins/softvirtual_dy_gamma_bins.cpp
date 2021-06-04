/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 11/05/2021
Soft+Virtual contributions for the DY process q qb -> gamma* -> l l up to N3LO QCD, including binning
*********************************************************************
********************************************************************* */

// pdf functions
#include "pdffunctions.h"

#include "dy_functions_bins.h"

#include "constants.h"

// Header for the routines alphaS(muR)
#include "alphaS.h"


double intpow(const double& x,int m){
        double res=1.0;
        for (int i=0;i<m;i++){
            res *= x;
        }
        return res;
    }

static const double eps = 1.e-10;

// Virtual delta(z) contribution up to N3LO
double delta(const double X[], const double s,
	     const double muf0, const double xmuf, const double mur0, const double xmur,
	     const double q2min, const double q2max, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, Q2, Den;
  double delterms;
  double udeltermsn3lo,ddeltermsn3lo;
  double z;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1,log2,log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[1];
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
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;

  delterms = 0.0;
  
  switch(k)
    {
    case 0: // LO
      delterms = 1.0;
      break;
    case 1: // NLO
      asopi = as_n3loxs(mur, 1, asopimz);
      log1 = log(Q2/muf2);
      delterms = 1.0 + 2*asopi/9.0*(-24.0 + 2*constants::Pi2 + 9*log1);
      break;
    case 2: // NNLO
      asopi  = as_n3loxs(mur, 2, asopimz);
      asopi2 = asopi*asopi;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      delterms = 1.0 +
	2*asopi/9.0*(-24.0 + 2*constants::Pi2 + 9*log1) +
	asopi2/6480.0*
	(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3)) +
	constants::b0*logmu1*2*asopi2/9.0*(-24.0 + 2*constants::Pi2 + 9*log1);
      break;
    case 3: // N3LO
      asopi  = as_n3loxs(mur, 3, asopimz);
      asopi2 = asopi*asopi;
      asopi3 = asopi*asopi2;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      delterms = 1.0 +
	2*asopi/9.0*(-24.0 + 2*constants::Pi2 + 9*log1) +
	asopi2/6480.0*
	(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3)) +
	constants::b0*logmu1*2*asopi2/9.0*(-24.0 + 2*constants::Pi2 + 9*log1) +
	asopi3*(684831*constants::Pi4 - 187912*constants::Pi6 - 3240*log3*(33 + 32*constants::Pi2 - 4096*constants::Zeta3) + 
		270*constants::Pi2*(47045 + 195068*constants::Zeta3) - 108*log2*
		(36675 + 35200*constants::Pi2 + 1408*constants::Pi4 + 243000*constants::Zeta3) + 
		45*(-2117173 + 461868*constants::Zeta3 + 9038016*constants::Zeta3*constants::Zeta3 - 12599856*constants::Zeta5) - 
		180*log1*(267*constants::Pi4 + 8*constants::Pi2*(-5773 + 47448*constants::Zeta3) - 
			  72*(2547 + 4373*constants::Zeta3 + 39654*constants::Zeta5)))/2.09952e6 +
	asopi3*(2*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)/9.0*(-24.0 + 2*constants::Pi2 + 9*log1) +
		2*constants::b0*logmu1*
		(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
		 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
		 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3))/6480.0);
      // For csum*uec/uec^2 terms
      udeltermsn3lo = 0.5*asopi3*(150*constants::Pi2 - constants::Pi4 + 60*(6 + 7*constants::Zeta3 - 40*constants::Zeta5))/1296.0;
      ddeltermsn3lo = -asopi3*(150*constants::Pi2 - constants::Pi4 + 60*(6 + 7*constants::Zeta3 - 40*constants::Zeta5))/1296.0;
      break;
    }

  z   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*z*log(tau)*fac;

  res = fac*delterms;
  res = res*dlumqqb(z,tau/z,muf2,pdf)/z;

  if(k==3)
    {
      res = res + fac/z*(udeltermsn3lo*dlumuub(z,tau/z,muf2,pdf)+ddeltermsn3lo*dlumddb(z,tau/z,muf2,pdf));
    }
  res = res/Den;
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


// PlusConst term, electric charge stripped out and included in dlumqqb
double PlusConst(const double X[], const double s,
		 const double muf0, const double xmuf, const double mur0, const double xmur,
		 const double q2min, const double q2max, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, Q2, Den;
  double plusterms[6];
  double x2;
  double fac;
  double res;
  double muf2,mur,mur2;
  double log1,log2,log3;
  double logmu1,logmu2;
  double asopi,asopi2,asopi3;

  Q2 = q2min*(1.0+eps)+(q2max*(1.0-eps)-q2min*(1.0+eps))*X[1];
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
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;

  x2   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*x2*log(tau)*fac;

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      asopi = as_n3loxs(mur, 1, asopimz);
      log1 = log(Q2/muf2);
      plusterms[0] = 8*asopi/3.0*log1;
      plusterms[1] = 16*asopi/3.0;

      res = plusterms[0]*log(1.0-tau) + plusterms[1]*intpow(log(1.0-tau),2)/2.0;
      break;
    case 2: // NNLO
      asopi = as_n3loxs(mur, 2, asopimz);
      asopi2 = asopi*asopi;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      plusterms[0] = 8*asopi/3.0*log1 +
	asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3) +
	asopi2*constants::b0*logmu1*8.0/3.0*log1;
      plusterms[1] = 16*asopi/3.0 -
	4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) +
	16*asopi2*constants::b0*logmu1/3.0;
      plusterms[2] = 4*asopi2/9.0*(-23 + 48*log1);
      plusterms[3] = 128*asopi2/9.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0;
      break;
    case 3: // N3LO
      asopi = as_n3loxs(mur, 3, asopimz);
      asopi2 = asopi*asopi;
      asopi3 = asopi*asopi2;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 8*asopi/3.0*log1 +
	asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3) +
	asopi2*constants::b0*logmu1*8.0/3.0*log1 +
	asopi3/87480.0*
	(8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	 540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	 120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	 18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	 5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5)) +
	asopi3*(8*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)*log1/3.0 +
		2*constants::b0*logmu1*
		(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3)/81.0);
      plusterms[1] = 16*asopi/3.0 -
	4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) +
	16*asopi2*constants::b0*logmu1/3.0 -
	asopi3/2430.0*
	(45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	 1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	 20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3)) +
	asopi3*(16*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)/3.0 -
		8*constants::b0*logmu1*
		(41 + 25*constants::Pi2 - 3*log1 - 48*log2)/27.0);
      plusterms[2] = 4*asopi2/9.0*(-23 + 48*log1) -
	2*asopi3/81.0*
	(5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 + 756*log2 - 384*log3 - 19896*constants::Zeta3) +
	8*asopi3*constants::b0*logmu1*(-23 + 48*log1)/9.0;
      plusterms[3] = 128*asopi2/9.0 -
	4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2) + 256*asopi3*constants::b0*logmu1/9.0;;
      plusterms[4] = 160*asopi3/81.0*(-23 + 24*log1);
      plusterms[5] = 512*asopi3/27.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0 +
	plusterms[4]*intpow(log(1.0-tau),5)/5.0 +
	plusterms[5]*intpow(log(1.0-tau),6)/6.0;
      break;
    }
  res = fac*res;
  res = res*dlumqqb(x2,tau/x2,muf2,pdf)/x2;

  res = res/Den;
  return res;
}

// PlusInt1 term, electric charge stripped out and included in dlumqqb
double PlusInt1(const double X[], const double s,
		const double muf0, const double xmuf, const double mur0, const double xmur,
		const double q2min, const double q2max, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, Q2, Den;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2,mur,mur2;
  double log1,log2,log3;
  double logmu1,logmu2;
  double asopi,asopi2,asopi3;

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
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      asopi = as_n3loxs(mur, 1, asopimz);
      log1 = log(Q2/muf2);
      plusterms[0] = 8*asopi/3.0*log1;
      plusterms[1] = 16*asopi/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      asopi = as_n3loxs(mur, 2, asopimz);
      asopi2 = asopi*asopi;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      plusterms[0] = 8*asopi/3.0*log1 +
	asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3) +
	asopi2*constants::b0*logmu1*8.0/3.0*log1;
      plusterms[1] =  16*asopi/3.0 -
	4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) +
	16*asopi2*constants::b0*logmu1/3.0;
      plusterms[2] = 4*asopi2/9.0*(-23 + 48*log1);
      plusterms[3] = 128*asopi2/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      asopi = as_n3loxs(mur, 3, asopimz);
      asopi2 = asopi*asopi;
      asopi3 = asopi*asopi2;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 8*asopi/3.0*log1 +
	asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3) +
	asopi2*constants::b0*logmu1*8.0/3.0*log1 +
	asopi3/87480.0*
	(8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	 540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	 120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	 18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	 5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5)) +
	asopi3*(8*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)*log1/3.0 +
		2*constants::b0*logmu1*
		(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3)/81.0);
      plusterms[1] = 16*asopi/3.0 -
	4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) +
	16*asopi2*constants::b0*logmu1/3.0 -
	asopi3/2430.0*
	(45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	 1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	 20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3)) +
	asopi3*(16*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)/3.0 -
		8*constants::b0*logmu1*
		(41 + 25*constants::Pi2 - 3*log1 - 48*log2)/27.0);
      plusterms[2] = 4*asopi2/9.0*(-23 + 48*log1) -
	2*asopi3/81.0*
	(5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 + 756*log2 - 384*log3 - 19896*constants::Zeta3) +
	8*asopi3*constants::b0*logmu1*(-23 + 48*log1)/9.0;
      plusterms[3] = 128*asopi2/9.0 -
	4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2) + 256*asopi3*constants::b0*logmu1/9.0;;
      plusterms[4] = 160*asopi3/81.0*(-23 + 24*log1);
      plusterms[5] = 512*asopi3/27.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = fac*res;
  res = res*( dlumqqb(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumqqb(x2,tau/x2,muf2,pdf)/x2 );

  res = res/Den;
  return res;
}


// PlusInt2 term, electric charge stripped out and included in dlumqqb
double PlusInt2(const double X[], const double s,
		const double muf0, const double xmuf, const double mur0, const double xmur,
		const double q2min, const double q2max, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, Q2, Den;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2,mur,mur2;
  double log1,log2,log3;
  double logmu1,logmu2;
  double asopi,asopi2,asopi3;

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
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau*exp(-(eps+(1.0-2.0*eps)*X[1])*log(x1));
  fac = intpow(1.0-2.0*eps,2)*x1*log(tau)*x2*log(x1)*fac;

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      asopi = as_n3loxs(mur, 1, asopimz);
      log1 = log(Q2/muf2);
      plusterms[0] = 8*asopi/3.0*log1;
      plusterms[1] = 16*asopi/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      asopi = as_n3loxs(mur, 2, asopimz);
      asopi2 = asopi*asopi;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      plusterms[0] = 8*asopi/3.0*log1 +
	asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3) +
	asopi2*constants::b0*logmu1*8.0/3.0*log1;
      plusterms[1] = 16*asopi/3.0 -
	4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) +
	16*asopi2*constants::b0*logmu1/3.0;
      plusterms[2] = 4*asopi2/9.0*(-23 + 48*log1);
      plusterms[3] = 128*asopi2/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      asopi = as_n3loxs(mur, 3, asopimz);
      asopi2 = asopi*asopi;
      asopi3 = asopi*asopi2;
      log1 = log(Q2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 8*asopi/3.0*log1 +
	asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3) +
	asopi2*constants::b0*logmu1*8.0/3.0*log1 +
	asopi3/87480.0*
	(8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	 540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	 120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	 18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	 5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5)) +
	asopi3*(8*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)*log1/3.0 +
		2*constants::b0*logmu1*
		(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3)/81.0);
      plusterms[1] = 16*asopi/3.0 -
	4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) +
	16*asopi2*constants::b0*logmu1/3.0 -
	asopi3/2430.0*
	(45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	 1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	 20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3)) +
	asopi3*(16*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1)/3.0 -
		8*constants::b0*logmu1*
		(41 + 25*constants::Pi2 - 3*log1 - 48*log2)/27.0);
      plusterms[2] = 4*asopi2/9.0*(-23 + 48*log1) -
	2*asopi3/81.0*
	(5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 + 756*log2 - 384*log3 - 19896*constants::Zeta3) +
	8*asopi3*constants::b0*logmu1*(-23 + 48*log1)/9.0;
      plusterms[3] = 128*asopi2/9.0 -
	4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2) + 256*asopi3*constants::b0*logmu1/9.0;;
      plusterms[4] = 160*asopi3/81.0*(-23 + 24*log1);
      plusterms[5] = 512*asopi3/27.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = -fac*res;
  res = res*dlumqqb(x2,tau/x2,muf2,pdf)/x2;

  res = res/Den;
  return res;
}
