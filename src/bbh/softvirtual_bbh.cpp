/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 18/02/2021
Soft+Virtual contributions for the bbH production process up to N3LO QCD
*********************************************************************
********************************************************************* */

// pdf functions
#include "pdffunctions_bbh.h"

#include "bbh_functions.h"

#include "constants.h"


double intpow(const double& x,int m){
        double res=1.0;
        for (int i=0;i<m;i++){
            res *= x;
        }
        return res;
    }

static const double eps = 1.e-12;

//static const double MH2 = constants::MH*constants::MH;


// Virtual delta(z) contribution up to N3LO
double delta_bbh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double delterms;
  double z;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MH2 = constants::MH*constants::MH;

  tau = MH2/s;
  muf2 = muf*muf;

  delterms = 0.0;
  
  switch(k)
    {
    case 0: // LO
      delterms = 1.0;
      break;
    case 1: // NLO
      delterms = 4.0/9.0*(-3 + constants::Pi2);
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      delterms = 1.0/1620.0*
	(10950 + 1240*constants::Pi2 - 19*constants::Pi4 - 
	 960*constants::Pi2*log2 -8640*constants::Zeta3 - 
	 60*log1*(81 + 12*constants::Pi2 - 366*constants::Zeta3));
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      delterms =
	(-298626*constants::Pi4 - 46978*constants::Pi6 + 25920*log3*(23*constants::Pi2 + 128*constants::Zeta3) + 
	 540*constants::Pi2*(1596 + 20401*constants::Zeta3) - 108*log2*
	 (16180*constants::Pi2 + 352*constants::Pi4 + 3105*(-9 + 62*constants::Zeta3)) +
	 270*(195249 - 218630*constants::Zeta3 + 376584*constants::Zeta3*constants::Zeta3 - 476934*constants::Zeta5) + 
	 9*log1*(33*constants::Pi4 + 160*constants::Pi2*(217 - 11862*constants::Zeta3) + 
		 720*(-3717 + 9454*constants::Zeta3 + 19827*constants::Zeta5)))/5.2488e5;
      break;
    }

  z   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*z*log(tau);

  res = fac*delterms;
  res = res*tau*dlumbbb(z,tau/z,muf2,pdf)/z;
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
double PlusConst_bbh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double plusterms[6];
  double x2;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MH2 = constants::MH*constants::MH;

  tau = MH2/s;
  muf2 = muf*muf;

  x2   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*x2*log(tau);

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      log1 = log(MH2/muf2);
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]*log(1.0-tau) + plusterms[1]*intpow(log(1.0-tau),2)/2.0;
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 + (618 - 150*constants::Pi2)*log1 - 207*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(-103 + 25*constants::Pi2 + 69*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0;
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 1.0/87480.0*
      (8298*constants::Pi4 - 540*log3*(-529 + 256*constants::Pi2) - 
       600*constants::Pi2*(-3376 + 9045*constants::Zeta3) + 
       540*log2*(-4739 + 935*constants::Pi2 + 12000*constants::Zeta3) - 
       18*log1*(-380605 + 101800*constants::Pi2 + 1218*constants::Pi4 + 1141560*constants::Zeta3) + 
       5*(-1203875 + 6499440*constants::Zeta3 + 6702912*constants::Zeta5));
      plusterms[1] = -1.0/2430.0*
	(-380605 + 33120*log3 + 101800*constants::Pi2 + 1218*constants::Pi4 + 
	 30*log2*(-5651 + 1312*constants::Pi2) + 1141560*constants::Zeta3 - 
	 20*log1*(-19729 + 4197*constants::Pi2 + 54720*constants::Zeta3));
      plusterms[2] = -2.0/81.0*
	(7171 + 2484*log2 - 384*log3 - 1679*constants::Pi2 + log1*(-7683 + 1584*constants::Pi2) - 19896*constants::Zeta3);
      plusterms[3] = -(4.0/81.0)*(-2561 + 1840*log1 - 768*log2 + 528*constants::Pi2);
      plusterms[4] = 160.0/81.0*(-23 + 24*log1);
      plusterms[5] = 512.0/27.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0 +
	plusterms[4]*intpow(log(1.0-tau),5)/5.0 +
	plusterms[5]*intpow(log(1.0-tau),6)/6.0;
      break;
    }
  res = fac*res;
  res = res*tau*dlumbbb(x2,tau/x2,muf2,pdf)/x2;
  return res;
}

// PlusInt1 term
double PlusInt1_bbh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MH2 = constants::MH*constants::MH;

  tau = MH2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      log1 = log(MH2/muf2);
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      //      std::cout << "DEBUG\t" << res << std::endl;
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 + (618 - 150*constants::Pi2)*log1 - 207*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(-103 + 25*constants::Pi2 + 69*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 1.0/87480.0*
      (8298*constants::Pi4 - 540*log3*(-529 + 256*constants::Pi2) - 
       600*constants::Pi2*(-3376 + 9045*constants::Zeta3) + 
       540*log2*(-4739 + 935*constants::Pi2 + 12000*constants::Zeta3) - 
       18*log1*(-380605 + 101800*constants::Pi2 + 1218*constants::Pi4 + 1141560*constants::Zeta3) + 
       5*(-1203875 + 6499440*constants::Zeta3 + 6702912*constants::Zeta5));
      plusterms[1] = -1.0/2430.0*
	(-380605 + 33120*log3 + 101800*constants::Pi2 + 1218*constants::Pi4 + 
	 30*log2*(-5651 + 1312*constants::Pi2) + 1141560*constants::Zeta3 - 
	 20*log1*(-19729 + 4197*constants::Pi2 + 54720*constants::Zeta3));
      plusterms[2] = -2.0/81.0*
	(7171 + 2484*log2 - 384*log3 - 1679*constants::Pi2 + log1*(-7683 + 1584*constants::Pi2) - 19896*constants::Zeta3);
      plusterms[3] = -(4.0/81.0)*(-2561 + 1840*log1 - 768*log2 + 528*constants::Pi2);
      plusterms[4] = 160.0/81.0*(-23 + 24*log1);
      plusterms[5] = 512.0/27.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = fac*res;
  res = res*tau*( dlumbbb(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumbbb(x2,tau/x2,muf2,pdf)/x2 );
  return res;
}


// PlusInt2 term
double PlusInt2_bbh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2;
  double log1,log2,log3;
  const double MH2 = constants::MH*constants::MH;

  tau = MH2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau*exp(-(eps+(1.0-2.0*eps)*X[1])*log(x1));
  fac = intpow(1.0-2.0*eps,2)*x1*log(tau)*x2*log(x1);

  res = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      log1 = log(MH2/muf2);
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 + (618 - 150*constants::Pi2)*log1 - 207*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(-103 + 25*constants::Pi2 + 69*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = 1.0/87480.0*
      (8298*constants::Pi4 - 540*log3*(-529 + 256*constants::Pi2) - 
       600*constants::Pi2*(-3376 + 9045*constants::Zeta3) + 
       540*log2*(-4739 + 935*constants::Pi2 + 12000*constants::Zeta3) - 
       18*log1*(-380605 + 101800*constants::Pi2 + 1218*constants::Pi4 + 1141560*constants::Zeta3) + 
       5*(-1203875 + 6499440*constants::Zeta3 + 6702912*constants::Zeta5));
      plusterms[1] = -1.0/2430.0*
	(-380605 + 33120*log3 + 101800*constants::Pi2 + 1218*constants::Pi4 + 
	 30*log2*(-5651 + 1312*constants::Pi2) + 1141560*constants::Zeta3 - 
	 20*log1*(-19729 + 4197*constants::Pi2 + 54720*constants::Zeta3));
      plusterms[2] = -2.0/81.0*
	(7171 + 2484*log2 - 384*log3 - 1679*constants::Pi2 + log1*(-7683 + 1584*constants::Pi2) - 19896*constants::Zeta3);
      plusterms[3] = -(4.0/81.0)*(-2561 + 1840*log1 - 768*log2 + 528*constants::Pi2);
      plusterms[4] = 160.0/81.0*(-23 + 24*log1);
      plusterms[5] = 512.0/27.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = -fac*res;
  res = res*tau*dlumbbb(x2,tau/x2,muf2,pdf)/x2;
  return res;
}
