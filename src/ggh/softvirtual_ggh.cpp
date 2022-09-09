/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 24/02/2021
Soft+Virtual contributions for the ggH production process up to N3LO QCD
in the heavy top-quark limit, without the Wilson coefficient
Based on arXiv:hep-ph/0302135, arXiv:1403.4616, and arXiv:1802.00833
*********************************************************************
********************************************************************* */

// pdf functions
#include "pdffunctions_ggh.h"

#include "ggh_functions.h"

#include "constants.h"

static const double eps = 1.e-12;


double intpow(const double& x,int m){
        double res=1.0;
        for (int i=0;i<m;i++){
            res *= x;
        }
        return res;
    }


// Virtual delta(z) contribution up to N3LO
double delta_ggh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
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
      delterms = constants::Pi2;
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      delterms = -3*log2*constants::Pi2 + 
	log1*(constants::nf*(11 + constants::Pi2)/6.0 +
	      (-54 - 11*constants::Pi2 + 342*constants::Zeta3)/4.0) + 
	(37665 + 4020*constants::Pi2 - 9*constants::Pi4 -
	 20*constants::nf*(247 + 10*constants::Pi2 - 30*constants::Zeta3) - 29700*constants::Zeta3)/720.0;      
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      delterms = log3*(constants::Pi2*(33 - 2*constants::nf)/4.0 + 72*constants::Zeta3) + 
      	log2*(20*constants::nf*constants::nf*(11 + constants::Pi2) +
      	      30*constants::nf*(-175 + 18*constants::Pi2 + 486*constants::Zeta3) - 
      	      9*(2075*constants::Pi2 + 216*constants::Pi4 + 2970*(-1 + 9*constants::Zeta3)))/720.0 + 
      	log1*(6*constants::nf*(21925 + 1248*constants::Pi2 + 15*constants::Pi4 - 18960*constants::Zeta3) - 
      	      2*constants::nf*constants::nf*(2289 + 80*constants::Pi2 - 240*constants::Zeta3) - 
      	      9*(8*constants::Pi2*(397 + 6534*constants::Zeta3) + 165*constants::Pi4 - 
      		 36*(-2071 + 6302*constants::Zeta3 + 13464*constants::Zeta5)))/1728.0 - 
      	(35*constants::nf*constants::nf*
      	 (4296*constants::Pi2 + 280*constants::Pi4 - 103753 + 5616*constants::Zeta3) - 
      	 105*constants::nf*
      	 (-1128767 + 3190*constants::Pi4 + 429096*constants::Zeta3 - 
      	  72*constants::Pi2*(188 + 981*constants::Zeta3) + 568872*constants::Zeta5) + 
      	 27*(106799*constants::Pi4 + 24036*constants::Pi6 -
      	     70*constants::Pi2*(16151 + 52866*constants::Zeta3) - 
      	     735*(30733 + 108*constants::Zeta3*(-351 + 472*constants::Zeta3)) + 28648620*constants::Zeta5)
      	 )/544320.0;      
      break;
    }

  z   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*z*log(tau);

  res = fac*delterms;
  res = res*tau*dlumgg(z,tau/z,muf2,pdf)/z;
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
double PlusConst_ggh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
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
      plusterms[0] = 6.0*log1;
      plusterms[1] = 12.0;

      res = plusterms[0]*log(1.0-tau) + plusterms[1]*intpow(log(1.0-tau),2)/2.0;
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      plusterms[0] = (9*log2*(-33 + 2*constants::nf) -
		      6*log1*(-201 + 10*constants::nf + 45*constants::Pi2) + 
		      2*(-606 + 28*constants::nf + 99*constants::Pi2 -
			 6*constants::nf*constants::Pi2 + 3159*constants::Zeta3))/36.0;
      plusterms[1] = 67 - 33*log1 + 36*log2 + constants::nf*(2*log1 -10.0/3.0) - 15*constants::Pi2;
      plusterms[2] = -33.0 + 108*log1 + 2*constants::nf;
      plusterms[3] = 72.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0;
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = log3*((33.0 - 2*constants::nf)*(33.0 - 2*constants::nf)/72.0 - 18*constants::Pi2) + 
      	log2*(-769.0/4.0 + 231*constants::Pi2/8.0 -
      	      constants::nf*(-1695.0 + 20*constants::nf + 126*constants::Pi2)/72.0 + 945*constants::Zeta3) + 
      	log1*(30569.0/48.0 - 5345*constants::nf/72.0 - 114*constants::Pi2 + 47*constants::nf*constants::Pi2/6.0 -
      	      231*constants::Pi4/20.0 + constants::nf*constants::nf*(25.0 - 6*constants::Pi2)/54.0 + 
      	      9*(-176.0 + 9*constants::nf)*constants::Zeta3) -
      	297029.0/864.0 + 253*constants::Pi4/80.0 +
      	constants::Pi2*(8563.0/72.0 - 2175*constants::Zeta3/4.0) + 8941*constants::Zeta3/4.0 + 
      	constants::nf*constants::nf*(-58.0 + 45*constants::Pi2 + 135*constants::Zeta3)/243.0 -
      	constants::nf*(-207895.0 + 67350*constants::Pi2 + 1062*constants::Pi4 + 796860*constants::Zeta3)/6480.0  +
      	5022*constants::Zeta5;
      plusterms[1] = log3*(-99.0 + 6*constants::nf) +
      	log2*(1971.0/4.0 + constants::nf*(-93.0 + constants::nf)/3.0 - 162*constants::Pi2) +
      	log1*(-1011.0 - 10*constants::nf*constants::nf/9.0 + 429*constants::Pi2/2.0 + 
      	      constants::nf*(545.0/6.0 - 13*constants::Pi2) + 4860*constants::Zeta3) +
      	30569.0/24.0 - 5345*constants::nf/36.0  - 228*constants::Pi2 +
      	47*constants::nf*constants::Pi2/3.0 - 231*constants::Pi4/10.0  + 
      	constants::nf*constants::nf*(25.0 - 6*constants::Pi2)/27.0 + 18*(-176.0 + 9*constants::nf)*constants::Zeta3;
      plusterms[2] =  108*log3 - log2*(891.0/2.0 - 27*constants::nf) +
      	log1*(2775.0/2.0 + 2*constants::nf*(-123.0 + constants::nf)/3.0 - 378*constants::Pi2)
      	-1051  - 10*constants::nf*constants::nf/9.0 + 561*constants::Pi2/2.0  + 
      	constants::nf*(469.0/6.0 - 17*constants::Pi2) + 4887*constants::Zeta3;
      plusterms[3] = 432*log2 - (660.0 - 40.0*constants::nf)*log1 +
      	925.0 + 4*(-123.0 + constants::nf)*constants::nf/9.0- 252*constants::Pi2;
      plusterms[4] = 540*log1 - 330.0 + 20*constants::nf;
      plusterms[5] = 216.0;
      
      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0 +
	plusterms[4]*intpow(log(1.0-tau),5)/5.0 +
	plusterms[5]*intpow(log(1.0-tau),6)/6.0;
      break;
    }
  res = fac*res;
  res = res*tau*dlumgg(x2,tau/x2,muf2,pdf)/x2;
  return res;
}

// PlusInt1 term
double PlusInt1_ggh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
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
      plusterms[0] = 6.0*log1;
      plusterms[1] = 12.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      plusterms[0] = (9*log2*(-33 + 2*constants::nf) -
		      6*log1*(-201 + 10*constants::nf + 45*constants::Pi2) + 
		      2*(-606 + 28*constants::nf + 99*constants::Pi2 -
			 6*constants::nf*constants::Pi2 + 3159*constants::Zeta3))/36.0;
      plusterms[1] = 67 - 33*log1 + 36*log2 + constants::nf*(2*log1 -10.0/3.0) - 15*constants::Pi2;
      plusterms[2] = -33.0 + 108*log1 + 2*constants::nf;
      plusterms[3] = 72.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
      plusterms[0] = log3*((33.0 - 2*constants::nf)*(33.0 - 2*constants::nf)/72.0 - 18*constants::Pi2) + 
      	log2*(-769.0/4.0 + 231*constants::Pi2/8.0 -
      	      constants::nf*(-1695.0 + 20*constants::nf + 126*constants::Pi2)/72.0 + 945*constants::Zeta3) + 
      	log1*(30569.0/48.0 - 5345*constants::nf/72.0 - 114*constants::Pi2 + 47*constants::nf*constants::Pi2/6.0 -
      	      231*constants::Pi4/20.0 + constants::nf*constants::nf*(25.0 - 6*constants::Pi2)/54.0 + 
      	      9*(-176.0 + 9*constants::nf)*constants::Zeta3) -
      	297029.0/864.0 + 253*constants::Pi4/80.0 +
      	constants::Pi2*(8563.0/72.0 - 2175*constants::Zeta3/4.0) + 8941*constants::Zeta3/4.0 + 
      	constants::nf*constants::nf*(-58.0 + 45*constants::Pi2 + 135*constants::Zeta3)/243.0 -
      	constants::nf*(-207895.0 + 67350*constants::Pi2 + 1062*constants::Pi4 + 796860*constants::Zeta3)/6480.0  +
      	5022*constants::Zeta5;
      plusterms[1] = log3*(-99.0 + 6*constants::nf) +
      	log2*(1971.0/4.0 + constants::nf*(-93.0 + constants::nf)/3.0 - 162*constants::Pi2) +
      	log1*(-1011.0 - 10*constants::nf*constants::nf/9.0 + 429*constants::Pi2/2.0 + 
      	      constants::nf*(545.0/6.0 - 13*constants::Pi2) + 4860*constants::Zeta3) +
      	30569.0/24.0 - 5345*constants::nf/36.0  - 228*constants::Pi2 +
      	47*constants::nf*constants::Pi2/3.0 - 231*constants::Pi4/10.0  + 
      	constants::nf*constants::nf*(25.0 - 6*constants::Pi2)/27.0 + 18*(-176.0 + 9*constants::nf)*constants::Zeta3;
      plusterms[2] =  108*log3 - log2*(891.0/2.0 - 27*constants::nf) +
      	log1*(2775.0/2.0 + 2*constants::nf*(-123.0 + constants::nf)/3.0 - 378*constants::Pi2)
      	-1051  - 10*constants::nf*constants::nf/9.0 + 561*constants::Pi2/2.0  + 
      	constants::nf*(469.0/6.0 - 17*constants::Pi2) + 4887*constants::Zeta3;
      plusterms[3] = 432*log2 - (660.0 - 40.0*constants::nf)*log1 +
      	925.0 + 4*(-123.0 + constants::nf)*constants::nf/9.0- 252*constants::Pi2;
      plusterms[4] = 540*log1 - 330.0 + 20*constants::nf;
      plusterms[5] = 216.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = fac*res;
  res = res*tau*( dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumgg(x2,tau/x2,muf2,pdf)/x2 );
  return res;
}


// PlusInt2 term
double PlusInt2_ggh(const double X[], const double s, const double muf, const int k, LHAPDF::PDF const* const pdf)
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
      plusterms[0] = 6.0*log1;
      plusterms[1] = 12.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      plusterms[0] = (9*log2*(-33 + 2*constants::nf) -
		      6*log1*(-201 + 10*constants::nf + 45*constants::Pi2) + 
		      2*(-606 + 28*constants::nf + 99*constants::Pi2 -
			 6*constants::nf*constants::Pi2 + 3159*constants::Zeta3))/36.0;
      plusterms[1] = 67 - 33*log1 + 36*log2 + constants::nf*(2*log1 -10.0/3.0) - 15*constants::Pi2;
      plusterms[2] = -33.0 + 108*log1 + 2*constants::nf;
      plusterms[3] = 72.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
      log1 = log(MH2/muf2);
      log2 = log1*log1;
      log3 = log1*log2;
            plusterms[0] = log3*((33.0 - 2*constants::nf)*(33.0 - 2*constants::nf)/72.0 - 18*constants::Pi2) + 
      	log2*(-769.0/4.0 + 231*constants::Pi2/8.0 -
      	      constants::nf*(-1695.0 + 20*constants::nf + 126*constants::Pi2)/72.0 + 945*constants::Zeta3) + 
      	log1*(30569.0/48.0 - 5345*constants::nf/72.0 - 114*constants::Pi2 + 47*constants::nf*constants::Pi2/6.0 -
      	      231*constants::Pi4/20.0 + constants::nf*constants::nf*(25.0 - 6*constants::Pi2)/54.0 + 
      	      9*(-176.0 + 9*constants::nf)*constants::Zeta3) -
      	297029.0/864.0 + 253*constants::Pi4/80.0 +
      	constants::Pi2*(8563.0/72.0 - 2175*constants::Zeta3/4.0) + 8941*constants::Zeta3/4.0 + 
      	constants::nf*constants::nf*(-58.0 + 45*constants::Pi2 + 135*constants::Zeta3)/243.0 -
      	constants::nf*(-207895.0 + 67350*constants::Pi2 + 1062*constants::Pi4 + 796860*constants::Zeta3)/6480.0  +
      	5022*constants::Zeta5;
      plusterms[1] = log3*(-99.0 + 6*constants::nf) +
      	log2*(1971.0/4.0 + constants::nf*(-93.0 + constants::nf)/3.0 - 162*constants::Pi2) +
      	log1*(-1011.0 - 10*constants::nf*constants::nf/9.0 + 429*constants::Pi2/2.0 + 
      	      constants::nf*(545.0/6.0 - 13*constants::Pi2) + 4860*constants::Zeta3) +
      	30569.0/24.0 - 5345*constants::nf/36.0  - 228*constants::Pi2 +
      	47*constants::nf*constants::Pi2/3.0 - 231*constants::Pi4/10.0  + 
      	constants::nf*constants::nf*(25.0 - 6*constants::Pi2)/27.0 + 18*(-176.0 + 9*constants::nf)*constants::Zeta3;
      plusterms[2] =  108*log3 - log2*(891.0/2.0 - 27*constants::nf) +
      	log1*(2775.0/2.0 + 2*constants::nf*(-123.0 + constants::nf)/3.0 - 378*constants::Pi2)
      	-1051  - 10*constants::nf*constants::nf/9.0 + 561*constants::Pi2/2.0  + 
      	constants::nf*(469.0/6.0 - 17*constants::Pi2) + 4887*constants::Zeta3;
      plusterms[3] = 432*log2 - (660.0 - 40.0*constants::nf)*log1 +
      	925.0 + 4*(-123.0 + constants::nf)*constants::nf/9.0- 252*constants::Pi2;
      plusterms[4] = 540*log1 - 330.0 + 20*constants::nf;
      plusterms[5] = 216.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);
      break;
    }
  res = -fac*res;
  res = res*tau*dlumgg(x2,tau/x2,muf2,pdf)/x2;
  return res;
}
