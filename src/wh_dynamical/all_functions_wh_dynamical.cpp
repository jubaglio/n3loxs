/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 20/05/2021
Regular hard terms for all W- H subprocesses (DY-type) up to N3LO QCD (dynamical scale)
*********************************************************************
********************************************************************* */

// pdf functions
#include "pdffunctions_w.h"

#include "dy_w_kernels.h"
#include "dy_functions_wh_dyn.h"

#include "constants.h"

// Header for the routines alphaS(muR)
#include "alphaS.h"

static const double eps = 1.e-8;


double intpow(const double& x,int m){
        double res=1.0;
        for (int i=0;i<m;i++){
            res *= x;
        }
        return res;
    }


double gammalo_wh(const double Q2)
{
  double lambda, result, Born;

  lambda = (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MW*constants::MW/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MW*constants::MW/(Q2*Q2);

  Born = constants::gevtopb*constants::MW*constants::MW*constants::MW*constants::MW
    /(48*constants::Pi*constants::Nc*constants::vev*constants::vev*constants::vev*constants::vev);
  
  result = Born*(Q2*lambda+12.0*constants::MW*constants::MW)*
    sqrt(lambda)/((Q2-constants::MW*constants::MW)*(Q2-constants::MW*constants::MW));

  return result;
}


// Virtual delta(z) contribution up to N3LO
double delta(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, Q2;
  double delterms;
  double z;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1, log2, log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[1])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);

  if(k>0) {
    asopi  = as_n3loxs(mur, k, asopimz);
    log1 = log(Q2/muf2);
  }
  if(k>1) {
    asopi2 = asopi*asopi;
    log2 = log1*log1;
    logmu1 = log(mur2/muf2);
  }
  if(k>2) {
    asopi3 = asopi2*asopi;
    log3 = log1*log2;
    logmu2 = logmu1*logmu1;
  }

  z   = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  fac = -(1.0-2.0*eps)*z*log(tau)*fac;

  if(k==0) {
    res = fac*gammalo_wh(Q2)*dlumdub(z,tau/z,muf2,pdf)/z;
  } else if(k==1) {
    res = fac*gammalo_wh(Q2)*
      (1.0 + 2*asopi/9.0*(-24.0 + 2*constants::Pi2 + 9*log1))*
      dlumdub(z,tau/z,muf2,pdf)/z;
  } else if(k==2) {
    res = fac*gammalo_wh(Q2)*
      (1.0 + 2*(asopi+constants::b0*logmu1*asopi2)/9.0*(-24.0 + 2*constants::Pi2 + 9*log1) +
       asopi2/6480.0*
       (-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3))
       )*dlumdub(z,tau/z,muf2,pdf)/z;
  } else {
    res = fac*gammalo_wh(Q2)*
      (1.0 + 2*(asopi + constants::b0*logmu1*asopi2 +
		asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))/9.0*(-24.0 + 2*constants::Pi2 + 9*log1) +
       (asopi2 + 2*asopi3*constants::b0*logmu1)/6480.0*
       (-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3)) +
       asopi3*(684831*constants::Pi4 - 187912*constants::Pi6 - 3240*log3*(33 + 32*constants::Pi2 - 4096*constants::Zeta3) + 
	       270*constants::Pi2*(47045 + 195068*constants::Zeta3) - 108*log2*
	       (36675 + 35200*constants::Pi2 + 1408*constants::Pi4 + 243000*constants::Zeta3) +
	       45*(-2117173 + 461868*constants::Zeta3 + 9038016*constants::Zeta3*constants::Zeta3 - 12599856*constants::Zeta5) - 
	       180*log1*(267*constants::Pi4 + 8*constants::Pi2*(-5773 + 47448*constants::Zeta3) - 
			 72*(2547 + 4373*constants::Zeta3 + 39654*constants::Zeta5)))/2.09952e6
       )*dlumdub(z,tau/z,muf2,pdf)/z;
  }

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
double PlusConst(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, Q2;
  double x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1, log2, log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[1])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  if(k>0) {
    asopi  = as_n3loxs(mur, k, asopimz);
    log1 = log(Q2/muf2);
  }
  if(k>1) {
    asopi2 = asopi*asopi;
    log2 = log1*log1;
    logmu1 = log(mur2/muf2);
  }
  if(k>2) {
    asopi3 = asopi2*asopi;
    log3 = log1*log2;
    logmu2 = logmu1*logmu1;
  }
  
  x2   = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  fac = -(1.0-2.0*eps)*x2*log(tau)*fac;

  if(k==0) {
    std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
    exit(1);
  } else if(k==1) {
    res = asopi*(8/3.0*log1*log(1.0-tau) + 16/3.0*intpow(log(1.0-tau),2)/2.0);
  } else if(k==2) {
    res = ((asopi + asopi2*constants::b0*logmu1)*8.0/3.0*log1 +
	   asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
			225*log2 + 3438*constants::Zeta3)
	   )*log(1.0-tau) +
      ((asopi + asopi2*constants::b0*logmu1)*16.0/3.0 -
       4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2)
       )*intpow(log(1.0-tau),2)/2.0 +
      4*asopi2/9.0*(-23 + 48*log1)*intpow(log(1.0-tau),3)/3.0 +
      128*asopi2/9.0*intpow(log(1.0-tau),4)/4.0;
  } else {
    res = ((asopi +
	    asopi2*constants::b0*logmu1 +
	    asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*8.0/3.0*log1 +
	   (asopi2 +
	    2*asopi3*constants::b0*logmu1)/81.0*
	   (-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
	    225*log2 + 3438*constants::Zeta3) +
	   asopi3/87480.0*
	   (8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	    540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	    120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	    18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	    5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5))
	   )*log(1.0-tau) +
      ((asopi +
	asopi2*constants::b0*logmu1 +
	asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*16.0/3.0 -
       4*(asopi2 +
	  2*asopi3*constants::b0*logmu1)/27.0*
       (41 + 25*constants::Pi2 - 3*log1 - 48*log2) -
       asopi3/2430.0*
       (45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3))
       )*intpow(log(1.0-tau),2)/2.0 +
      (4*(asopi2 +
	  2*asopi3*constants::b0*logmu1)/9.0*(-23 + 48*log1) -
       2*asopi3/81.0*
       (5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 +
	756*log2 - 384*log3 - 19896*constants::Zeta3)
       )*intpow(log(1.0-tau),3)/3.0 +
      (128*(asopi2 +
	    2*asopi3*constants::b0*logmu1)/9.0 -
       4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2)
       )*intpow(log(1.0-tau),4)/4.0 +
      160*asopi3/81.0*(-23 + 24*log1)*intpow(log(1.0-tau),5)/5.0 +
      512*asopi3/27.0*intpow(log(1.0-tau),6)/6.0;
  }

  res = res*fac*gammalo_wh(Q2)*dlumdub(x2,tau/x2,muf2,pdf)/x2;

  return res;
}

// PlusInt1 term
double PlusInt1(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, Q2;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1,log2,log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  if(k>0) {
    asopi  = as_n3loxs(mur, k, asopimz);
    log1 = log(Q2/muf2);
  }
  if(k>1) {
    asopi2 = asopi*asopi;
    log2 = log1*log1;
    logmu1 = log(mur2/muf2);
  }
  if(k>2) {
    asopi3 = asopi2*asopi;
    log3 = log1*log2;
    logmu2 = logmu1*logmu1;
  }

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  if(k==0) {
    std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
    exit(1);
  } else if(k==1) {
    res = asopi*(8/3.0*log1/(1.0-x1) + 16/3.0*log(1.0-x1)/(1.0-x1));
  } else if(k==2) {
    res = ((asopi + asopi2*constants::b0*logmu1)*8.0/3.0*log1 +
	   asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
			225*log2 + 3438*constants::Zeta3)
	   )/(1.0-x1) +
      ((asopi + asopi2*constants::b0*logmu1)*16.0/3.0 -
       4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2)
       )*log(1.0-x1)/(1.0-x1) +
      4*asopi2/9.0*(-23 + 48*log1)*intpow(log(1.0-x1),2)/(1.0-x1) +
      128*asopi2/9.0*intpow(log(1.0-x1),3)/(1.0-x1);
  } else {
    res = ((asopi +
	    asopi2*constants::b0*logmu1 +
	    asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*8.0/3.0*log1 +
	   (asopi2 +
	    2*asopi3*constants::b0*logmu1)/81.0*
	   (-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
	    225*log2 + 3438*constants::Zeta3) +
	   asopi3/87480.0*
	   (8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	    540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	    120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	    18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	    5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5))
	   )/(1.0-x1) +
      ((asopi +
	asopi2*constants::b0*logmu1 +
	asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*16.0/3.0 -
       4*(asopi2 +
	  2*asopi3*constants::b0*logmu1)/27.0*
       (41 + 25*constants::Pi2 - 3*log1 - 48*log2) -
       asopi3/2430.0*
       (45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3))
       )*log(1.0-x1)/(1.0-x1) +
      (4*(asopi2 +
	  2*asopi3*constants::b0*logmu1)/9.0*(-23 + 48*log1) -
       2*asopi3/81.0*
       (5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 +
	756*log2 - 384*log3 - 19896*constants::Zeta3)
       )*intpow(log(1.0-x1),2)/(1.0-x1) +
      (128*(asopi2 +
	    2*asopi3*constants::b0*logmu1)/9.0 -
       4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2)
       )*intpow(log(1.0-x1),3)/(1.0-x1) +
      160*asopi3/81.0*(-23 + 24*log1)*intpow(log(1.0-x1),4)/(1.0-x1) +
      512*asopi3/27.0*intpow(log(1.0-x1),5)/(1.0-x1);
  }

  res = res*fac*gammalo_wh(Q2)*( dlumdub(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumdub(x2,tau/x2,muf2,pdf)/x2 );

  return res;
}

// PlusInt2 term
double PlusInt2(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauwh, Q2;
  double plusterms[6];
  double x1,x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1,log2,log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  if(k>0) {
    asopi  = as_n3loxs(mur, k, asopimz);
    log1 = log(Q2/muf2);
  }
  if(k>1) {
    asopi2 = asopi*asopi;
    log2 = log1*log1;
    logmu1 = log(mur2/muf2);
  }
  if(k>2) {
    asopi3 = asopi2*asopi;
    log3 = log1*log2;
    logmu2 = logmu1*logmu1;
  }

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau*exp(-(eps+(1.0-2.0*eps)*X[1])*log(x1));
  fac = intpow(1.0-2.0*eps,2)*x1*log(tau)*x2*log(x1)*fac;

  if(k==0) {
    std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
    exit(1);
  } else if(k==1) {
    res = asopi*(8/3.0*log1/(1.0-x1) + 16/3.0*log(1.0-x1)/(1.0-x1));
  } else if(k==2) {
    res = ((asopi + asopi2*constants::b0*logmu1)*8.0/3.0*log1 +
	   asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
			225*log2 + 3438*constants::Zeta3)
	   )/(1.0-x1) +
      ((asopi + asopi2*constants::b0*logmu1)*16.0/3.0 -
       4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2)
       )*log(1.0-x1)/(1.0-x1) +
      4*asopi2/9.0*(-23 + 48*log1)*intpow(log(1.0-x1),2)/(1.0-x1) +
      128*asopi2/9.0*intpow(log(1.0-x1),3)/(1.0-x1);
  } else {
    res = ((asopi +
	    asopi2*constants::b0*logmu1 +
	    asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*8.0/3.0*log1 +
	   (asopi2 +
	    2*asopi3*constants::b0*logmu1)/81.0*
	   (-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
	    225*log2 + 3438*constants::Zeta3) +
	   asopi3/87480.0*
	   (8298*constants::Pi4 - 540*(263 + 256*constants::Pi2)*log3 + 
	    540*log2*(2281 + 335*constants::Pi2 + 12000*constants::Zeta3) - 
	    120*constants::Pi2*(-11912 + 45225*constants::Zeta3) -
	    18*log1*(152465 + 28600*constants::Pi2 + 1218*constants::Pi4 + 612360*constants::Zeta3) +
	    5*(-398627 + 3529008*constants::Zeta3 + 6702912*constants::Zeta5))
	   )/(1.0-x1) +
      ((asopi +
	asopi2*constants::b0*logmu1 +
	asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*16.0/3.0 -
       4*(asopi2 +
	  2*asopi3*constants::b0*logmu1)/27.0*
       (41 + 25*constants::Pi2 - 3*log1 - 48*log2) -
       asopi3/2430.0*
       (45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3))
       )*log(1.0-x1)/(1.0-x1) +
      (4*(asopi2 +
	  2*asopi3*constants::b0*logmu1)/9.0*(-23 + 48*log1) -
       2*asopi3/81.0*
       (5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 +
	756*log2 - 384*log3 - 19896*constants::Zeta3)
       )*intpow(log(1.0-x1),2)/(1.0-x1) +
      (128*(asopi2 +
	    2*asopi3*constants::b0*logmu1)/9.0 -
       4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2)
       )*intpow(log(1.0-x1),3)/(1.0-x1) +
      160*asopi3/81.0*(-23 + 24*log1)*intpow(log(1.0-x1),4)/(1.0-x1) +
      512*asopi3/27.0*intpow(log(1.0-x1),5)/(1.0-x1);
  }

  res = -res*fac*gammalo_wh(Q2)*dlumdub(x2,tau/x2,muf2,pdf)/x2;

  return res;
}


//////////////////////////////////////////////////
///////////////// d-ubar channel /////////////////
//////////////////////////////////////////////////

// NLO d-ubar regular term
double dub_regular_nlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 1, asopimz);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi*dub_regular_kernel_nlo(x1, log1);

  res = res*dlumdub(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


// NNLO d-ubar regular term
double dub_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2;
  double res_nlo, res1, res2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res_nlo = (asopi + asopi2*constants::b0*logmu1)*dub_regular_kernel_nlo(x1, log1);

  auto result0 = dub_regular_kernel_nnlo(x1, log1);
  res1 = asopi2*std::get<0>(result0);
  res2 = asopi2*std::get<1>(result0);
  
  res = ((res_nlo + res1)*dlumdub(x2,tau/x1/x2,muf2,pdf) +
	 res2*dlumdub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO d-ubar regular term
double dub_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  double res_nlo, res1, res2, res1_n3lo, res2_n3lo;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res_nlo = (asopi +
	     asopi2*constants::b0*logmu1 +
	     asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*dub_regular_kernel_nlo(x1, log1);
  
  auto result0 = dub_regular_kernel_nnlo(x1, log1);
  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  res2 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  auto result1 = dub_regular_kernel_n3lo(x1, log1);
  res1_n3lo = asopi3*std::get<0>(result1);
  res2_n3lo = asopi3*std::get<1>(result1);

  res = ((res_nlo + res1 + res1_n3lo)*dlumdub(x2,tau/x1/x2,muf2,pdf) +
	 (res2 + res2_n3lo)*dlumdub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
///////////////// g-ubar channel /////////////////
//////////////////////////////////////////////////

// NLO g-ub + d-g regular term
double gub_regular_nlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 1, asopimz);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi*gub_regular_kernel_nlo(x1, log1);

  res = res*dlumgub(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// NNLO g-ub + d-g regular term
double gub_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double res1, res2;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi + asopi2*constants::b0*logmu1)*gub_regular_kernel_nlo(x1, log1);
  res2 = asopi2*gub_regular_kernel_nnlo(x1, log1);

  res = (res1 + res2)*dlumgub(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO g-ub + d-g regular term
double gub_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  double res_nlo, res_nnlo, res1, res2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res_nlo = (asopi +
	     asopi2*constants::b0*logmu1 +
	     asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*gub_regular_kernel_nlo(x1, log1);
  res_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*gub_regular_kernel_nnlo(x1, log1);
  
  auto result0 = gub_regular_kernel_n3lo(x1, log1);
  res1 = asopi3*std::get<0>(result0);
  res2 = asopi3*std::get<1>(result0);

  res = ((res_nlo + res_nnlo + res1)*dlumgub(x2,tau/x1/x2,muf2,pdf) +
	 res2*dlumgub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
////////////////// g-g channel ///////////////////
//////////////////////////////////////////////////

// NNLO g-g regular term
double gg_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*gg_regular_kernel_nnlo(x1, log1);
  
  res = res*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO g-g regular term
double gg_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double res1, res2;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*gg_regular_kernel_nnlo(x1, log1);
  res2 = asopi3*gg_regular_kernel_n3lo(x1, log1);
  
  res = (res1 + res2)*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
//////////////// g-dbar channel //////////////////
//////////////////////////////////////////////////

// N3LO g-dbar + dbar-g regular term
double gdb_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1, res2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = gdb_regular_kernel_n3lo(x1, log1);
  res1 = asopi3*std::get<0>(result0);
  res2 = asopi3*std::get<1>(result0);
  
  res = (res1*dlumgu(x2,tau/x1/x2,muf2,pdf) +
	 res2*dlumgu2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
//////////////// c-ubar channel //////////////////
//////////////////////////////////////////////////

// NNLO c-ubar + u-cbar regular term
double cub_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*cub_regular_kernel_nnlo(x1, log1);

  res = res*(dlumcub(x2,tau/x1/x2,muf2,pdf)+dlumcub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO c-ubar + u-cbar regular term
double cub_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*cub_regular_kernel_nnlo(x1, log1);
  res2 = asopi3*cub_regular_kernel_n3lo(x1, log1);
  
  res = (res1 + res2)*(dlumcub(x2,tau/x1/x2,muf2,pdf)+dlumcub2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
//////////////// q-qbar channel //////////////////
//////////////////////////////////////////////////

// NNLO q-qbar regular term
double qqb_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1,res2;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  res1 = asopi2*std::get<0>(result0);
  res2 = asopi2*std::get<1>(result0);

  res = (res1*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	 res2*dlumqqb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO q-qbar regular term
double qqb_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1,res2, res3, res4;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  res2 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  auto result1 = qqb_regular_kernel_n3lo(x1, log1);
  res3 = asopi3*std::get<0>(result1);
  res4 = asopi3*std::get<1>(result1);

  res = ((res1 + res3)*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	 (res2 + res4)*dlumqqb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
//////////// q-q / qbar-qbar channel /////////////
//////////////////////////////////////////////////

// NNLO q-q + qbar-qbar regular term
double qq_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*qq_regular_kernel_nnlo(x1, log1);

  res = res*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO q-q regular term
double qq_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*qq_regular_kernel_nnlo(x1, log1);
  res2 = asopi3*qq_regular_kernel_n3lo(x1, log1);

  res = (res1 + res2)*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
////////////////// u-d channel ///////////////////
//////////////////////////////////////////////////

// NNLO u-d regular term
double qqprime_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1,res2;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqprime_regular_kernel_nnlo(x1, log1);
  res1 = asopi2*std::get<0>(result0);
  res2 = asopi2*std::get<1>(result0);

  res = (res1*dlumqqprime(x2,tau/x1/x2,muf2,pdf) +
	 res2*dlumqqprime2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO u-d regular term
double qqprime_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double res3, res4;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qqprime_regular_kernel_nnlo(x1, log1);
  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  res2 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  auto result1 = qqprime_regular_kernel_n3lo(x1, log1);
  res3 = asopi3*std::get<0>(result1);
  res4 = asopi3*std::get<1>(result1);

  res = ((res1 + res3)*dlumqqprime(x2,tau/x1/x2,muf2,pdf) +
	 (res2 + res4)*dlumqqprime2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
/////////////// ubar-dbar channel ////////////////
//////////////////////////////////////////////////

// NNLO ubar-dbar regular term
double qbqprimeb_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double res1,res2;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qbqprimeb_regular_kernel_nnlo(x1, log1);
  res1 = asopi2*std::get<0>(result0);
  res2 = asopi2*std::get<1>(result0);

  res = (res1*dlumqbqprimeb(x2,tau/x1/x2,muf2,pdf) +
	 res2*dlumqbqprimeb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO ubar-dbar regular term
double qbqprimeb_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double res3, res4;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  auto result0 = qbqprimeb_regular_kernel_nnlo(x1, log1);
  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  res2 = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  auto result1 = qbqprimeb_regular_kernel_n3lo(x1, log1);
  res3 = asopi3*std::get<0>(result1);
  res4 = asopi3*std::get<1>(result1);

  res = ((res1 + res3)*dlumqbqprimeb(x2,tau/x1/x2,muf2,pdf) +
	 (res2 + res4)*dlumqbqprimeb2(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
////////////////// d-s channel ///////////////////
//////////////////////////////////////////////////

// NNLO d-s regular term
double ds_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*ds_regular_kernel_nnlo(x1, log1);

  res = res*dlumds(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO d-s regular term
double ds_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*ds_regular_kernel_nnlo(x1, log1);
  res2 = asopi3*ds_regular_kernel_n3lo(x1, log1);

  res = (res1 + res2)*dlumds(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}


//////////////////////////////////////////////////
/////////////// ubar-cbar channel ////////////////
//////////////////////////////////////////////////

// NNLO ub-cb regular term
double ubcb_regular_nnlo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res = asopi2*ubcb_regular_kernel_nnlo(x1, log1);

  res = res*dlumubcb(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}

// N3LO ub-cb regular term
double ubcb_regular_n3lo(const double X[], const double s, const double xmuf, const double xmur, const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauwh, Q2;
  double x1, x2;
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double res1, res2;
  double asopi, asopi2, asopi3;
  const double MHW2 = (constants::MW+constants::MH)*(constants::MW+constants::MH);

  tauwh = MHW2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauwh));
  fac = -(1.0-2.0*eps)*tau*log(tauwh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1 = (asopi2 + 2*asopi3*constants::b0*logmu1)*ubcb_regular_kernel_nnlo(x1, log1);
  res2 = asopi3*ubcb_regular_kernel_n3lo(x1, log1);

  res = (res1 + res2)*dlumubcb(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_wh(Q2);

  return res;
}
