/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 20/10/2021
Soft+Virtual + Regular hard terms for all Z H subprocesses 
(DY-type) up to N3LO QCD (dynamical scale)
*********************************************************************
********************************************************************* */

#include "ncdy_kernels.h"

// pdf functions
#include "pdffunctions.h"

#include "zh_dyn_functions.h"

#include "constants.h"
#include "ncdy_couplings.h"

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

double gammalo_zh(const double Q2)
{
  double lambda, result, Born;

  lambda = (1.0-constants::MZ*constants::MZ/Q2-constants::MH*constants::MH/Q2)*
    (1.0-constants::MZ*constants::MZ/Q2-constants::MH*constants::MH/Q2)
    -4*constants::MH*constants::MH*constants::MZ*constants::MZ/(Q2*Q2);

  Born = constants::gevtopb*constants::MZ*constants::MZ*constants::MZ*constants::MZ
    /(12*constants::Pi*constants::Nc*constants::vev*constants::vev*constants::vev*constants::vev);
  
  result = Born*(Q2*lambda+12.0*constants::MZ*constants::MZ)*
    sqrt(lambda)/((Q2-constants::MZ*constants::MZ)*(Q2-constants::MZ*constants::MZ));

  return result;
}


// Virtual delta(z) contribution up to N3LO
double delta_zh(const double X[], const double s, const double xmuf, const double xmur,
		const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauzh, Q2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double qucsumZv = -9.0/64.0 + 3.0*ncdycouplings::sw*ncdycouplings::sw/16.0 +
    ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/2.0;
  double qdcsumZv = 9.0/16.0 - ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2;
  double wc3;
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  //
  double z;
  double fac;
  double delqu, delqd;
  double res;
  double muf2, mur, mur2;
  double log1, log2, log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[1])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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
    wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf2))/(32.0*constants::Nc);
  }
  if(k>2) {
    asopi3 = asopi2*asopi;
    log3 = log1*log2;
    logmu2 = logmu1*logmu1;
    wc3 = (-1.0+constants::Nc*constants::Nc)*
      (-459.0-656*constants::Nc+3757*constants::Nc*constants::Nc-748*constants::Nc*constants::nf- 
       12*(27.0 + 49*constants::Nc*constants::Nc + 8*constants::Nc*constants::nf)*log(mt2/muf2) + 
       72*constants::Nc*(11*constants::Nc - 2*constants::nf)*log(mt2/muf2)*log(mt2/muf2) - 
       432*(3.0 + 4*constants::Nc*constants::Nc)*constants::Zeta3)/(4608*constants::Nc*constants::Nc);
  }

  z   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1.0-2.0*eps)*z*log(tau)*fac;

  delqu = dlumuub(z,tau/z,muf2,pdf);
  delqd = dlumddb(z,tau/z,muf2,pdf);

  if(k==0) {
    res = fac*gammalo_zh(Q2)*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd)/z;
  } else if(k==1) {
    res = fac*gammalo_zh(Q2)*
      (1.0 + 2*asopi/9.0*(-24 + 2*constants::Pi2 + 9*log1))*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd)/z;
  } else if(k==2) {
    res = fac*gammalo_zh(Q2)*
      ((1.0 + 2*(asopi+constants::b0*logmu1*asopi2)/9.0*(-24 + 2*constants::Pi2 + 9*log1) +
	asopi2/6480.0*
	(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3)))*
       ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd) +
       asopi2*(-27 + constants::Pi2 + 9*log1)/9.0*(axucsum*delqu+axdcsum*delqd) +
       asopi2*wc2*(wcu2*delqu+wcd2*delqd))/z;
  } else {
    res = fac*gammalo_zh(Q2)*
      ((1.0 + 2*(asopi+asopi2*constants::b0*logmu1 +
		 asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))/9.0*(-24 + 2*constants::Pi2 + 9*log1) +
	(asopi2 + 2*asopi3*constants::b0*logmu1)/6480.0*
	(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3)) +
	asopi3*(684831*constants::Pi4 - 187912*constants::Pi6 - 3240*log3*(33 + 32*constants::Pi2 - 4096*constants::Zeta3) + 
		270*constants::Pi2*(47045 + 195068*constants::Zeta3) - 108*log2*
		(36675 + 35200*constants::Pi2 + 1408*constants::Pi4 + 243000*constants::Zeta3) + 
		45*(-2117173 + 461868*constants::Zeta3 + 9038016*constants::Zeta3*constants::Zeta3 - 12599856*constants::Zeta5) - 
		180*log1*(267*constants::Pi4 + 8*constants::Pi2*(-5773 + 47448*constants::Zeta3) - 
			  72*(2547 + 4373*constants::Zeta3 + 39654*constants::Zeta5)))/2.09952e6)*
       ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd) +
       ((asopi2 + 2*asopi3*constants::b0*logmu1)*(-27 + constants::Pi2 + 9*log1)/9.0 +
	asopi3*(-27465 + 1898*constants::Pi2 - 5*constants::Pi4 + 12*log1*(519 + 26*constants::Pi2) +
		108*log2 + 10044*constants::Zeta3)/1296.0)*(axucsum*delqu+axdcsum*delqd) +
       asopi3*(150*constants::Pi2 - constants::Pi4 + 60*(6 + 7*constants::Zeta3 - 40*constants::Zeta5))/1296.0*(qucsumZv*delqu+qdcsumZv*delqd) +
       ((asopi2 + 2*asopi3*constants::b0*logmu1 +
	 2*asopi3/9.0*(-24 + 2*constants::Pi2 + 9*log1))*wc2 + asopi3*wc3)*(wcu2*delqu+wcd2*delqd))/z;
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

// PlusConst term, electric charge stripped out and included in dlumqqb
double PlusConst_zh(const double X[], const double s, const double xmuf, const double xmur,
		    const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauzh, Q2;
  double x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2;
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delqu, delqd;
  double res;
  double muf2, mur, mur2;
  double log1, log2, log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[1])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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
    wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf2))/(32.0*constants::Nc);
  }

  x2   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*x2*log(tau)*fac;

  delqu = dlumuub(x2,tau/x2,muf2,pdf);
  delqd = dlumddb(x2,tau/x2,muf2,pdf);

  if(k==0) {
    std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
    exit(1);
  } else if(k==1) {
    res = asopi*(8/3.0*log1*log(1.0-tau) + 16/3.0*intpow(log(1.0-tau),2)/2.0)*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd)/x2;    
  } else if(k==2) {
    res = (((asopi + asopi2*constants::b0*logmu1)*8.0/3.0*log1 +
	    asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
			 225*log2 + 3438*constants::Zeta3)
	    )*log(1.0-tau) +
	   ((asopi + asopi2*constants::b0*logmu1)*16/3.0 -
	    4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2)
	    )*intpow(log(1.0-tau),2)/2.0 +
	   4*asopi2/9.0*(-23 + 48*log1)*intpow(log(1.0-tau),3)/3.0 +
	   128*asopi2/9.0*intpow(log(1.0-tau),4)/4.0)*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd)/x2;
  } else {
    res = ((((asopi +
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
	      asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*16/3.0 -
	     (asopi2 +
	      2*asopi3*constants::b0*logmu1)*
	     4.0/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) -
	     asopi3/2430.0*
	     (45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	      1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	      20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3))
	     )*intpow(log(1.0-tau),2)/2.0 +
	    ((asopi2 +
	      2*asopi3*constants::b0*logmu1)*4.0/9.0*(-23 + 48*log1) -
	     2*asopi3/81.0*
	     (5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 +
	      756*log2 - 384*log3 - 19896*constants::Zeta3)
	     )*intpow(log(1.0-tau),3)/3.0 +
	    ((asopi2 +
	      2*asopi3*constants::b0*logmu1)*128/9.0 -
	     4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2)
	     )*intpow(log(1.0-tau),4)/4.0 +
	    160*asopi3/81.0*(-23 + 24*log1)*intpow(log(1.0-tau),5)/5.0 +
	    512*asopi3/27.0*intpow(log(1.0-tau),6)/6.0)*
	   ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd) +
	   // For csum*uec/q^2 terms
	   asopi3*(8/27.0*(9*log2 + log1*(-27 + constants::Pi2))*log(1.0-tau) +
		   16/27.0*(-27 + 9*log1 + constants::Pi2)*intpow(log(1.0-tau),2)/2.0)*
	   (axucsum*delqu+axdcsum*delqd) +
	   // For Wilson coefficient terms
	   asopi3*(8/3.0*log1*log(1.0-tau) +
		   16/3.0*intpow(log(1.0-tau),2)/2.0)*wc2*(wcu2*delqu+wcd2*delqd))/x2;
  }

  res = fac*gammalo_zh(Q2)*res;

  return res;
}

// PlusInt1 term, electric charge stripped out and included in dlumqqb
double PlusInt1_zh(const double X[], const double s, const double xmuf, const double xmur,
		   const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauzh, Q2;
  double x1,x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2;
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delqu, delqd;
  double resv, res;
  double resax, resueccsumax;
  double muf2, mur, mur2;
  double log1, log2, log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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
    wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf2))/(32.0*constants::Nc);
  }

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  delqu = ( dlumuub(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumuub(x2,tau/x2,muf2,pdf)/x2 );
  delqd = ( dlumddb(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumddb(x2,tau/x2,muf2,pdf)/x2 );

  if(k==0) {
    std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
    exit(1);
  } else if(k==1) {
    res = asopi*(8/3.0*log1/(1.0-x1) + 16/3.0*log(1.0-x1)/(1.0-x1))*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd);    
  } else if(k==2) {
    res = (((asopi + asopi2*constants::b0*logmu1)*8.0/3.0*log1 +
	    asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
			 225*log2 + 3438*constants::Zeta3)
	    )/(1.0-x1) +
	   ((asopi + asopi2*constants::b0*logmu1)*16/3.0 -
	    4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2)
	    )*log(1.0-x1)/(1.0-x1) +
	   4*asopi2/9.0*(-23 + 48*log1)*intpow(log(1.0-x1),2)/(1.0-x1) +
	   128*asopi2/9.0*intpow(log(1.0-x1),3)/(1.0-x1))*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd);
  } else {
    res = (((asopi +
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
	     asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*16/3.0 -
	    (asopi2 +
	     2*asopi3*constants::b0*logmu1)*
	    4.0/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) -
	    asopi3/2430.0*
	    (45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	     1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	     20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3))
	    )*log(1.0-x1)/(1.0-x1) +
	   ((asopi2 +
	     2*asopi3*constants::b0*logmu1)*4.0/9.0*(-23 + 48*log1) -
	    2*asopi3/81.0*
	    (5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 +
	     756*log2 - 384*log3 - 19896*constants::Zeta3)
	    )*intpow(log(1.0-x1),2)/(1.0-x1) +
	   ((asopi2 +
	     2*asopi3*constants::b0*logmu1)*128/9.0 -
	    4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2)
	    )*intpow(log(1.0-x1),3)/(1.0-x1) +
	   160*asopi3/81.0*(-23 + 24*log1)*intpow(log(1.0-x1),4)/(1.0-x1) +
	   512*asopi3/27.0*intpow(log(1.0-x1),5)/(1.0-x1))*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd) +
      // For csum*uec/q^2 terms
      asopi3*(8/27.0*(9*log2 + log1*(-27 + constants::Pi2))/(1.0-x1) +
	      16/27.0*(-27 + 9*log1 + constants::Pi2)*log(1.0-x1)/(1.0-x1))*
      (axucsum*delqu+axdcsum*delqd) +
      // For Wilson coefficient terms
      asopi3*(8/3.0*log1/(1.0-x1) +
	      16/3.0*log(1.0-x1)/(1.0-x1))*wc2*(wcu2*delqu+wcd2*delqd);
    
  }

  res = fac*gammalo_zh(Q2)*res;

  return res;
}

// PlusInt2 term, electric charge stripped out and included in dlumqqb
double PlusInt2_zh(const double X[], const double s, const double xmuf, const double xmur,
		   const double asopimz, const int k, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau, tauzh, Q2;
  double x1,x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2;
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delqu, delqd;
  double res;
  double muf2, mur, mur2;
  double log1, log2, log3;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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
    wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf2))/(32.0*constants::Nc);
  }

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau*exp(-(eps+(1.0-2.0*eps)*X[1])*log(x1));
  fac = intpow(1.0-2.0*eps,2)*x1*log(tau)*x2*log(x1)*fac;

  delqu = dlumuub(x2,tau/x2,muf2,pdf);
  delqd = dlumddb(x2,tau/x2,muf2,pdf);

  if(k==0) {
    std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
    exit(1);
  } else if(k==1) {
    res = asopi*(8/3.0*log1/(1.0-x1) + 16/3.0*log(1.0-x1)/(1.0-x1))*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd)/x2;    
  } else if(k==2) {
    res = (((asopi + asopi2*constants::b0*logmu1)*8.0/3.0*log1 +
	    asopi2/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 +
			 225*log2 + 3438*constants::Zeta3)
	    )/(1.0-x1) +
	   ((asopi + asopi2*constants::b0*logmu1)*16/3.0 -
	    4*asopi2/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2)
	    )*log(1.0-x1)/(1.0-x1) +
	   4*asopi2/9.0*(-23 + 48*log1)*intpow(log(1.0-x1),2)/(1.0-x1) +
	   128*asopi2/9.0*intpow(log(1.0-x1),3)/(1.0-x1))*
      ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd)/x2;
  } else {
    res = ((((asopi +
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
	      asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*16/3.0 -
	     (asopi2 +
	      2*asopi3*constants::b0*logmu1)*
	     4.0/27.0*(41 + 25*constants::Pi2 - 3*log1 - 48*log2) -
	     asopi3/2430.0*
	     (45160*constants::Pi2 + 1218*constants::Pi4 + 30*(-1727 + 1312*constants::Pi2)*log2 - 
	      1440*log3 + 65*(625 + 15768*constants::Zeta3) - 
	      20*log1*(-6715 + 3297*constants::Pi2 + 54720*constants::Zeta3))
	     )*log(1.0-x1)/(1.0-x1) +
	    ((asopi2 +
	      2*asopi3*constants::b0*logmu1)*4.0/9.0*(-23 + 48*log1) -
	     2*asopi3/81.0*
	     (5515 - 1679*constants::Pi2 + 33*(-103 + 48*constants::Pi2)*log1 +
	      756*log2 - 384*log3 - 19896*constants::Zeta3)
	     )*intpow(log(1.0-x1),2)/(1.0-x1) +
	    ((asopi2 +
	      2*asopi3*constants::b0*logmu1)*128/9.0 -
	     4*asopi3/81.0*(-1409 + 528*constants::Pi2 + 1264*log1 - 768*log2)
	     )*intpow(log(1.0-x1),3)/(1.0-x1) +
	    160*asopi3/81.0*(-23 + 24*log1)*intpow(log(1.0-x1),4)/(1.0-x1) +
	    512*asopi3/27.0*intpow(log(1.0-x1),5)/(1.0-x1))*
	   ((qu2Zv+axu2)*delqu+(qd2Zv+axd2)*delqd) +
	   // For csum*uec/q^2 terms
	   asopi3*(8/27.0*(9*log2 + log1*(-27 + constants::Pi2))/(1.0-x1) +
		   16/27.0*(-27 + 9*log1 + constants::Pi2)*log(1.0-x1)/(1.0-x1))*
	   (axucsum*delqu+axdcsum*delqd) +
	   // For Wilson coefficient terms
	   asopi3*(8/3.0*log1/(1.0-x1) +
		   16/3.0*log(1.0-x1)/(1.0-x1))*wc2*(wcu2*delqu+wcd2*delqd))/x2;
  }

  res = -fac*gammalo_zh(Q2)*res;

  return res;
}


//////////////////////////////////////////////////
///////////////// q-qbar channel /////////////////
//////////////////////////////////////////////////

// NLO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_zh_nlo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  //
  double fac;
  double resv, resax, res;
  double muf2, mur;
  double log1;
  double asopi;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 1, asopimz);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv = fac*gammalo_zh(Q2)*qqb_regular_kernel_nlo(x1, log1);
  resax = fac*gammalo_zh(Q2)*qqb_axial_regular_kernel_nlo(x1, log1);

  res = asopi*
    ((resv*qu2Zv+resax*axu2)*dlumuub(x2,tau/x1/x2,muf2,pdf)+
     (resv*qd2Zv+resax*axd2)*dlumddb(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}


// NNLO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_zh_nnlo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2ax = 5.0/16.0;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  //
  double fac;
  double delqu, delqd;
  double rescsum2, resuec2;
  double rescsum2ax, resuec2ax, resueccsumax;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2;
  double resv_nlo, resax_nlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv_nlo  = (asopi + asopi2*constants::b0*logmu1)*qqb_regular_kernel_nlo(x1, log1);
  resax_nlo = (asopi + asopi2*constants::b0*logmu1)*qqb_axial_regular_kernel_nlo(x1, log1);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  resuec2  = asopi2*std::get<0>(result0);
  rescsum2 = asopi2*std::get<1>(result0);
  auto result1 = qqb_axial_regular_kernel_nnlo(x1, log1);
  resuec2ax    = asopi2*std::get<0>(result1);
  rescsum2ax   = asopi2*std::get<1>(result1);
  resueccsumax = asopi2*std::get<2>(result1);

  delqu = dlumuub(x2,tau/x1/x2,muf2,pdf);
  delqd = dlumddb(x2,tau/x1/x2,muf2,pdf);

  res =
    ((resv_nlo+resuec2)*(qu2Zv*delqu+qd2Zv*delqd)+
     rescsum2*csum2Zv*(9.0/4.0*delqu+9.0*delqd)+
     (resax_nlo+resuec2ax)*(axu2*delqu+axd2*delqd)+
     rescsum2ax*csum2ax*(9.0/4.0*delqu+9.0*delqd)+
     resueccsumax*(9.0/4.0*csumQuax*delqu+9.0*csumQdax*delqd))/x1/x2;

  res = res*fac*gammalo_zh(Q2);
  
  return res;
}

// N3LO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_zh_n3lo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  double csumQuZv = (3.0+4.0*ncdycouplings::sw*ncdycouplings::sw)*(-3.0+8.0*ncdycouplings::sw*ncdycouplings::sw)/144.0;
  double csumQdZv = 1.0/16.0-ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2ax = 5.0/16.0;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  double csumsqax =  1.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2;
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delqu, delqd;
  double rescsum2_n3lo, resuec2_n3lo, resueccsum_n3lo;
  double resuec2ax_n3lo, rescsum2ax_n3lo, resueccsumax_n3lo, rescsumsqax_n3lo;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  double resv_nlo, resax_nlo, resuec2_nnlo, rescsum2_nnlo;
  double resuec2ax_nnlo, rescsum2ax_nnlo, resueccsumax_nnlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf2))/(32.0*constants::Nc);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv_nlo  = (asopi +
	       asopi2*constants::b0*logmu1 +
	       asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*qqb_regular_kernel_nlo(x1, log1);
  resax_nlo = (asopi +
	       asopi2*constants::b0*logmu1 +
	       asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*qqb_axial_regular_kernel_nlo(x1, log1);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  resuec2_nnlo  = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  rescsum2_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);
  auto result1 = qqb_axial_regular_kernel_nnlo(x1, log1);
  resuec2ax_nnlo    = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result1);
  rescsum2ax_nnlo   = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result1);
  resueccsumax_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<2>(result1);

  auto result2 = qqb_regular_kernel_n3lo(x1, log1);
  resuec2_n3lo    = asopi3*std::get<0>(result2);
  rescsum2_n3lo   = asopi3*std::get<1>(result2);
  resueccsum_n3lo = asopi3*std::get<2>(result2);
  auto result3 = qqb_axial_regular_kernel_n3lo(x1, log1);
  resuec2ax_n3lo    = asopi3*std::get<0>(result3);
  rescsum2ax_n3lo   = asopi3*std::get<1>(result3);
  resueccsumax_n3lo = asopi3*std::get<2>(result3);
  rescsumsqax_n3lo  = asopi3*std::get<3>(result3);
  reswc = asopi3*qqb_axial_regular_kernel_nlo(x1, log1);

  delqu = dlumuub(x2,tau/x1/x2,muf2,pdf);
  delqd = dlumddb(x2,tau/x1/x2,muf2,pdf);

  res =
    ((resv_nlo+resuec2_nnlo+resuec2_n3lo)*(qu2Zv*delqu+qd2Zv*delqd)+
     (rescsum2_nnlo+rescsum2_n3lo)*csum2Zv*(9.0/4.0*delqu+9.0*delqd)+
     resueccsum_n3lo*(9.0/4.0*csumQuZv*delqu+9.0*csumQdZv*delqd)+
     (resax_nlo+resuec2ax_nnlo+resuec2ax_n3lo)*(axu2*delqu+axd2*delqd)+
     (rescsum2ax_nnlo+rescsum2ax_n3lo)*csum2ax*(9.0/4.0*delqu+9.0*delqd)+
     (resueccsumax_nnlo+resueccsumax_n3lo)*(9.0/4.0*csumQuax*delqu+9.0*csumQdax*delqd)+
     rescsumsqax_n3lo*csumsqax*(9.0/4.0*delqu+9.0*delqd)+
     reswc*wc2*(wcu2*delqu+wcd2*delqd))/x1/x2;

  res = res*fac*gammalo_zh(Q2);
  
  return res;
}


//////////////////////////////////////////////////
/////////////////// g-q channel //////////////////
//////////////////////////////////////////////////

// NLO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_zh_nlo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  //
  double fac;
  double resv, resax, res;
  double muf2, mur;
  double log1;
  double asopi;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 1, asopimz);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv = fac*gammalo_zh(Q2)*gq_regular_kernel_nlo(x1, log1);
  resax = fac*gammalo_zh(Q2)*gq_axial_regular_kernel_nlo(x1, log1);

  res = asopi*
    ((resv*qu2Zv+resax*axu2)*dlumgu(x2,tau/x1/x2,muf2,pdf)+
     (resv*qd2Zv+resax*axd2)*dlumgd(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}

// NNLO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_zh_nnlo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  //
  double fac;
  double resv, res;
  double resuec2ax, resueccsumax;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2;
  double resv_nlo, resax_nlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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

  resv_nlo = fac*gammalo_zh(Q2)*(asopi + asopi2*constants::b0*logmu1)*gq_regular_kernel_nlo(x1, log1);
  resax_nlo = fac*gammalo_zh(Q2)*(asopi + asopi2*constants::b0*logmu1)*gq_axial_regular_kernel_nlo(x1, log1);

  resv = fac*gammalo_zh(Q2)*asopi2*gq_regular_kernel_nnlo(x1, log1);
  auto result1 = gq_axial_regular_kernel_nnlo(x1, log1);
  resuec2ax    = fac*gammalo_zh(Q2)*asopi2*std::get<0>(result1);
  resueccsumax = fac*gammalo_zh(Q2)*asopi2*std::get<1>(result1);

  res =
    (((resv_nlo+resv)*qu2Zv+axu2*(resax_nlo+resuec2ax)+9.0/4.0*csumQuax*resueccsumax)*dlumgu(x2,tau/x1/x2,muf2,pdf)+
     ((resv_nlo+resv)*qd2Zv+axd2*(resax_nlo+resuec2ax)+9.0*csumQdax*resueccsumax)*dlumgd(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}

// N3LO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_zh_n3lo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  double csumQuZv = (3.0+4.0*ncdycouplings::sw*ncdycouplings::sw)*(-3.0+8.0*ncdycouplings::sw*ncdycouplings::sw)/144.0;
  double csumQdZv = 1.0/16.0-ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2ax = 5.0/16.0;
  double csumsqax =  1.0/16.0;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2;
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delgu, delgd;
  double rescsum2_n3lo, resuec2_n3lo, resueccsum_n3lo;
  double resuec2ax_n3lo, resueccsumax_n3lo, rescsum2ax_n3lo, resueccsumsqax_n3lo;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1, logmu2;
  double asopi, asopi2, asopi3;
  double resv_nlo, resax_nlo;
  double resv_nnlo, resuec2ax_nnlo, resueccsumax_nnlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur2 = xmur*xmur*Q2;
  mur  = xmur*sqrt(Q2);
  logmu1 = log(mur2/muf2);
  logmu2 = logmu1*logmu1;
  asopi  = as_n3loxs(mur, 3, asopimz);
  asopi2 = asopi*asopi;
  asopi3 = asopi*asopi2;

  wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf2))/(32.0*constants::Nc);

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv_nlo = (asopi +
	      asopi2*constants::b0*logmu1 +
	      asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*gq_regular_kernel_nlo(x1, log1);
  resax_nlo = (asopi +
	       asopi2*constants::b0*logmu1 +
	       asopi3*(constants::b0*constants::b0*logmu2 + constants::b1*logmu1))*gq_axial_regular_kernel_nlo(x1, log1);

  resv_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*gq_regular_kernel_nnlo(x1, log1);
  auto result0 = gq_axial_regular_kernel_nnlo(x1, log1);
  resuec2ax_nnlo    = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  resueccsumax_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  auto result1 = gq_regular_kernel_n3lo(x1, log1);
  resuec2_n3lo    = asopi3*std::get<0>(result1);
  rescsum2_n3lo   = asopi3*std::get<1>(result1);
  resueccsum_n3lo = asopi3*std::get<2>(result1);
  auto result2 = gq_axial_regular_kernel_n3lo(x1, log1);
  resuec2ax_n3lo      = asopi3*std::get<0>(result2);
  resueccsumax_n3lo   = asopi3*std::get<1>(result2);
  rescsum2ax_n3lo     = asopi3*std::get<2>(result2);
  resueccsumsqax_n3lo = asopi3*std::get<3>(result2);
  reswc = asopi3*gq_axial_regular_kernel_nlo(x1, log1);

  delgu = dlumgu(x2,tau/x1/x2,muf2,pdf);
  delgd = dlumgd(x2,tau/x1/x2,muf2,pdf);

  res =
    ((resv_nlo+resv_nnlo+resuec2_n3lo)*(qu2Zv*delgu+qd2Zv*delgd)+
     rescsum2_n3lo*csum2Zv*(9.0/4.0*delgu+9.0*delgd)+
     resueccsum_n3lo*(9.0/4.0*csumQuZv*delgu+9.0*csumQdZv*delgd)+
     (resax_nlo+resuec2ax_nnlo+resuec2ax_n3lo)*(axu2*delgu+axd2*delgd)+
     (resueccsumax_nnlo+resueccsumax_n3lo)*(9.0/4.0*csumQuax*delgu+9.0*csumQdax*delgd)+
     (rescsum2ax_n3lo*csum2ax+resueccsumsqax_n3lo*csumsqax)*
     (9.0/4.0*delgu+9.0*delgd)+
     reswc*wc2*(wcu2*delgu+wcd2*delgd))/x1/x2;

  res = res*fac*gammalo_zh(Q2);
  
  return res;
}


//////////////////////////////////////////////////
////////////////// g-g channel ///////////////////
//////////////////////////////////////////////////

// NNLO g-g regular term
double gg_regular_zh_nnlo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  // for Za* part:
  double csum2ax = 5.0/16.0;
  //
  double fac;
  double resv, resax, res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv = fac*gammalo_zh(Q2)*gg_regular_kernel_nnlo(x1, log1);
  resax = fac*gammalo_zh(Q2)*gg_axial_regular_kernel_nnlo(x1, log1);

  res = asopi2*(resv*csum2Zv+resax*csum2ax)*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  return res;
}

// N3LO g-g regular term
double gg_regular_zh_n3lo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  double csumsqZv = 1.0/16.0 + ncdycouplings::sw*ncdycouplings::sw/6.0 +
    ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double csum2ax = 5.0/16.0;
  double csumsqax =  1.0/16.0;
  //
  double fac;
  double res;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;
  double resv_nnlo, resax_nnlo;
  double rescsum2,rescsumsq;
  double rescsum2ax,rescsumsqax;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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

  resv_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*gg_regular_kernel_nnlo(x1, log1);
  resax_nnlo = (asopi2 + 2*asopi3*constants::b0*logmu1)*gg_axial_regular_kernel_nnlo(x1, log1);

  auto result0 = gg_regular_kernel_n3lo(x1, log1);
  rescsum2 = asopi3*std::get<0>(result0);
  rescsumsq = asopi3*std::get<1>(result0);
  auto result1 = gg_axial_regular_kernel_n3lo(x1, log1);
  rescsum2ax = asopi3*std::get<0>(result1);
  rescsumsqax = asopi3*std::get<1>(result1);

  res =
    ((resv_nnlo+rescsum2)*csum2Zv+
     rescsumsq*csumsqZv+
     (resax_nnlo+rescsum2ax)*csum2ax+
     rescsumsqax*csumsqax)*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;

  res = res*fac*gammalo_zh(Q2);

  return res;
}


//////////////////////////////////////////////////
///////////////// q-q channel ////////////////////
//////////////////////////////////////////////////

// NNLO q-q regular term, electric charge stripped out and included in dlumqq
double qq_regular_zh_nnlo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  //
  double fac;
  double resv, resax, res;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  resv = fac*gammalo_zh(Q2)*qq_regular_kernel_nnlo(x1, log1);
  resax = fac*gammalo_zh(Q2)*qq_axial_regular_kernel_nnlo(x1, log1);

  res = asopi2*
    ((resv*qu2Zv+resax*axu2)*dlumuu(x2,tau/x1/x2,muf2,pdf)+
     (resv*qd2Zv+resax*axd2)*dlumdd(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}

// N3LO q-q regular term, electric charge stripped out and included in dlumqq
double qq_regular_zh_n3lo(const double X[], const double s, const double xmuf, const double xmur,
			  const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  //
  double fac;
  double resv_n3lo, res;
  double resuec2ax_n3lo, resueccsumax_n3lo;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;
  double resv_nnlo, resax_nnlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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

  resv_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*qq_regular_kernel_nnlo(x1, log1);
  resax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*qq_axial_regular_kernel_nnlo(x1, log1);

  resv_n3lo = fac*gammalo_zh(Q2)*asopi3*qq_regular_kernel_n3lo(x1, log1);
  auto result1 = qq_axial_regular_kernel_n3lo(x1, log1);
  resuec2ax_n3lo    = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result1);
  resueccsumax_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result1);

  res =
    (((resv_nnlo+resv_n3lo)*qu2Zv+(resax_nnlo+resuec2ax_n3lo)*axu2+
      resueccsumax_n3lo*9.0/4.0*csumQuax)*dlumuu(x2,tau/x1/x2,muf2,pdf)+
     ((resv_nnlo+resv_n3lo)*qd2Zv+(resax_nnlo+resuec2ax_n3lo)*axd2+
      resueccsumax_n3lo*9.0*csumQdax)*dlumdd(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}


//////////////////////////////////////////////////
////////////////// q-Q channel ///////////////////
//////////////////////////////////////////////////

// NNLO q-Q + u-d regular term, electric charge stripped out and included in dlumqQq
double qQq_regular_zh_nnlo(const double X[], const double s, const double xmuf, const double xmur,
			   const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd;
  double axuaxd = ncdycouplings::axu*ncdycouplings::axd;
  //
  double fac;
  double res;
  double res1v;
  double resqv, resudv;
  double res1ax;
  double resqax, resudax;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1v = fac*gammalo_zh(Q2)*qQq_regular_kernel_nnlo(x1, log1);
  auto result0 = ud_regular_kernel_nnlo(x1, log1);
  resqv = fac*gammalo_zh(Q2)*std::get<0>(result0);
  resudv = fac*gammalo_zh(Q2)*std::get<1>(result0);

  res1ax = fac*gammalo_zh(Q2)*qQq_axial_regular_kernel_nnlo(x1, log1);
  auto result1 = ud_axial_regular_kernel_nnlo(x1, log1);
  resqax = fac*gammalo_zh(Q2)*std::get<0>(result1);
  resudax = fac*gammalo_zh(Q2)*std::get<1>(result1);

  res = asopi2*
    ((res1v*9.0/4.0*qu2Zv+
      res1ax*9.0/4.0*axu2)*dlumuUq(x2,tau/x1/x2,muf2,pdf)+
     (res1v*9.0*qd2Zv+
      res1ax*9.0*axd2)*dlumdDq(x2,tau/x1/x2,muf2,pdf)+
     (resqv*(qu2Zv+qd2Zv)+
      resudv*quqdZv+
      resqax*(axu2+axd2)+
      resudax*axuaxd)*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}

// N3LO q-Q + u-d regular term, electric charge stripped out and included in dlumqQq
double qQq_regular_zh_n3lo(const double X[], const double s, const double xmuf, const double xmur,
			   const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  double axuaxd = ncdycouplings::axu*ncdycouplings::axd;
  //
  double fac;
  double res;
  double res1v_n3lo;
  double resqv_n3lo, resudv_n3lo;
  double resuec2ax_n3lo, resueccsumax_n3lo;
  double resqax_n3lo, resqcsumax_n3lo, resudax_n3lo;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;
  double res1v_nnlo, resqv_nnlo, resudv_nnlo;
  double res1ax_nnlo, resqax_nnlo, resudax_nnlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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

  res1v_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*qQq_regular_kernel_nnlo(x1, log1);
  auto result0 = ud_regular_kernel_nnlo(x1, log1);
  resqv_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  resudv_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  res1ax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*qQq_axial_regular_kernel_nnlo(x1, log1);
  auto result1 = ud_axial_regular_kernel_nnlo(x1, log1);
  resqax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result1);
  resudax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result1);

  res1v_n3lo = fac*gammalo_zh(Q2)*asopi3*qQq_regular_kernel_n3lo(x1, log1);
  auto result2 = ud_regular_kernel_n3lo(x1, log1);
  resqv_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result2);
  resudv_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result2);

  auto result3 = qQq_axial_regular_kernel_n3lo(x1, log1);
  resuec2ax_n3lo    = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result3);
  resueccsumax_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result3);
  auto result4 = ud_axial_regular_kernel_n3lo(x1, log1);
  resqax_n3lo     = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result4);
  resqcsumax_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result4);
  resudax_n3lo    = fac*gammalo_zh(Q2)*asopi3*std::get<2>(result4);

  res = (((res1v_nnlo+res1v_n3lo)*9.0/4.0*qu2Zv+
	  (res1ax_nnlo+resuec2ax_n3lo)*9.0/4.0*axu2+
	  resueccsumax_n3lo*9.0/4.0*csumQuax)*dlumuUq(x2,tau/x1/x2,muf2,pdf)+
	 ((res1v_nnlo+res1v_n3lo)*9.0*qd2Zv+
	  (res1ax_nnlo+resuec2ax_n3lo)*9.0*axd2+
	  resueccsumax_n3lo*9.0*csumQdax)*dlumdDq(x2,tau/x1/x2,muf2,pdf)+
	 ((resqv_nnlo+resqv_n3lo)*(qu2Zv+qd2Zv)+
	  (resudv_nnlo+resudv_n3lo)*quqdZv+
	  (resqax_nnlo+resqax_n3lo)*(axu2+axd2)+
	  resqcsumax_n3lo*(csumQuax+csumQdax)+
	  (resudax_nnlo+resudax_n3lo)*axuaxd)*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}


//////////////////////////////////////////////////
//////////////// q-Qbar channel //////////////////
//////////////////////////////////////////////////

// NNLO q-Qbar + u-dbar regular term, electric charge stripped out and included in dlumqQqb
double qQqb_regular_zh_nnlo(const double X[], const double s, const double xmuf, const double xmur,
			    const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd;
  double axuaxd = ncdycouplings::axu*ncdycouplings::axd;
  //
  double fac;
  double res;
  double res1v;
  double resqv, resudv;
  double res1ax;
  double resqax, resudax;
  double muf2, mur;
  double log1;
  double asopi, asopi2;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
  Q2 = s*tau;

  muf2 = xmuf*xmuf*Q2;
  mur  = xmur*sqrt(Q2);
  asopi  = as_n3loxs(mur, 2, asopimz);
  asopi2 = asopi*asopi;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau)*fac;

  log1 = log(Q2/muf2);

  res1v = fac*gammalo_zh(Q2)*qQqb_regular_kernel_nnlo(x1, log1);
  auto result0 = udb_regular_kernel_nnlo(x1, log1);
  resqv = fac*gammalo_zh(Q2)*std::get<0>(result0);
  resudv = fac*gammalo_zh(Q2)*std::get<1>(result0);
  res1ax = fac*gammalo_zh(Q2)*qQqb_axial_regular_kernel_nnlo(x1, log1);

  auto result1 = udb_axial_regular_kernel_nnlo(x1, log1);
  resqax = fac*gammalo_zh(Q2)*std::get<0>(result1);
  resudax = fac*gammalo_zh(Q2)*std::get<1>(result1);

  res = asopi2*
    (((res1v*9.0/4.0*qu2Zv+
       res1ax*9.0/4.0*axu2)*dlumuUqb(x2,tau/x1/x2,muf2,pdf)+
      (res1v*9.0*qd2Zv+
       res1ax*9.0*axd2)*dlumdDqb(x2,tau/x1/x2,muf2,pdf))+
     (resqv*(qu2Zv+qd2Zv)+
      resudv*quqdZv+
      resqax*(axu2+axd2)+
      resudax*axuaxd)*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}

// N3LO q-Qbar + u-dbar regular term, electric charge stripped out and included in dlumqQqb
double qQqb_regular_zh_n3lo(const double X[], const double s, const double xmuf, const double xmur,
			    const double asopimz, LHAPDF::PDF const* const pdf)
{
  double tau, tauzh, Q2;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  double axuaxd = ncdycouplings::axu*ncdycouplings::axd;
  //
  double fac;
  double res;
  double res1v_n3lo;
  double resqv_n3lo, resudv_n3lo;
  double resuec2ax_n3lo, resueccsumax_n3lo;
  double resqax_n3lo, resqcsumax_n3lo, resudax_n3lo;
  double muf2, mur, mur2;
  double log1;
  double logmu1;
  double asopi, asopi2, asopi3;
  double res1v_nnlo, resqv_nnlo, resudv_nnlo;
  double res1ax_nnlo, resqax_nnlo, resudax_nnlo;
  const double MHZ2 = (constants::MZ+constants::MH)*(constants::MZ+constants::MH);

  tauzh = MHZ2/s;
  tau = exp((eps+(1.0-2.0*eps)*X[2])*log(tauzh));
  fac = -(1.0-2.0*eps)*tau*log(tauzh);
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

  res1v_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*qQqb_regular_kernel_nnlo(x1, log1);
  auto result0 = udb_regular_kernel_nnlo(x1, log1);
  resqv_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result0);
  resudv_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result0);

  res1ax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*qQqb_axial_regular_kernel_nnlo(x1, log1);
  auto result1 = udb_axial_regular_kernel_nnlo(x1, log1);
  resqax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<0>(result1);
  resudax_nnlo = fac*gammalo_zh(Q2)*(asopi2 + 2*asopi3*constants::b0*logmu1)*std::get<1>(result1);

  res1v_n3lo = fac*gammalo_zh(Q2)*asopi3*qQqb_regular_kernel_n3lo(x1, log1);
  auto result2 = udb_regular_kernel_n3lo(x1, log1);
  resqv_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result2);
  resudv_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result2);

  auto result3 = qQqb_axial_regular_kernel_n3lo(x1, log1);
  resuec2ax_n3lo    = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result3);
  resueccsumax_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result3);

  auto result4 = udb_axial_regular_kernel_n3lo(x1, log1);
  resqax_n3lo     = fac*gammalo_zh(Q2)*asopi3*std::get<0>(result4);
  resqcsumax_n3lo = fac*gammalo_zh(Q2)*asopi3*std::get<1>(result4);
  resudax_n3lo    = fac*gammalo_zh(Q2)*asopi3*std::get<2>(result4);

  res = (((res1v_nnlo+res1v_n3lo)*9.0/4.0*qu2Zv+
	  (res1ax_nnlo+resuec2ax_n3lo)*9.0/4.0*axu2+
	  resueccsumax_n3lo*9.0/4.0*csumQuax)*dlumuUqb(x2,tau/x1/x2,muf2,pdf)+
	 ((res1v_nnlo+res1v_n3lo)*9.0*qd2Zv+
	  (res1ax_nnlo+resuec2ax_n3lo)*9.0*axd2+
	  resueccsumax_n3lo*9.0*csumQdax)*dlumdDqb(x2,tau/x1/x2,muf2,pdf)+
	 ((resqv_nnlo+resqv_n3lo)*(qu2Zv+qd2Zv)+
	  (resudv_nnlo+resudv_n3lo)*quqdZv+
	  (resqax_nnlo+resqax_n3lo)*(axu2+axd2)+
	  resqcsumax_n3lo*(csumQuax+csumQdax)+
	  (resudax_nnlo+resudax_n3lo)*axuaxd)*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;

  return res;
}
