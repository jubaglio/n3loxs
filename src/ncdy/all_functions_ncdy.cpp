/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 16/09/2021
Soft+Virtual + Regular hard terms for all neutral DY subprocesses up to N3LO QCD
*********************************************************************
********************************************************************* */

#include "ncdy_kernels.h"

// pdf functions
#include "pdffunctions.h"

#include "ncdy_functions.h"

#include "born_normalization_ncdy.h"

#include "constants.h"
#include "ncdy_couplings.h"

static const double eps = 1.e-12;

// proc = 0: gamma* only
//      = 1: full process (gamma* + Z*)


// Virtual delta(z) contribution up to N3LO
double delta(const double X[], const double s, const double Q2, const double muf, const int k, const int proc, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double delterms,deltermsn3lo;
  double delterms_axial, deltermsueccsum_axial;
  // for gamma* part:
  const double qucsum = 0.5;
  const double qdcsum = -1.0;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double qucsumZv = -9.0/64.0 + 3.0*ncdycouplings::sw*ncdycouplings::sw/16.0 +
    ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/2.0;
  double qdcsumZv = 9.0/16.0 - ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  double qucsumint = 3.0/16.0 + ncdycouplings::sw*ncdycouplings::sw;
  double qdcsumint = -2*ncdycouplings::sw*ncdycouplings::sw;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf/muf))/(32.0*constants::Nc);
  double wc3 = (-1.0+constants::Nc*constants::Nc)*
    (-459.0-656*constants::Nc+3757*constants::Nc*constants::Nc-748*constants::Nc*constants::nf- 
     12*(27.0 + 49*constants::Nc*constants::Nc + 8*constants::Nc*constants::nf)*log(mt2/muf/muf) + 
     72*constants::Nc*(11*constants::Nc - 2*constants::nf)*log(mt2/muf/muf)*log(mt2/muf/muf) - 
     432*(3.0 + 4*constants::Nc*constants::Nc)*constants::Zeta3)/(4608*constants::Nc*constants::Nc);
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  //
  double z;
  double fac;
  double res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  z   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*z*log(tau);

  log1 = log(Q2/muf2);

  auto result0 = delta_kernel(log1, k);
  delterms = std::get<0>(result0);
  deltermsn3lo = std::get<1>(result0);

  if(proc==1) {
    auto result1 = delta_axial_kernel(log1, k);
    delterms_axial = std::get<0>(result1);
    deltermsueccsum_axial = std::get<1>(result1);
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*fac*tau*
	((delterms+qucsum*deltermsn3lo)*dlumuub(z,tau/z,muf2,pdf)+
	 (delterms+qdcsum*deltermsn3lo)*dlumddb(z,tau/z,muf2,pdf))/z;
    }
  else
    {
      if(k==0 || k==1) {
	res = fac*tau*
	  ((BornDYphot(Q2)*(delterms+qucsum*deltermsn3lo)+
	    BornDYint(Q2)*(qu2int*delterms+qucsumint*deltermsn3lo)+
	    BornDYZ(Q2)*(qu2Zv*delterms+qucsumZv*deltermsn3lo+
			 axu2*delterms_axial+axucsum*deltermsueccsum_axial))*dlumuub(z,tau/z,muf2,pdf)+
	   (BornDYphot(Q2)*(delterms+qdcsum*deltermsn3lo)+
	    BornDYint(Q2)*(qd2int*delterms+qdcsumint*deltermsn3lo)+
	    BornDYZ(Q2)*(qd2Zv*delterms+qucsumZv*deltermsn3lo+
			 axd2*delterms_axial+axdcsum*deltermsueccsum_axial))*dlumddb(z,tau/z,muf2,pdf))/z;
      } else if(k==2) {
	res = fac*tau*
	  ((BornDYphot(Q2)*(delterms+qucsum*deltermsn3lo)+
	    BornDYint(Q2)*(qu2int*delterms+qucsumint*deltermsn3lo)+
	    BornDYZ(Q2)*(qu2Zv*delterms+qucsumZv*deltermsn3lo+
			 axu2*delterms_axial+axucsum*deltermsueccsum_axial+
			 wc2*wcu2))*dlumuub(z,tau/z,muf2,pdf)+
	   (BornDYphot(Q2)*(delterms+qdcsum*deltermsn3lo)+
	    BornDYint(Q2)*(qd2int*delterms+qdcsumint*deltermsn3lo)+
	    BornDYZ(Q2)*(qd2Zv*delterms+qucsumZv*deltermsn3lo+
			 axd2*delterms_axial+axdcsum*deltermsueccsum_axial+
			 wc2*wcd2))*dlumddb(z,tau/z,muf2,pdf))/z;
      } else if(k==3) {
	res = fac*tau*
	  ((BornDYphot(Q2)*(delterms+qucsum*deltermsn3lo)+
	    BornDYint(Q2)*(qu2int*delterms+qucsumint*deltermsn3lo)+
	    BornDYZ(Q2)*(qu2Zv*delterms+qucsumZv*deltermsn3lo+
			 axu2*delterms_axial+axucsum*deltermsueccsum_axial+
			 (wc3+wc2*(-24 + 2*constants::Pi2 + 9*log1)*2.0/9.0)*wcu2))*dlumuub(z,tau/z,muf2,pdf)+
	   (BornDYphot(Q2)*(delterms+qdcsum*deltermsn3lo)+
	    BornDYint(Q2)*(qd2int*delterms+qdcsumint*deltermsn3lo)+
	    BornDYZ(Q2)*(qd2Zv*delterms+qdcsumZv*deltermsn3lo+
			 axd2*delterms_axial+axdcsum*deltermsueccsum_axial+
			 (wc3+wc2*(-24 + 2*constants::Pi2 + 9*log1)*2.0/9.0)*wcd2))*dlumddb(z,tau/z,muf2,pdf))/z;
      }
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
double PlusConst(const double X[], const double s, const double Q2, const double muf, const int k, const int proc, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf/muf))/(32.0*constants::Nc);
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double resv, res;
  double resax, resueccsumax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x2   = exp((eps+(1-2*eps)*X[0])*log(tau));
  fac = -(1-2*eps)*x2*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*PlusConst_kernel(tau, log1, k);
  if(proc==1) {
    auto result0 = PlusConst_axial_kernel(tau, log1, k);
    resax = fac*tau*std::get<0>(result0);
    resueccsumax = fac*tau*std::get<1>(result0);
    if(k==3) {
      auto result1 = PlusConst_axial_kernel(tau,log1, 1);
      reswc = fac*tau*std::get<0>(result1);
    }
  }

  if(proc==0)
    {
      res = resv*BornDYphot(Q2)*dlumqqb(x2,tau/x2,muf2,pdf)/x2;
    }
  else
    {
      if(k==3) {
	res =
	  ((resv*
	    (BornDYphot(Q2)+
	     BornDYZ(Q2)*qu2Zv+
	     BornDYint(Q2)*qu2int)+
	    BornDYZ(Q2)*(axu2*resax+axucsum*resueccsumax+reswc*wc2*wcu2))*dlumuub(x2,tau/x2,muf2,pdf)+
	   (resv*
	    (BornDYphot(Q2)+
	     BornDYZ(Q2)*qd2Zv+
	     BornDYint(Q2)*qd2int)+
	    BornDYZ(Q2)*(axd2*resax+axdcsum*resueccsumax+reswc*wc2*wcd2))*dlumddb(x2,tau/x2,muf2,pdf))/x2;
      } else {
	res =
	  ((resv*
	    (BornDYphot(Q2)+
	     BornDYZ(Q2)*qu2Zv+
	     BornDYint(Q2)*qu2int)+
	    BornDYZ(Q2)*(axu2*resax+axucsum*resueccsumax))*dlumuub(x2,tau/x2,muf2,pdf)+
	   (resv*
	    (BornDYphot(Q2)+
	     BornDYZ(Q2)*qd2Zv+
	     BornDYint(Q2)*qd2int)+
	    BornDYZ(Q2)*(axd2*resax+axdcsum*resueccsumax))*dlumddb(x2,tau/x2,muf2,pdf))/x2;
      }
    }

  return res;
}

// PlusInt1 term, electric charge stripped out and included in dlumqqb
double PlusInt1(const double X[], const double s, const double Q2, const double muf, const int k, const int proc, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double x1,x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf/muf))/(32.0*constants::Nc);
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double resv, res;
  double resax, resueccsumax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*PlusInt_kernel(x1, log1, k);

  if(proc==1) {
    auto result0 = PlusInt_axial_kernel(x1, log1, k);
    resax = fac*tau*std::get<0>(result0);
    resueccsumax = fac*tau*std::get<1>(result0);
    if(k==3) {
      auto result1 = PlusInt_axial_kernel(x1,log1, 1);
      reswc = fac*tau*std::get<0>(result1);
    }
  }

  if(proc==0)
    {
      res = resv*BornDYphot(Q2)*( dlumqqb(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumqqb(x2,tau/x2,muf2,pdf)/x2 );
    }
  else
    {
      if(k==3) {
	res =
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qu2Zv+
	    BornDYint(Q2)*qu2int)+
	   BornDYZ(Q2)*(axu2*resax+axucsum*resueccsumax+reswc*wc2*wcu2))*
	  ( dlumuub(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumuub(x2,tau/x2,muf2,pdf)/x2 )+
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qd2Zv+
	    BornDYint(Q2)*qd2int)+
	   BornDYZ(Q2)*(axd2*resax+axdcsum*resueccsumax+reswc*wc2*wcd2))*
	  ( dlumddb(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumddb(x2,tau/x2,muf2,pdf)/x2 );
      } else {
	res =
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qu2Zv+
	    BornDYint(Q2)*qu2int)+
	   BornDYZ(Q2)*(axu2*resax+axucsum*resueccsumax))*
	  ( dlumuub(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumuub(x2,tau/x2,muf2,pdf)/x2 )+
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qd2Zv+
	    BornDYint(Q2)*qd2int)+
	   BornDYZ(Q2)*(axd2*resax+axdcsum*resueccsumax))*
	  ( dlumddb(x2,tau/x1/x2,muf2,pdf)/x1/x2 - dlumddb(x2,tau/x2,muf2,pdf)/x2 );
      }
    }

  return res;
}

// PlusInt2 term, electric charge stripped out and included in dlumqqb
double PlusInt2(const double X[], const double s, const double Q2, const double muf, const int k, const int proc, LHAPDF::PDF const* const pdf)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double tau;
  double x1,x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  const double axucsum = -9.0/64.0;
  const double axdcsum = 9.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf/muf))/(32.0*constants::Nc);
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double resv, res;
  double resax, resueccsumax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau*exp(-(eps+(1.0-2.0*eps)*X[1])*log(x1));
  fac = intpow(1.0-2.0*eps,2)*x1*log(tau)*x2*log(x1);

  log1 = log(Q2/muf2);

  resv = -fac*tau*PlusInt_kernel(x1, log1, k);
  if(proc==1) {
    auto result0 = PlusInt_axial_kernel(x1, log1, k);
    resax = -fac*tau*std::get<0>(result0);
    resueccsumax = -fac*tau*std::get<1>(result0);
    if(k==3) {
      auto result1 = PlusInt_axial_kernel(x1,log1, 1);
      reswc = -fac*tau*std::get<0>(result1);
    }
  }  

  if(proc==0)
    {
      res = resv*BornDYphot(Q2)*dlumqqb(x2,tau/x2,muf2,pdf)/x2;
    }
  else
    {
      if(k==3) {
	res =
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qu2Zv+
	    BornDYint(Q2)*qu2int)+
	   BornDYZ(Q2)*(axu2*resax+axucsum*resueccsumax+reswc*wc2*wcu2))*dlumuub(x2,tau/x2,muf2,pdf)/x2 +
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qd2Zv+
	    BornDYint(Q2)*qd2int)+
	   BornDYZ(Q2)*(axd2*resax+axdcsum*resueccsumax+reswc*wc2*wcd2))*dlumddb(x2,tau/x2,muf2,pdf)/x2;
      } else {
	res =
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qu2Zv+
	    BornDYint(Q2)*qu2int)+
	   BornDYZ(Q2)*(axu2*resax+axucsum*resueccsumax))*dlumuub(x2,tau/x2,muf2,pdf)/x2 +
	  (resv*
	   (BornDYphot(Q2)+
	    BornDYZ(Q2)*qd2Zv+
	    BornDYint(Q2)*qd2int)+
	   BornDYZ(Q2)*(axd2*resax+axdcsum*resueccsumax))*dlumddb(x2,tau/x2,muf2,pdf)/x2;
      }
    }

  return res;
}


//////////////////////////////////////////////////
///////////////// q-qbar channel /////////////////
//////////////////////////////////////////////////

// NLO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_nlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  //
  double fac;
  double resv, resax, res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*qqb_regular_kernel_nlo(x1, log1);
  if(proc==1) {
    resax = fac*tau*qqb_axial_regular_kernel_nlo(x1, log1);
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*resv*dlumqqb(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res =
	((resv*
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qu2Zv+
	   BornDYint(Q2)*qu2int)+
	  resax*axu2*BornDYZ(Q2))*dlumuub(x2,tau/x1/x2,muf2,pdf)+
	 (resv*
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qd2Zv+
	   BornDYint(Q2)*qd2int)+
	  resax*axd2*BornDYZ(Q2))*dlumddb(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}


// NNLO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_nnlo(const double X[], const double s, const double Q2,
			const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  const double csum2 = 11.0/9.0;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  double csum2int = -7.0/6.0+22.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
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
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_nnlo(x1, log1);
  resuec2  = std::get<0>(result0);
  rescsum2 = std::get<1>(result0);
  if(proc==1) {
    auto result1 = qqb_axial_regular_kernel_nnlo(x1, log1);
    resuec2ax    = std::get<0>(result1);
    rescsum2ax   = std::get<1>(result1);
    resueccsumax = std::get<2>(result1);
  }

  delqu = dlumuub(x2,tau/x1/x2,muf2,pdf);
  delqd = dlumddb(x2,tau/x1/x2,muf2,pdf);

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(resuec2*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	 rescsum2*csum2*(9.0/4.0*delqu+9.0*delqd))/x1/x2;
    }
  else
    {
      res =
	(resuec2*
	 ((BornDYphot(Q2)+
	   BornDYZ(Q2)*qu2Zv+
	   BornDYint(Q2)*qu2int)*delqu+
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qd2Zv+
	   BornDYint(Q2)*qd2int)*delqd)+
	 rescsum2*
	 ((BornDYphot(Q2)*csum2+
	   BornDYZ(Q2)*csum2Zv+
	   BornDYint(Q2)*csum2int)*(9.0/4.0*delqu+9.0*delqd))+
	 resuec2ax*BornDYZ(Q2)*(axu2*delqu+axd2*delqd)+
	 rescsum2ax*BornDYZ(Q2)*csum2ax*(9.0/4.0*delqu+9.0*delqd)+
	 resueccsumax*BornDYZ(Q2)*(9.0/4.0*csumQuax*delqu+9.0*csumQdax*delqd))/x1/x2;
    }

  res = res*fac*tau;
  
  return res;
}

// N3LO q-qbar regular term, electric charge stripped out and included in dlumqqb
double qqb_regular_n3lo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  const double csum2 = 11.0/9.0;
  const double csumQuoverQu2 = 0.5;
  const double csumQdoverQd2 = -1.0;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  double csumQuZv = (3.0+4.0*ncdycouplings::sw*ncdycouplings::sw)*(-3.0+8.0*ncdycouplings::sw*ncdycouplings::sw)/144.0;
  double csumQdZv = 1.0/16.0-ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  double csum2int = -7.0/6.0+22.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  double csumQuint = 1.0/12.0 + 4*ncdycouplings::sw*ncdycouplings::sw/9.0;
  double csumQdint = -2.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2ax = 5.0/16.0;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  double csumsqax =  1.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf/muf))/(32.0*constants::Nc);
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delqu, delqd;
  double rescsum2, resuec2, resueccsum;
  double resuec2ax, rescsum2ax, resueccsumax, rescsumsqax;
  double res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -pow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  auto result0 = qqb_regular_kernel_n3lo(x1, log1);
  resuec2    = std::get<0>(result0);
  rescsum2   = std::get<1>(result0);
  resueccsum = std::get<2>(result0);
  if(proc==1) {
    auto result1 = qqb_axial_regular_kernel_n3lo(x1, log1);
    resuec2ax    = std::get<0>(result1);
    rescsum2ax   = std::get<1>(result1);
    resueccsumax = std::get<2>(result1);
    rescsumsqax  = std::get<3>(result1);
    reswc = qqb_axial_regular_kernel_nlo(x1, log1);
  }

  delqu = dlumuub(x2,tau/x1/x2,muf2,pdf);
  delqd = dlumddb(x2,tau/x1/x2,muf2,pdf);

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(resuec2*dlumqqb(x2,tau/x1/x2,muf2,pdf) +
	 rescsum2*csum2*(9.0/4.0*delqu+9.0*delqd) +
	 resueccsum*(csumQuoverQu2*delqu+csumQdoverQd2*delqd))/x1/x2;
    }
  else
    {
      res =
	(resuec2*
	 ((BornDYphot(Q2)+
	   BornDYZ(Q2)*qu2Zv+
	   BornDYint(Q2)*qu2int)*delqu+
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qd2Zv+
	   BornDYint(Q2)*qd2int)*delqd)+
	 rescsum2*
	 ((BornDYphot(Q2)*csum2+
	   BornDYZ(Q2)*csum2Zv+
	   BornDYint(Q2)*csum2int)*(9.0/4.0*delqu+9.0*delqd))+
	 resueccsum*
	 ((BornDYphot(Q2)*csumQuoverQu2+
	   BornDYZ(Q2)*9.0/4.0*csumQuZv+
	   BornDYint(Q2)*9.0/4.0*csumQuint)*delqu+
	  (BornDYphot(Q2)*csumQdoverQd2+
	   BornDYZ(Q2)*9.0*csumQdZv+
	   BornDYint(Q2)*9.0*csumQdint)*delqd)+
	 resuec2ax*BornDYZ(Q2)*(axu2*delqu+axd2*delqd)+
	 rescsum2ax*csum2ax*BornDYZ(Q2)*(9.0/4.0*delqu+9.0*delqd)+
	 resueccsumax*BornDYZ(Q2)*(9.0/4.0*csumQuax*delqu+9.0*csumQdax*delqd)+
	 rescsumsqax*csumsqax*BornDYZ(Q2)*(9.0/4.0*delqu+9.0*delqd)+
	 reswc*wc2*BornDYZ(Q2)*(wcu2*delqu+wcd2*delqd))/x1/x2;
    }

  res = res*fac*tau;
  
  return res;
}


//////////////////////////////////////////////////
/////////////////// g-q channel //////////////////
//////////////////////////////////////////////////

// NLO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_nlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  //
  double fac;
  double resv, resax, res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*gq_regular_kernel_nlo(x1, log1);
  if(proc==1) {
    resax = fac*tau*gq_axial_regular_kernel_nlo(x1, log1);
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*resv*dlumgq(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res =
	((resv*
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qu2Zv+
	   BornDYint(Q2)*qu2int)+
	  resax*BornDYZ(Q2)*axu2)*dlumgu(x2,tau/x1/x2,muf2,pdf)+
	 (resv*
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qd2Zv+
	   BornDYint(Q2)*qd2int)+
	  resax*BornDYZ(Q2)*axd2)*dlumgd(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}

// NNLO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_nnlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  //
  double fac;
  double resv, res;
  double resuec2ax, resueccsumax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*gq_regular_kernel_nnlo(x1, log1);
  if(proc==1) {
    auto result1 = gq_axial_regular_kernel_nnlo(x1, log1);
    resuec2ax    = std::get<0>(result1);
    resuec2ax    = fac*tau*resuec2ax;
    resueccsumax = std::get<1>(result1);
    resueccsumax  = fac*tau*resueccsumax;
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*resv*dlumgq(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res =
	((resv*(BornDYphot(Q2)+
		BornDYZ(Q2)*qu2Zv+
		BornDYint(Q2)*qu2int)+
	  BornDYZ(Q2)*(axu2*resuec2ax+9.0/4.0*csumQuax*resueccsumax))*dlumgu(x2,tau/x1/x2,muf2,pdf)+
	 (resv*(BornDYphot(Q2)+
		BornDYZ(Q2)*qd2Zv+
		BornDYint(Q2)*qd2int)+
	  BornDYZ(Q2)*(axd2*resuec2ax+9.0*csumQdax*resueccsumax))*dlumgd(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}

// N3LO g-q regular term, electric charge stripped out and included in dlumgq
double gq_regular_n3lo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  const double csum2 = 11.0/9.0;
  const double csumQuoverQu2 = 0.5;
  const double csumQdoverQd2 = -1.0;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  double csumQuZv = (3.0+4.0*ncdycouplings::sw*ncdycouplings::sw)*(-3.0+8.0*ncdycouplings::sw*ncdycouplings::sw)/144.0;
  double csumQdZv = 1.0/16.0-ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  double csum2int = -7.0/6.0+22.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  double csumQuint = 1.0/12.0 + 4*ncdycouplings::sw*ncdycouplings::sw/9.0;
  double csumQdint = -2.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csum2ax = 5.0/16.0;
  double csumsqax =  1.0/16.0;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  // for the Wilson coefficient part (Za* non-decoupling top loop):
  const double mt2=constants::Mt*constants::Mt;
  double wc2 = 3*(1.0-constants::Nc*constants::Nc)*(1.0+2*log(mt2/muf/muf))/(32.0*constants::Nc);
  double wcu2 = 2*ncdycouplings::axu*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qu/ncdycouplings::qu;
  double wcd2 = 2*ncdycouplings::axd*(3*ncdycouplings::axd + 2*ncdycouplings::axu)/ncdycouplings::qd/ncdycouplings::qd;
  double reswc;
  //
  double fac;
  double delgu, delgd;
  double rescsum2, resuec2, resueccsum;
  double resuec2ax, resueccsumax, rescsum2ax, resueccsumsqax;
  double res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  auto result0 = gq_regular_kernel_n3lo(x1, log1);
  resuec2    = std::get<0>(result0);
  rescsum2   = std::get<1>(result0);
  resueccsum = std::get<2>(result0);
  if(proc==1) {
    auto result1 = gq_axial_regular_kernel_n3lo(x1, log1);
    resuec2ax      = std::get<0>(result1);
    resueccsumax   = std::get<1>(result1);
    rescsum2ax     = std::get<2>(result1);
    resueccsumsqax = std::get<3>(result1);
    reswc = gq_axial_regular_kernel_nlo(x1, log1);
  }

  delgu = dlumgu(x2,tau/x1/x2,muf2,pdf);
  delgd = dlumgd(x2,tau/x1/x2,muf2,pdf);

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(resuec2*dlumgq(x2,tau/x1/x2,muf2,pdf) +
	 rescsum2*csum2*(9.0/4.0*delgu+9.0*delgd) +
	 resueccsum*(csumQuoverQu2*delgu+csumQdoverQd2*delgd))/x1/x2;
    }
  else
    {
      res =
	(resuec2*
	 ((BornDYphot(Q2)+
	   BornDYZ(Q2)*qu2Zv+
	   BornDYint(Q2)*qu2int)*delgu+
	  (BornDYphot(Q2)+
	   BornDYZ(Q2)*qd2Zv+
	   BornDYint(Q2)*qd2int)*delgd)+
	 rescsum2*
	 ((BornDYphot(Q2)*csum2+
	   BornDYZ(Q2)*csum2Zv+
	   BornDYint(Q2)*csum2int)*(9.0/4.0*delgu+9.0*delgd))+
	 resueccsum*
	 ((BornDYphot(Q2)*csumQuoverQu2+
	   BornDYZ(Q2)*9.0/4.0*csumQuZv+
	   BornDYint(Q2)*9.0/4.0*csumQuint)*delgu+
	  (BornDYphot(Q2)*csumQdoverQd2+
	   BornDYZ(Q2)*9.0*csumQdZv+
	   BornDYint(Q2)*9.0*csumQdint)*delgd)+
	 resuec2ax*BornDYZ(Q2)*(axu2*delgu+axd2*delgd)+
	 resueccsumax*BornDYZ(Q2)*(9.0/4.0*csumQuax*delgu+9.0*csumQdax*delgd)+
	 (rescsum2ax*csum2ax+resueccsumsqax*csumsqax)*BornDYZ(Q2)*
	 (9.0/4.0*delgu+9.0*delgd)+
	 reswc*wc2*BornDYZ(Q2)*(wcu2*delgu+wcd2*delgd))/x1/x2;
    }

  res = res*fac*tau;
  
  return res;
}


//////////////////////////////////////////////////
////////////////// g-g channel ///////////////////
//////////////////////////////////////////////////

// NNLO g-g regular term
double gg_regular_nnlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  const double csum2 = 11.0/9.0;
  // for Zv* part:
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  // for gamma*-Zv* inteference part:
  double csum2int = -7.0/6.0+22.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double csum2ax = 5.0/16.0;
  //
  double fac;
  double resv, resax, res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*gg_regular_kernel_nnlo(x1, log1);
  if(proc==1) {
    resax = fac*tau*gg_axial_regular_kernel_nnlo(x1, log1);
  }

  if(proc==0)
    {
      res = resv*BornDYphot(Q2)*csum2*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res = (resv*
	     (BornDYphot(Q2)*csum2+
	      BornDYZ(Q2)*csum2Zv+
	      BornDYint(Q2)*csum2int)+
	     resax*BornDYZ(Q2)*csum2ax)*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }

  return res;
}

// N3LO g-g regular term
double gg_regular_n3lo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  double csum2 = 11.0/9.0;
  double csumsq = 1.0/9.0;
  // for Zv* part:
  double csum2Zv = 5.0/16.0-(7.0*ncdycouplings::sw*ncdycouplings::sw)/6.0+
    (11.0*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw)/9.0;
  double csumsqZv = 1.0/16.0 + ncdycouplings::sw*ncdycouplings::sw/6.0 +
    ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for gamma*-Zv* inteference part:
  double csum2int = -7.0/6.0+22.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  double csumsqint = 1.0/6.0 + 2.0*ncdycouplings::sw*ncdycouplings::sw/9.0;
  // for Za* part:
  double csum2ax = 5.0/16.0;
  double csumsqax =  1.0/16.0;
  //
  double fac;
  double res;
  double muf2;
  double log1;
  double rescsum2,rescsumsq;
  double rescsum2ax,rescsumsqax;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  auto result0 = gg_regular_kernel_n3lo(x1, log1);
  rescsum2 = std::get<0>(result0);
  rescsumsq = std::get<1>(result0);
  if(proc==1) {
    auto result1 = gg_axial_regular_kernel_n3lo(x1, log1);
    rescsum2ax = std::get<0>(result1);
    rescsumsqax = std::get<1>(result1);
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*(csum2*rescsum2+csumsq*rescsumsq)*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res =
	(rescsum2*
	 (BornDYphot(Q2)*csum2+
	  BornDYZ(Q2)*csum2Zv+
	  BornDYint(Q2)*csum2int)+
	 rescsumsq*
	 (BornDYphot(Q2)*csumsq+
	  BornDYZ(Q2)*csumsqZv+
	  BornDYint(Q2)*csumsqint)+
	 rescsum2ax*BornDYZ(Q2)*csum2ax+
	 rescsumsqax*BornDYZ(Q2)*csumsqax)*dlumgg(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }

  res = res*fac*tau;

  return res;
}


//////////////////////////////////////////////////
///////////////// q-q channel ////////////////////
//////////////////////////////////////////////////

// NNLO q-q regular term, electric charge stripped out and included in dlumqq
double qq_regular_nnlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  //
  double fac;
  double resv, resax, res;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*qq_regular_kernel_nnlo(x1, log1);
  if(proc==1) {
    resax = fac*tau*qq_axial_regular_kernel_nnlo(x1, log1);
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*resv*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res = 
	((resv*(BornDYphot(Q2)+
		BornDYZ(Q2)*qu2Zv+
		BornDYint(Q2)*qu2int)+
	  resax*BornDYZ(Q2)*axu2)*dlumuu(x2,tau/x1/x2,muf2,pdf)+
	 (resv*(BornDYphot(Q2)+
		BornDYZ(Q2)*qd2Zv+
		BornDYint(Q2)*qd2int)+
	  resax*BornDYZ(Q2)*axd2)*dlumdd(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}

// N3LO q-q regular term, electric charge stripped out and included in dlumqq
double qq_regular_n3lo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu/ncdycouplings::qu/ncdycouplings::qu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd/ncdycouplings::qd/ncdycouplings::qd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::vecu/ncdycouplings::qu;
  double qd2int = 2*ncdycouplings::vecd/ncdycouplings::qd;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu/ncdycouplings::qu/ncdycouplings::qu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd/ncdycouplings::qd/ncdycouplings::qd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  //
  double fac;
  double resv, res;
  double resuec2ax, resueccsumax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  resv = fac*tau*qq_regular_kernel_n3lo(x1, log1);
  if(proc==1) {
    auto result1 = qq_axial_regular_kernel_n3lo(x1, log1);
    resuec2ax    = std::get<0>(result1);
    resuec2ax    = fac*tau*resuec2ax;
    resueccsumax = std::get<1>(result1);
    resueccsumax = fac*tau*resueccsumax;
  }
  
  if(proc==0)
    {
      res = BornDYphot(Q2)*resv*dlumqq(x2,tau/x1/x2,muf2,pdf)/x1/x2;
    }
  else
    {
      res =
	((resv*(BornDYphot(Q2)+
		BornDYZ(Q2)*qu2Zv+
		BornDYint(Q2)*qu2int)+
	  resuec2ax*axu2*BornDYZ(Q2)+
	  resueccsumax*9.0/4.0*csumQuax*BornDYZ(Q2))*dlumuu(x2,tau/x1/x2,muf2,pdf)+
	 (resv*(BornDYphot(Q2)+
		BornDYZ(Q2)*qd2Zv+
		BornDYint(Q2)*qd2int)+
	  resuec2ax*axd2*BornDYZ(Q2)+
	  resueccsumax*9.0*csumQdax*BornDYZ(Q2))*dlumdd(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}


//////////////////////////////////////////////////
////////////////// q-Q channel ///////////////////
//////////////////////////////////////////////////

// NNLO q-Q + u-d regular term, electric charge stripped out and included in dlumqQq
double qQq_regular_nnlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  double qu2 = ncdycouplings::qu*ncdycouplings::qu;
  double qd2 = ncdycouplings::qd*ncdycouplings::qd;
  double quqd = ncdycouplings::qu*ncdycouplings::qd;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::qu*ncdycouplings::vecu;
  double qd2int = 2*ncdycouplings::qd*ncdycouplings::vecd;
  double quqdint = ncdycouplings::qu*ncdycouplings::vecd + ncdycouplings::qd*ncdycouplings::vecu;
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
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  res1v = fac*tau*qQq_regular_kernel_nnlo(x1, log1);
  auto result0 = ud_regular_kernel_nnlo(x1, log1);
  resqv = std::get<0>(result0);
  resqv = fac*tau*resqv;
  resudv = std::get<1>(result0);
  resudv = fac*tau*resudv;

  if(proc==1) {
    res1ax = fac*tau*qQq_axial_regular_kernel_nnlo(x1, log1);
    auto result1 = ud_axial_regular_kernel_nnlo(x1, log1);
    resqax = std::get<0>(result1);
    resqax = fac*tau*resqax;
    resudax = std::get<1>(result1);
    resudax = fac*tau*resudax;
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(res1v*dlumqQq(x2,tau/x1/x2,muf2,pdf)+
	 ((qu2+qd2)*resqv+quqd*resudv)*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }
  else
    {
      res = (((res1v*(BornDYphot(Q2)+
		      BornDYZ(Q2)*9.0/4.0*qu2Zv+
		      BornDYint(Q2)*9.0/4.0*qu2int)+
	       res1ax*BornDYZ(Q2)*9.0/4.0*axu2)*dlumuUq(x2,tau/x1/x2,muf2,pdf)+
	      (res1v*(BornDYphot(Q2)+
		      BornDYZ(Q2)*9.0*qd2Zv+
		      BornDYint(Q2)*9.0*qd2int)+
	       res1ax*BornDYZ(Q2)*9.0*axd2)*dlumdDq(x2,tau/x1/x2,muf2,pdf))+
	     (resqv*
	      (BornDYphot(Q2)*(qu2+qd2)+
	       BornDYZ(Q2)*(qu2Zv+qd2Zv)+
	       BornDYint(Q2)*(qu2int+qd2int))+
	      resudv*
	      (BornDYphot(Q2)*quqd+
	       BornDYZ(Q2)*quqdZv+
	       BornDYint(Q2)*quqdint)+
	      resqax*BornDYZ(Q2)*(axu2+axd2)+
	      resudax*BornDYZ(Q2)*axuaxd)*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}

// N3LO q-Q + u-d regular term, electric charge stripped out and included in dlumqQq
double qQq_regular_n3lo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  double qu2 = ncdycouplings::qu*ncdycouplings::qu;
  double qd2 = ncdycouplings::qd*ncdycouplings::qd;
  double quqd = ncdycouplings::qu*ncdycouplings::qd;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::qu*ncdycouplings::vecu;
  double qd2int = 2*ncdycouplings::qd*ncdycouplings::vecd;
  double quqdint = ncdycouplings::qu*ncdycouplings::vecd + ncdycouplings::qd*ncdycouplings::vecu;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  double axuaxd = ncdycouplings::axu*ncdycouplings::axd;
  //
  double fac;
  double res;
  double res1v;
  double resqv, resudv;
  double resuec2ax, resueccsumax;
  double resqax, resqcsumax, resudax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  res1v = fac*tau*qQq_regular_kernel_n3lo(x1, log1);
  auto result0 = ud_regular_kernel_n3lo(x1, log1);
  resqv = std::get<0>(result0);
  resqv = fac*tau*resqv;
  resudv = std::get<1>(result0);
  resudv = fac*tau*resudv;

  if(proc==1) {
    auto result1 = qQq_axial_regular_kernel_n3lo(x1, log1);
    resuec2ax    = std::get<0>(result1);
    resuec2ax    = fac*tau*resuec2ax;
    resueccsumax = std::get<1>(result1);
    resueccsumax = fac*tau*resueccsumax;
    auto result2 = ud_axial_regular_kernel_n3lo(x1, log1);
    resqax     = std::get<0>(result2);
    resqax     = fac*tau*resqax;
    resqcsumax = std::get<1>(result2);
    resqcsumax = fac*tau*resqcsumax;
    resudax    = std::get<2>(result2);
    resudax    = fac*tau*resudax;
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(res1v*dlumqQq(x2,tau/x1/x2,muf2,pdf)+
	 ((qu2+qd2)*resqv+quqd*resudv)*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }
  else
    {
      res = ((res1v*(BornDYphot(Q2)+
		     BornDYZ(Q2)*9.0/4.0*qu2Zv+
		     BornDYint(Q2)*9.0/4.0*qu2int)+
	      resuec2ax*BornDYZ(Q2)*9.0/4.0*axu2+
	      resueccsumax*BornDYZ(Q2)*9.0/4.0*csumQuax)*dlumuUq(x2,tau/x1/x2,muf2,pdf)+
	     (res1v*(BornDYphot(Q2)+
		     BornDYZ(Q2)*9.0*qd2Zv+
		     BornDYint(Q2)*9.0*qd2int)+
	      resuec2ax*BornDYZ(Q2)*9.0*axd2+
	      resueccsumax*BornDYZ(Q2)*9.0*csumQdax)*dlumdDq(x2,tau/x1/x2,muf2,pdf)+
	     (resqv*
	      (BornDYphot(Q2)*(qu2+qd2)+
	       BornDYZ(Q2)*(qu2Zv+qd2Zv)+
	       BornDYint(Q2)*(qu2int+qd2int))+
	      resudv*
	      (BornDYphot(Q2)*quqd+
	       BornDYZ(Q2)*quqdZv+
	       BornDYint(Q2)*quqdint)+
	      resqax*BornDYZ(Q2)*(axu2+axd2)+
	      resqcsumax*BornDYZ(Q2)*(csumQuax+csumQdax)+
	      resudax*BornDYZ(Q2)*axuaxd)*dlumud(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}


//////////////////////////////////////////////////
//////////////// q-Qbar channel //////////////////
//////////////////////////////////////////////////

// NNLO q-Qbar + u-dbar regular term, electric charge stripped out and included in dlumqQqb
double qQqb_regular_nnlo(const double X[], const double s, const double Q2,
		       const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  double qu2 = ncdycouplings::qu*ncdycouplings::qu;
  double qd2 = ncdycouplings::qd*ncdycouplings::qd;
  double quqd = ncdycouplings::qu*ncdycouplings::qd;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::qu*ncdycouplings::vecu;
  double qd2int = 2*ncdycouplings::qd*ncdycouplings::vecd;
  double quqdint = ncdycouplings::qu*ncdycouplings::vecd + ncdycouplings::qd*ncdycouplings::vecu;
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
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  res1v = fac*tau*qQqb_regular_kernel_nnlo(x1, log1);
  auto result0 = udb_regular_kernel_nnlo(x1, log1);
  resqv = std::get<0>(result0);
  resqv = fac*tau*resqv;
  resudv = std::get<1>(result0);
  resudv = fac*tau*resudv;
  if(proc==1) {
    res1ax = fac*tau*qQqb_axial_regular_kernel_nnlo(x1, log1);
    auto result1 = udb_axial_regular_kernel_nnlo(x1, log1);
    resqax = std::get<0>(result1);
    resqax = fac*tau*resqax;
    resudax = std::get<1>(result1);
    resudax = fac*tau*resudax;
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(res1v*dlumqQqb(x2,tau/x1/x2,muf2,pdf)+
	 ((qu2+qd2)*resqv+quqd*resudv)*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }
  else
    {
      res = (((res1v*(BornDYphot(Q2)+
		      BornDYZ(Q2)*9.0/4.0*qu2Zv+
		      BornDYint(Q2)*9.0/4.0*qu2int)+
	       res1ax*BornDYZ(Q2)*9.0/4.0*axu2)*dlumuUqb(x2,tau/x1/x2,muf2,pdf)+
	      (res1v*(BornDYphot(Q2)+
		      BornDYZ(Q2)*9.0*qd2Zv+
		      BornDYint(Q2)*9.0*qd2int)+
	       res1ax*BornDYZ(Q2)*9.0*axd2)*dlumdDqb(x2,tau/x1/x2,muf2,pdf))+
	     (resqv*
	      (BornDYphot(Q2)*(qu2+qd2)+
	       BornDYZ(Q2)*(qu2Zv+qd2Zv)+
	       BornDYint(Q2)*(qu2int+qd2int))+
	      resudv*
	      (BornDYphot(Q2)*quqd+
	       BornDYZ(Q2)*quqdZv+
	       BornDYint(Q2)*quqdint)+
	      resqax*BornDYZ(Q2)*(axu2+axd2)+
	      resudax*BornDYZ(Q2)*axuaxd)*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}

// N3LO q-Qbar + u-dbar regular term, electric charge stripped out and included in dlumqQqb
double qQqb_regular_n3lo(const double X[], const double s, const double Q2,
			 const double muf, const int proc, LHAPDF::PDF const* const pdf)
{
  double tau;
  double x1, x2;
  // for gamma* part:
  double qu2 = ncdycouplings::qu*ncdycouplings::qu;
  double qd2 = ncdycouplings::qd*ncdycouplings::qd;
  double quqd = ncdycouplings::qu*ncdycouplings::qd;
  // for Zv* part:
  double qu2Zv = ncdycouplings::vecu*ncdycouplings::vecu;
  double qd2Zv = ncdycouplings::vecd*ncdycouplings::vecd;
  double quqdZv = ncdycouplings::vecu*ncdycouplings::vecd;
  // for gamma*-Zv* inteference part:
  double qu2int = 2*ncdycouplings::qu*ncdycouplings::vecu;
  double qd2int = 2*ncdycouplings::qd*ncdycouplings::vecd;
  double quqdint = ncdycouplings::qu*ncdycouplings::vecd + ncdycouplings::qd*ncdycouplings::vecu;
  // for Za* part:
  double axu2 = ncdycouplings::axu*ncdycouplings::axu;
  double axd2 = ncdycouplings::axd*ncdycouplings::axd;
  double csumQuax = -1.0/16.0;
  double csumQdax =  1.0/16.0;
  double axuaxd = ncdycouplings::axu*ncdycouplings::axd;
  //
  double fac;
  double res;
  double res1v;
  double resqv, resudv;
  double resuec2ax, resueccsumax;
  double resqax, resqcsumax, resudax;
  double muf2;
  double log1;

  tau = Q2/s;
  muf2 = muf*muf;

  x1 = exp((eps+(1.0-2.0*eps)*X[0])*log(tau));
  x2 = tau/x1 + (1.0-tau/x1)*(eps+(1.0-2.0*eps)*X[1]);
  fac = -intpow(1.0-2.0*eps,2)*x1*(1.0-tau/x1)*log(tau);

  log1 = log(Q2/muf2);

  res1v = fac*tau*qQqb_regular_kernel_n3lo(x1, log1);
  auto result0 = udb_regular_kernel_n3lo(x1, log1);
  resqv = std::get<0>(result0);
  resqv = fac*tau*resqv;
  resudv = std::get<1>(result0);
  resudv = fac*tau*resudv;
  if(proc==1) {
    auto result1 = qQqb_axial_regular_kernel_n3lo(x1, log1);
    resuec2ax    = std::get<0>(result1);
    resuec2ax    = fac*tau*resuec2ax;
    resueccsumax = std::get<1>(result1);
    resueccsumax = fac*tau*resueccsumax;
    auto result2 = udb_axial_regular_kernel_n3lo(x1, log1);
    resqax     = std::get<0>(result2);
    resqax     = fac*tau*resqax;
    resqcsumax = std::get<1>(result2);
    resqcsumax = fac*tau*resqcsumax;
    resudax    = std::get<2>(result2);
    resudax    = fac*tau*resudax;
  }

  if(proc==0)
    {
      res = BornDYphot(Q2)*
	(res1v*dlumqQqb(x2,tau/x1/x2,muf2,pdf)+
	 ((qu2+qd2)*resqv+quqd*resudv)*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }
  else
    {
      res = ((res1v*(BornDYphot(Q2)+
		     BornDYZ(Q2)*9.0/4.0*qu2Zv+
		     BornDYint(Q2)*9.0/4.0*qu2int)+
	      resuec2ax*BornDYZ(Q2)*9.0/4.0*axu2+
	      resueccsumax*BornDYZ(Q2)*9.0/4.0*csumQuax)*dlumuUqb(x2,tau/x1/x2,muf2,pdf)+
	     (res1v*(BornDYphot(Q2)+
		     BornDYZ(Q2)*9.0*qd2Zv+
		     BornDYint(Q2)*9.0*qd2int)+
	      resuec2ax*BornDYZ(Q2)*9.0*axd2+
	      resueccsumax*BornDYZ(Q2)*9.0*csumQdax)*dlumdDqb(x2,tau/x1/x2,muf2,pdf)+
	     (resqv*
	      (BornDYphot(Q2)*(qu2+qd2)+
	       BornDYZ(Q2)*(qu2Zv+qd2Zv)+
	       BornDYint(Q2)*(qu2int+qd2int))+
	      resudv*
	      (BornDYphot(Q2)*quqd+
	       BornDYZ(Q2)*quqdZv+
	       BornDYint(Q2)*quqdint)+
	      resqax*BornDYZ(Q2)*(axu2+axd2)+
	      resqcsumax*BornDYZ(Q2)*(csumQuax+csumQdax)+
	      resudax*BornDYZ(Q2)*axuaxd)*dlumudb(x2,tau/x1/x2,muf2,pdf))/x1/x2;
    }

  return res;
}
