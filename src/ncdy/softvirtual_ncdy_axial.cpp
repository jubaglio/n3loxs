/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 10/09/2021
Soft+Virtual functions for the DY process q qb -> Z -> l l up to N3LO QCD (axial part)
*********************************************************************
********************************************************************* */

#include <iostream>
#include <fstream>

#include "ncdy_kernels.h"
#include "constants.h"


// Virtual delta(z) function up to N3LO
std::tuple<double, double> delta_axial_kernel(const double log1, const int k)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double delterms, deltermsueccsum;
  double log2,log3;

  delterms = 0.0;
  deltermsueccsum = 0.0;
  
  switch(k)
    {
    case 0: // LO
      delterms = 1.0;
      break;
    case 1: // NLO
      delterms = 2.0/9.0*(-24 + 2*constants::Pi2 + 9*log1);
      break;
    case 2: // NNLO
      log2 = log1*log1;
      delterms = 1.0/6480.0*
	(-58095 + 3760*constants::Pi2 - 76*constants::Pi4 - 
	 60*(-9 + 64*constants::Pi2)*log2 + 23760*constants::Zeta3 + 
	 180*log1*(37 + 16*constants::Pi2 + 488*constants::Zeta3));
      // For csum*uec/q^2 terms
      deltermsueccsum = (-27 + constants::Pi2 + 9*log1)/9.0;
      break;
    case 3: // N3LO
      log2 = log1*log1;
      log3 = log1*log2;
      delterms =
	(684831*constants::Pi4 - 187912*constants::Pi6 - 3240*log3*(33 + 32*constants::Pi2 - 4096*constants::Zeta3) + 
	 270*constants::Pi2*(47045 + 195068*constants::Zeta3) - 108*log2*
	 (36675 + 35200*constants::Pi2 + 1408*constants::Pi4 + 243000*constants::Zeta3) + 
	 45*(-2117173 + 461868*constants::Zeta3 + 9038016*constants::Zeta3*constants::Zeta3 - 12599856*constants::Zeta5) - 
	 180*log1*(267*constants::Pi4 + 8*constants::Pi2*(-5773 + 47448*constants::Zeta3) - 
		   72*(2547 + 4373*constants::Zeta3 + 39654*constants::Zeta5)))/2.09952e6;
      // For csum*uec/q^2 terms
      deltermsueccsum =
	(-27465 + 1898*constants::Pi2 - 5*constants::Pi4 + 12*log1*(519 + 26*constants::Pi2) +
	 108*log2 + 10044*constants::Zeta3)/1296.0;
      break;
    }

  return std::make_tuple(delterms, deltermsueccsum);
}

////////////////////////////////////////////////////
// Soft PlusDistributions functions up to N3LO


// PlusConst function
std::tuple<double, double> PlusConst_axial_kernel(const double tau, const double log1, const int k)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double plusterms[6];
  double plustermsueccsum[2];
  double res, resueccsum;
  double log2,log3;

  res = 0.0;
  resueccsum = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]*log(1.0-tau) + plusterms[1]*intpow(log(1.0-tau),2)/2.0;
      break;
    case 2: // NNLO
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(41 + 25*constants::Pi2 - 3*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]*log(1.0-tau) +
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0;
      break;
    case 3: // N3LO
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
	plusterms[1]*intpow(log(1.0-tau),2)/2.0 +
	plusterms[2]*intpow(log(1.0-tau),3)/3.0 +
	plusterms[3]*intpow(log(1.0-tau),4)/4.0 +
	plusterms[4]*intpow(log(1.0-tau),5)/5.0 +
	plusterms[5]*intpow(log(1.0-tau),6)/6.0;

      // For csum*uec/q^2 terms
      plustermsueccsum[0] = 8.0/27.0*(9*log2 + log1*(-27 + constants::Pi2));
      plustermsueccsum[1] = 16.0/27.0*(-27 + 9*log1 + constants::Pi2);
      
      resueccsum = plustermsueccsum[0]*log(1.0-tau) +
	plustermsueccsum[1]*intpow(log(1.0-tau),2)/2.0;
      break;
    }

  return std::make_tuple(res, resueccsum);
}

// PlusInt function
std::tuple<double, double> PlusInt_axial_kernel(const double x1, const double log1, const int k)
{
/* *******************************************************************
***  Declaration of variables 
********************************************************************* */
  double plusterms[6];
  double plustermsueccsum[2];
  double res, resueccsum;
  double log2,log3;

  res = 0.0;
  resueccsum = 0.0;

  switch(k)
    {
    case 0:
      std::cout << "Error: no PlusDistributions at LO! Program will exit..." << std::endl;
      exit(1);
      break;
    case 1: // NLO
      plusterms[0] = 8.0/3.0*log1;
      plusterms[1] = 16.0/3.0;

      res = plusterms[0]/(1.0-x1) + plusterms[1]*log(1.0-x1)/(1.0-x1);
      break;
    case 2: // NNLO
      log2 = log1*log1;
      plusterms[0] = 1.0/81.0*(-932 + 138*constants::Pi2 - 6*(41 + 25*constants::Pi2)*log1 + 225*log2 + 3438*constants::Zeta3);
      plusterms[1] =  -(4.0/27.0)*(41 + 25*constants::Pi2 - 3*log1 - 48*log2);
      plusterms[2] = 4.0/9.0*(-23 + 48*log1);
      plusterms[3] = 128.0/9.0;

      res = plusterms[0]/(1.0-x1) +
	plusterms[1]*log(1.0-x1)/(1.0-x1) +
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1);
      break;
    case 3: // N3LO
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
	plusterms[2]*intpow(log(1.0-x1),2)/(1.0-x1) +
	plusterms[3]*intpow(log(1.0-x1),3)/(1.0-x1) +
	plusterms[4]*intpow(log(1.0-x1),4)/(1.0-x1) +
	plusterms[5]*intpow(log(1.0-x1),5)/(1.0-x1);

      // For csum*uec/q^2 terms
      plustermsueccsum[0] = 8.0/27.0*(9*log2 + log1*(-27 + constants::Pi2));
      plustermsueccsum[1] = 16.0/27.0*(-27 + 9*log1 + constants::Pi2);
      
      resueccsum = plustermsueccsum[0]/(1.0-x1) +
	plustermsueccsum[1]*log(1.0-x1)/(1.0-x1);
      break;
    }

  return std::make_tuple(res, resueccsum);
}
