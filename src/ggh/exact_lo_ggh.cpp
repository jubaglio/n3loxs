/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 28/02/2021
Exact 1-loop LO form factor for g g --> H process
See eq.(50) of arXiv:1612.07651
*********************************************************************
********************************************************************* */

#include "constants.h"

#include "exact_lo_ggh.h"

#include <complex>
#include <cmath>

#include <iostream>


// function for exact LO form factor
double oneloopfac(const double m_q)
{
  double tau_q;
  const double shat = constants::MH*constants::MH;
  double res;
  std::complex<double> aux(1.0, 0.0);
  std::complex<double> res0(1.0, 0.0);
  const std::complex<double> ipi(0.0, constants::Pi);
  

  tau_q = 4*m_q*m_q/shat;

  if(tau_q>=1.0)
    {
      aux = asin(1.0/sqrt(tau_q))*asin(1.0/sqrt(tau_q));
    }
  else
    {
      aux = -0.25*(log((1.0 + sqrt(1.0 - tau_q))/(1.0 - sqrt(1.0 - tau_q))) - ipi)*
	(log((1.0 + sqrt(1.0 - tau_q))/(1.0 - sqrt(1.0 - tau_q))) - ipi);
    }
  
  res0 = 1.5*tau_q*(1.0 + (1.0 - tau_q)*aux);
 
  res = abs(res0)*abs(res0);

  return res;
}
