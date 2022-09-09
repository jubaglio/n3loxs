/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Programming Function: 10/09/2021
Born normalization functions for the DY process q qb -> gamma* / Z -> l l (vector part)
*********************************************************************
********************************************************************* */

#include "constants.h"
#include "ncdy_couplings.h"
#include "born_normalization_ncdy.h"

double BW(const double Q2)
{
  return Q2*Q2/((Q2-constants::MZ*constants::MZ)*(Q2-constants::MZ*constants::MZ)+constants::MZ*constants::MZ*constants::GammaZ*constants::GammaZ);
} 

double BornDYphot(const double Q2)
{
  double bornprod, borndecay;
  double ql2;
  ql2 = ncdycouplings::ql*ncdycouplings::ql;
  bornprod  = constants::gevtopb*constants::Pi*ncdycouplings::ee2/(Q2*constants::Nc);
  borndecay = ncdycouplings::ee2*ql2/(12*constants::Pi*constants::Pi);
  return bornprod*borndecay;
}

double BornDYint(const double Q2)
{
  double bornprod, borndecay;
  double sw, cw, vecl, ql;
  sw = ncdycouplings::sw;
  cw = ncdycouplings::cw;
  vecl = ncdycouplings::vecl;
  ql = ncdycouplings::ql;
  bornprod  = constants::gevtopb*constants::Pi*ncdycouplings::ee2/(Q2*sw*cw*constants::Nc);
  borndecay = BW(Q2)*(1.0-constants::MZ*constants::MZ/Q2)*ncdycouplings::ee2*vecl*ql/(12*sw*cw*constants::Pi*constants::Pi);
  return bornprod*borndecay;
}

double BornDYZ(const double Q2)
{
  double bornprod, borndecay;
  double sw2, cw2, vecl2, axl2;
  sw2 = ncdycouplings::sw*ncdycouplings::sw;
  cw2 = ncdycouplings::cw*ncdycouplings::cw;
  vecl2 = ncdycouplings::vecl*ncdycouplings::vecl;
  axl2 = ncdycouplings::axl*ncdycouplings::axl;
  bornprod  = constants::gevtopb*constants::Pi*ncdycouplings::ee2/(Q2*sw2*cw2*constants::Nc);
  borndecay = BW(Q2)*ncdycouplings::ee2*(vecl2+axl2)/(12*sw2*cw2*constants::Pi*constants::Pi);
  return bornprod*borndecay;
}
