/**
 * \file  constants.h
 * \brief This file contains the constant inputs
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace constants
{

  // Riemann parameters:
  const double Pi = 3.1415926535897932384626433832795;
  const double Pi2 = constants::Pi*constants::Pi;
  const double Pi4 = constants::Pi2*constants::Pi2;
  const double Pi6 = constants::Pi4*constants::Pi2;
  const double Zeta3 = 1.2020569031595942853997381615114499907649862923405;
  const double Zeta5 = 1.0369277551433699263313654864570341680570809195019;

  // Mass parameters:
  const double MH = 125.0;   // H pole mass
  const double MW = 80.398;  // W pole mass
  const double MZ = 91.1876; // Z pole mass
  const double Md = 0.0017;  // down-quark MSbar mass
  const double Mu = 0.0041;  // up-quark MSbar mass
  const double Ms = 0.1;     // strange-quark MSbar mass
  const double Mc = 1.29;    // charm-quark MSbar mass
  const double Mb = 4.1;     // bottom-quark MSbar mass
  const double Mt = 172.5;   // top-quark pole mass
    
  // Coupling parameters:

  // alpha(MZ) via vev,MW,MZ:
  const double vev = 246.221;
  //const double ee2 = (1.0-MW*MW/(MZ*MZ))*4*MW*MW/(vev*vev);

  // alpha(0):
  const double ee2 = 4.0*constants::Pi/137.035999139;

  // QCD parameters:
  const double Nc = 3;
  const double nf = 5;

  const int asorder = 4;

  const double b0 = 11.0/4.0-constants::nf/6.0;
  const double b1 = 51.0/8.0-19.0*constants::nf/24.0;

  // Calculation parameters:
  const double gevtopb = 3.893793656e8;
  const double zsmall0  = 0.0769231;
  const double zlarge0  = 0.75;

}

#endif // CONSTANTS_H
