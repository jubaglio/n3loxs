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
  const double MH = 125.09;   // H pole mass
  const double MW = 80.398;  // W pole mass
  const double MZ = 91.1876; // Z pole mass
  const double Mb = 4.18;     // bottom-quark MSbar mass
  const double Mt = 172.5;   // top-quark pole mass
  const double mbpole = 4.58; // bottom-quark pole mass
    
  // Coupling parameters:

  // alpha(MZ) via vev,MW,MZ:
  const double vev = 246.221;
  //const double ee2 = (1.0-MW*MW/(MZ*MZ))*4*MW*MW/(vev*vev);

  // alpha(0):
  const double ee2 = 4.0*constants::Pi/137.035999139;

  // QCD parameters:
  const double Nc = 3;
  const double nf = 5;

  // QCD beta function
  const double b0 = 11.0/4.0-constants::nf/6.0;
  const double b1 = 51.0/8.0-19.0*constants::nf/24.0;
  const double b2 = 2857.0/128.0 - constants::nf*(5033.0/1152.0-constants::nf*325.0/3456.0);
  const double b3 = 149753.0/1536.0 + constants::nf*constants::nf*constants::nf*1093.0/186624.0 +
    constants::nf*constants::nf*(50065.0 + 12944.0*constants::Zeta3)/41472.0 -
    constants::nf*(1078361.0 + 39048.0*constants::Zeta3)/41472.0 + 891.0*constants::Zeta3/64.0; 

  // b-quark anomalous dimension
  const double byuk0 = -1.0;
  const double byuk1 = (-404.0 + (40.0*constants::nf)/3.0)/96.0;
  const double byuk2 = (-134892.0 + constants::nf*(8864.0 + 5760.0*constants::Zeta3) +
			560.0*constants::nf*constants::nf/3.0)/6912.0;
  const double byuk3 = -4603055.0/41472.0 + constants::nf*constants::nf*
    (-2621.0 + 72.0*constants::Pi*constants::Pi*constants::Pi*constants::Pi - 10800.0*constants::Zeta3)/31104.0 +
    constants::nf*constants::nf*constants::nf*(83.0 - 144.0*constants::Zeta3)/15552.0 -
    530*constants::Zeta3/27.0 + constants::nf*
    (91723 - 264*constants::Pi*constants::Pi*constants::Pi*constants::Pi + 102576.0*constants::Zeta3 - 55200.0*constants::Zeta5)/6912.0 +
    275*constants::Zeta5/8.0;

  // Calculation parameters:
  const double gevtopb = 3.893793656e8;
  const double zsmall0  = 0.0769231;
  const double zlarge0  = 0.75;

}

#endif // CONSTANTS_H
