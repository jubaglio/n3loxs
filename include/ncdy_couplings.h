/**
 * \file  ncdy_couplings.h
 * \brief This file contains the neutral DY couplings to quarks and charged leptons
 */

#ifndef NCDY_COUPLINGS
#define NCDY_COUPLINGS

namespace ncdycouplings
{
  // electric charges
  const double qu = 2.0/3.0;
  const double qd = -1.0/3.0;
  const double ql = -1.0;
  // Weak mixing angle
  extern double sw;
  extern double cw;
  // Z axial charge
  const double axu = 0.25;
  const double axd = -0.25;
  const double axl = -0.25;
  // Z vector charge
  extern double vecu;
  extern double vecd;
  extern double vecl;

  // Z decay width and QED constant
  extern double ee2;

}

#endif // NCDY_COUPLINGS
