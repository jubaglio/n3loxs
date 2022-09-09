/**
 * \file  alphaS_simplified.cpp
 * \brief This file contains the routines to calculate alphaS(mur, N)
 *           at a given order N in QCD, using a simplified iterative solution
 */

#include "constants.h"

#include "alphaS.h"

double as_n3loxs(double mur, int N, double as_mz)
{
  const int nsteps = 10000;
  double steps;
  double res;
  double k1;

  steps = log(mur*mur/constants::MZ/constants::MZ)/(1.0*nsteps);
  
  res = as_mz;

  // RK4 method, without any threshold

  switch(N)
    {
    case 0:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*res*res*constants::b0;
	  res = res + k1;
	}
      break;
    case 1:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*res*res*
	    (constants::b0 + res*constants::b1);
	  res = res + k1;
	}
      break;
    case 2:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*res*res*
	    (constants::b0 +
	     res*constants::b1 +
	     res*res*constants::b2);
	  res = res + k1;
	}
      break;
    case 3:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*res*res*
	    (constants::b0 +
	     res*constants::b1 +
	     res*res*constants::b2 +
	     res*res*res*constants::b3);
	  res = res + k1;
	}
      break;
    }

  return res;
}
