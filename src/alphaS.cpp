/**
 * \file  alphaS.cpp
 * \brief This file contains the routines to calculate alphaS(mur, N)
 *           as well as mb_msbar(muR,N) at a given order N in QCD
 */

#include "constants.h"

#include "alphaS.h"

double as_n3loxs(double mur, int N, double as_mz)
{
  const int nsteps = 100000;
  double steps;
  double res;
  double k1, k2, k3, k4;

  steps = log(mur*mur/constants::MZ/constants::MZ)/(1.0*nsteps);
  
  res = as_mz;

  // RK4 method, without any threshold

  switch(N)
    {
    case 0:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*constants::b0*res*res;
	  k2 = -steps*constants::b0*(res+k1/2.0)*(res+k1/2.0);
	  k3 = -steps*constants::b0*(res+k2/2.0)*(res+k2/2.0);
	  k4 = -steps*constants::b0*(res+k3)*(res+k3);
	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	}
      break;
    case 1:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*res*res*
	    (constants::b0 +
	     res*constants::b1);
	  k2 = -steps*(res+k1/2.0)*(res+k1/2.0)*
	    (constants::b0 +
	     (res+k1/2.0)*constants::b1);
	  k3 = -steps*(res+k2/2.0)*(res+k2/2.0)*
	    (constants::b0 +
	     (res+k2/2.0)*constants::b1);
	  k4 = -steps*(res+k3)*(res+k3)*
	    (constants::b0 +
	     (res+k3)*constants::b1);
	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	}
      break;
    case 2:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = -steps*res*res*
	    (constants::b0 +
	     res*constants::b1 +
	     res*res*constants::b2);
	  k2 = -steps*(res+k1/2.0)*(res+k1/2.0)*
	    (constants::b0 +
	     (res+k1/2.0)*constants::b1 +
	     (res+k1/2.0)*(res+k1/2.0)*constants::b2);
	  k3 = -steps*(res+k2/2.0)*(res+k2/2.0)*
	    (constants::b0 +
	     (res+k2/2.0)*constants::b1 + 
	     (res+k2/2.0)*(res+k2/2.0)*constants::b2);
	  k4 = -steps*(res+k3)*(res+k3)*
	    (constants::b0 +
	     (res+k3)*constants::b1 +
	     (res+k3)*(res+k3)*constants::b2);
	  res = res + (k1+2*k2+2*k3+k4)/6.0;
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
	  k2 = -steps*(res+k1/2.0)*(res+k1/2.0)*
	    (constants::b0 +
	     (res+steps*k1/2.0)*constants::b1 +
	     (res+steps*k1/2.0)*(res+k1/2.0)*constants::b2 +
	     (res+steps*k1/2.0)*(res+k1/2.0)*(res+k1/2.0)*constants::b3);
	  k3 = -steps*(res+k2/2.0)*(res+steps*k2/2.0)*
	    (constants::b0 +
	     (res+k2/2.0)*constants::b1 +
	     (res+k2/2.0)*(res+k2/2.0)*constants::b2 +
	     (res+k2/2.0)*(res+k2/2.0)*(res+k2/2.0)*constants::b3);
	  k4 = -steps*(res+k3)*(res+k3)*
	    (constants::b0 +
	     (res+k3)*constants::b1 +
	     (res+k3)*(res+k3)*constants::b2 +
	     (res+k3)*(res+k3)*(res+k3)*constants::b3);
	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	}
      break;
    }

  return res;
}

double mb_n3loxs(double mur, int N, double mu_init, double mb_init, double asopi_init)
{
  const int nsteps = 100000;
  double steps;
  double res, res_as;
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;

  steps = log(mur*mur/mu_init/mu_init)/(1.0*nsteps);

  res_as = asopi_init;
  res    = mb_init;

  // RK4 method, without any threshold

  switch(N)
    {
    case 0:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = steps*res*res_as*constants::byuk0;
	  l1 = -steps*res_as*res_as*constants::b0;
	  k2 = steps*(res+k1/2.0)*(res_as+l1/2.0)*constants::byuk0;
	  l2 = -steps*(res_as+l1/2.0)*(res_as+l1/2.0)*constants::b0;
	  k3 = steps*(res+k2/2.0)*(res_as+l2/2.0)*constants::byuk0;
	  l3 = -steps*(res_as+l2/2.0)*(res_as+l2/2.0)*constants::b0;
	  k4 = steps*(res+k3)*(res_as+l3)*constants::byuk0;
	  l4 = -steps*(res_as+l3)*(res_as+l3)*constants::b0;

	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	  res_as = res_as + (l1+2*l2+2*l3+l4)/6.0;
	}
      break;
    case 1:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = steps*res*res_as*
	    (constants::byuk0 +
	     res_as*constants::byuk1);
	  l1 = -steps*res_as*res_as*
	    (constants::b0 +
	     res_as*constants::b1);
	  k2 = steps*(res+k1/2.0)*(res_as+l1/2.0)*
	    (constants::byuk0 +
	     (res_as+l1/2.0)*constants::byuk1);
	  l2 = -steps*(res_as+l1/2.0)*(res_as+l1/2.0)*
	    (constants::b0 +
	     (res_as+l1/2.0)*constants::b1);
	  k3 = steps*(res+k2/2.0)*(res_as+l2/2.0)*
	    (constants::byuk0 +
	     (res_as+l2/2.0)*constants::byuk1);
	  l3 = -steps*(res_as+l2/2.0)*(res_as+l2/2.0)*
	    (constants::b0 +
	     (res_as+l2/2.0)*constants::b1);
	  k4 = steps*(res+k3)*(res_as+l3)*
	    (constants::byuk0 +
	     (res_as+l3)*constants::byuk1);
	  l4 = -steps*(res_as+l3)*(res_as+l3)*
	    (constants::b0 +
	     (res_as+l3)*constants::b1);

	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	  res_as = res_as + (l1+2*l2+2*l3+l4)/6.0;
	}
      break;
    case 2:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = steps*res*res_as*
	    (constants::byuk0 +
	     res_as*constants::byuk1 +
	     res_as*res_as*constants::byuk2);
	  l1 = -steps*res_as*res_as*
	    (constants::b0 +
	     res_as*constants::b1 +
	     res_as*res_as*constants::b2);
	  k2 = steps*(res+k1/2.0)*(res_as+l1/2.0)*
	    (constants::byuk0 +
	     (res_as+l1/2.0)*constants::byuk1 +
	     (res_as+l1/2.0)*(res_as+l1/2.0)*constants::byuk2);
	  l2 = -steps*(res_as+l1/2.0)*(res_as+l1/2.0)*
	    (constants::b0 +
	     (res_as+l1/2.0)*constants::b1 +
	     (res_as+l1/2.0)*(res_as+l1/2.0)*constants::b2);
	  k3 = steps*(res+k2/2.0)*(res_as+l2/2.0)*
	    (constants::byuk0 +
	     (res_as+l2/2.0)*constants::byuk1 +
	     (res_as+l2/2.0)*(res_as+l2/2.0)*constants::byuk2);
	  l3 = -steps*(res_as+l2/2.0)*(res_as+l2/2.0)*
	    (constants::b0 +
	     (res_as+l2/2.0)*constants::b1 +
	     (res_as+l2/2.0)*(res_as+l2/2.0)*constants::b2);
	  k4 = steps*(res+k3)*(res_as+l3)*
	    (constants::byuk0 +
	     (res_as+l3)*constants::byuk1 +
	     (res_as+l3)*(res_as+l3)*constants::byuk2);
	  l4 = -steps*(res_as+l3)*(res_as+l3)*
	    (constants::b0 +
	     (res_as+l3)*constants::b1 +
	     (res_as+l3)*(res_as+l3)*constants::b2);

	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	  res_as = res_as + (l1+2*l2+2*l3+l4)/6.0;
	}
      break;
    case 3:
      for(int i=0;i < nsteps; i++)
	{
	  k1 = steps*res*res_as*
	    (constants::byuk0 +
	     res_as*constants::byuk1 +
	     res_as*res_as*constants::byuk2 +
	     res_as*res_as*res_as*constants::byuk3);
	  l1 = -steps*res_as*res_as*
	    (constants::b0 +
	     res_as*constants::b1 +
	     res_as*res_as*constants::b2 + 
	     res_as*res_as*res_as*constants::b3);
	  k2 = steps*(res+k1/2.0)*(res_as+l1/2.0)*
	    (constants::byuk0 +
	     (res_as+l1/2.0)*constants::byuk1 +
	     (res_as+l1/2.0)*(res_as+l1/2.0)*constants::byuk2 +
	     (res_as+l1/2.0)*(res_as+l1/2.0)*(res_as+l1/2.0)*constants::byuk3);
	  l2 = -steps*(res_as+l1/2.0)*(res_as+l1/2.0)*
	    (constants::b0 +
	     (res_as+l1/2.0)*constants::b1 +
	     (res_as+l1/2.0)*(res_as+l1/2.0)*constants::b2 +
	     (res_as+l1/2.0)*(res_as+l1/2.0)*(res_as+l1/2.0)*constants::b3);
	  k3 = steps*(res+k2/2.0)*(res_as+l2/2.0)*
	    (constants::byuk0 +
	     (res_as+l2/2.0)*constants::byuk1 +
	     (res_as+l2/2.0)*(res_as+l2/2.0)*constants::byuk2 +
	     (res_as+l2/2.0)*(res_as+l2/2.0)*(res_as+l2/2.0)*constants::byuk3);
	  l3 = -steps*(res_as+l2/2.0)*(res_as+l2/2.0)*
	    (constants::b0 +
	     (res_as+l2/2.0)*constants::b1 +
	     (res_as+l2/2.0)*(res_as+l2/2.0)*constants::b2 +
	     (res_as+l2/2.0)*(res_as+l2/2.0)*(res_as+l2/2.0)*constants::b3);
	  k4 = steps*(res+k3)*(res_as+l3)*
	    (constants::byuk0 +
	     (res_as+l3)*constants::byuk1 +
	     (res_as+l3)*(res_as+l3)*constants::byuk2 +
	     (res_as+l3)*(res_as+l3)*(res_as+l3)*constants::byuk3);
	  l4 = -steps*(res_as+l3)*(res_as+l3)*
	    (constants::b0 +
	     (res_as+l3)*constants::b1 +
	     (res_as+l3)*(res_as+l3)*constants::b2 +
	     (res_as+l3)*(res_as+l3)*(res_as+l3)*constants::b3);

	  res = res + (k1+2*k2+2*k3+k4)/6.0;
	  res_as = res_as + (l1+2*l2+2*l3+l4)/6.0;
	}
      break;
    }

  return res;
  
}
