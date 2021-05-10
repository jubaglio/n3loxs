/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Coding: 21/02/2021
pdf luminosities for the ggH production process up to N3LO QCD
*********************************************************************
********************************************************************* */

// PDG index: d u s c b t
//            1 2 3 4 5 6
//            g
//            21

#include <iostream>
#include <fstream>
#include <cmath>

#include "pdffunctions_ggh.h"

#include "pdfpar.h"

static const double _almost_zero = 1e-12;


// g-g luminosity function
double dlumgg(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = (pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2( 0,x2,muf2))/x1/x2;
  return res;
}


// g-q luminosity function (including also g-qbar)
double dlumgq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double dl;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dl = 0.0;
  for (int i=1;i<6;i = i + 1)
    {
      //:  q-g + qbar-g 
      dl += pdf->xfxQ2( i,x1,muf2) * pdf->xfxQ2( 0,x2,muf2) +
	pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( 0,x2,muf2) +
	pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
    }

  res = dl/x1/x2;
  return res;
}


// q-qbar luminosity function
double dlumqqb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = 0.0;
  for (int i=1;i<6;i = i + 1)
    {
      //:  q qbar
      res +=
	pdf->xfxQ2( i,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2);
    }

  res = res/x1/x2;
  return res;
}


// q-q luminosity function
double dlumqq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = 0.0;
  for (int i=1;i<6;i = i + 1)
    {
      //:  q q + qbar qbar
      res +=
	pdf->xfxQ2( i,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
    }

  res = res/x1/x2;
  return res;
}

// q1-q2 luminosity function
double dlumq1q2(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = 0.0;
  for (int i=-5;i<6;i = i + 1)
    {
      for (int j=-5;j<6;j = j + 1)
	{
	  //:  q2 q2
	  if(i!=0 && j!=0 && i!=j && i!=-j)
	    {
	      res +=
		pdf->xfxQ2( i,x1,muf2) * pdf->xfxQ2( j*collider,x2,muf2);
	    }
	}
    }

  res = res/x1/x2;
  return res;
}
