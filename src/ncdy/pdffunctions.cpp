/* *********************************************************************
*********************************************************************
Author: Julien Baglio
E-mail: julien.baglio@cern.ch
Date of Coding: 15/09/2020
pdf luminosities for the DY process p p -> gamma* -> l l up to N3LO QCD
*********************************************************************
********************************************************************* */

// PDG index: d u s c b t
//            1 2 3 4 5 6
//            g
//            21

#include <iostream>
#include <fstream>
#include <cmath>

#include "pdffunctions.h"

#include "pdfpar.h"

static const double _almost_zero = 1e-12;

// u-ubar luminosity function including electric charges
double dlumuub(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qu,qu2;
  double dlu;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dlu = 0.0;
  for (int i=2;i<6;i = i + 2)
    {
      //:  u ubar + ubar u
      dlu +=
	pdf->xfxQ2(i,x1,muf2)  * pdf->xfxQ2(-i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(i*collider,x2,muf2);
    }

  qu = 2.0/3.0;
  qu2 = qu*qu;

  res = dlu*qu2/x1/x2;
  return res;
}

// d-dbar luminosity function including electric charges
double dlumddb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qd,qd2;
  double dld;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dld = 0.0;
  for (int i=1;i<6;i = i + 2)
    {
      //:  d dbar + dbar d
      dld +=
	pdf->xfxQ2(i,x1,muf2)  * pdf->xfxQ2(-i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(i*collider,x2,muf2);
    }

  qd = -1.0/3.0;
  qd2 = qd*qd;

  res = dld*qd2/x1/x2;
  return res;
}

// q-qbar luminosity function including electric charges
double dlumqqb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = dlumuub(x1, x2, muf2, pdf) + dlumddb(x1, x2, muf2, pdf);
  return res;
}


// g-u luminosity function including electric charges
double dlumgu(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qu,qu2;
  double dlu;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dlu = 0.0;
  for (int i=2;i<6;i = i + 2)
    {
      //:  g u (and g ubar)
      dlu += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(0,x2,muf2) +
	pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( 0,x2,muf2) +
	pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
    }

  qu = 2.0/3.0;
  qu2 = qu*qu;

  res = dlu*qu2/x1/x2;
  return res;
}

// g-d luminosity function including electric charges
double dlumgd(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qd,qd2;
  double dld;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dld = 0.0;
  for (int i=1;i<6;i = i + 2)
    {
      //:  g d (and g dbar)
      dld += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(0,x2,muf2) +
	pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( 0,x2,muf2) +
	pdf->xfxQ2( 0,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
    }

  qd = -1.0/3.0;
  qd2 = qd*qd;

  res = dld*qd2/x1/x2;
  return res;
}

// g-q luminosity function including electric charges
double dlumgq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = dlumgu(x1, x2, muf2, pdf) + dlumgd(x1, x2, muf2, pdf);
  return res;
}


// g-g luminosity function
double dlumgg(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = (pdf->xfxQ2(0,x1,muf2) * pdf->xfxQ2(0,x2,muf2))/x1/x2;
  return res;
}


// u-u luminosity function including electric charges
double dlumuu(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qu,qu2;
  double dlu;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dlu = 0.0;
  for (int i=2;i<6;i = i + 2)
    {
      //:  u u (and uubar uubar)
      dlu += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
    }

  qu = 2.0/3.0;
  qu2 = qu*qu;

  res = dlu*qu2/x1/x2;
  return res;
}

// d-d luminosity function including electric charges
double dlumdd(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qd,qd2;
  double dld;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dld = 0.0;
  for (int i=1;i<6;i = i + 2)
    {
      //:  d d (and dbar dbar)
      dld += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(i*collider,x2,muf2) +
	pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
    }

  qd = -1.0/3.0;
  qd2 = qd*qd;

  res = dld*qd2/x1/x2;
  return res;
}

// q-q luminosity function including electric charges
double dlumqq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = dlumuu(x1, x2, muf2, pdf) + dlumdd(x1, x2, muf2, pdf);
  return res;
}


// u-U non-singlet luminosity function including electric charges
double dlumuUq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qu,qu2;
  double dlu;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dlu = 0.0;
  for (int i=-4;i<0;i = i + 2)
    {
      for (int j=i+2;j<0;j = j + 2)
	{
   	  //:  u U (and ubar Ubar)
	  dlu += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(j*collider,x2,muf2) +
	    pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(-j*collider,x2,muf2) +
	    pdf->xfxQ2( j,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	    pdf->xfxQ2(-j,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
   	}
    }

  qu = 2.0/3.0;
  qu2 = qu*qu;

  res = dlu*qu2/x1/x2;
  return res;
}

// d-D non-singlet luminosity function including electric charges
double dlumdDq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qd,qd2;
  double dld;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dld = 0.0;
  for (int i=-5;i<0;i = i + 2)
    {
      for (int j=i+2;j<0;j = j + 2)
	{
   	  //:  d D (and dbar Dbar)
	  dld += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(j*collider,x2,muf2) +
	    pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(-j*collider,x2,muf2) +
	    pdf->xfxQ2( j,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	    pdf->xfxQ2(-j,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
   	}
    }

  qd = -1.0/3.0;
  qd2 = qd*qd;

  res = dld*qd2/x1/x2;
  return res;
}

// q-Q non-singlet luminosity function including electric charges
double dlumqQq(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = dlumuUq(x1, x2, muf2, pdf) + dlumdDq(x1, x2, muf2, pdf);
  return res;
}

// u-d non-singlet luminosity function (without electric charges)
double dlumud(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = 0.0;
  for (int i=-5;i<0;i += 1)
    {
      for (int j=i+1;j<0;j = j + 2)
	{
   	  //:  u d (and ubar dbar)
	  res += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(j*collider,x2,muf2) +
	    pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2(-j*collider,x2,muf2) +
	    pdf->xfxQ2( j,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	    pdf->xfxQ2(-j,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2);
   	}
    }
    
  res = res/x1/x2;
  return res;
}

// u-Ubar non-singlet luminosity function including electric charges
double dlumuUqb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qu,qu2;
  double dlu;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dlu = 0.0;
  for (int i=-4;i<0;i = i + 2)
    {
      for (int j=i+2;j<0;j = j + 2)
	{
   	  //:  u Ubar (and ubar U)
	  dlu += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(-j*collider,x2,muf2) +
	    pdf->xfxQ2(-j,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	    pdf->xfxQ2( j,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2) +
	    pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( j*collider,x2,muf2);
   	}
    }

  qu = 2.0/3.0;
  qu2 = qu*qu;

  res = dlu*qu2/x1/x2;
  return res;
}

// d-Dbar non-singlet luminosity function including electric charges
double dlumdDqb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double qd,qd2;
  double dld;
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  dld = 0.0;
  for (int i=-5;i<0;i = i + 2)
    {
      for (int j=i+2;j<0;j = j + 2)
	{
   	  //:  d Dbar (and dbar D)
	  dld += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(-j*collider,x2,muf2) +
	    pdf->xfxQ2(-j,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	    pdf->xfxQ2( j,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2) +
	    pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( j*collider,x2,muf2);
   	}
    }

  qd = -1.0/3.0;
  qd2 = qd*qd;

  res = dld*qd2/x1/x2;
  return res;
}

// q-Qbar non-singlet luminosity function including electric charges
double dlumqQqb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = dlumuUqb(x1, x2, muf2, pdf) + dlumdDqb(x1, x2, muf2, pdf);
  return res;
}

// u-dbar non-singlet luminosity function (without electric charges)
double dlumudb(const double x1, const double x2, const double muf2, LHAPDF::PDF const* const pdf)
{
  double res;

  int collider = parampdf.collidertype;

  if (x1>1.0-_almost_zero or x2>1.0-_almost_zero or x1<_almost_zero or x2<_almost_zero)
    {
      return 0.0;
      std::cout << "warning edges of pdfs!!!" << std::endl;
    }

  res = 0.0;
  for (int i=-5;i<0;i += 1)
    {
      for (int j=i+1;j<0;j = j + 2)
	{
   	  //:  u dbar (and ubar d)
	  res += pdf->xfxQ2(i,x1,muf2) * pdf->xfxQ2(-j*collider,x2,muf2) +
	    pdf->xfxQ2(-j,x1,muf2) * pdf->xfxQ2( i*collider,x2,muf2) +
	    pdf->xfxQ2( j,x1,muf2) * pdf->xfxQ2(-i*collider,x2,muf2) +
	    pdf->xfxQ2(-i,x1,muf2) * pdf->xfxQ2( j*collider,x2,muf2);
   	}
    }
    
  res = res/x1/x2;
  return res;
}
