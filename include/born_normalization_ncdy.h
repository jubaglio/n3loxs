/**
 *
 * \file    born_normalization_ncdy.h
 * \author  Julien Baglio
 * \date    September 2021
 *
 */

/**
 *
 * \brief   The header for the neutral-current DY Born normalization functions
 *
 */

double BW(const double);         // Breit-Wigner factor for Z propagator

double BornDYphot(const double); // Born DY factor for the gamma* contribution

double BornDYint(const double);  // Born DY factor for the gamma*-Zv* interference contribution

//double BornDYZv(const double);   // Born DY factor for the Zv* contribution

//double BornDYZa(const double);   // Born DY factor for the Zax* contribution

double BornDYZ(const double);   // Born DY factor for the Z* contribution
