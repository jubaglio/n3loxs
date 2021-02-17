# n3loxs
N3LO cross sections calculator

A tool suite to calculate up to N3LO in QCD various cross sections at hadron colliders:

* Neutral Drell-Yan p p(pbar) --> gamma* + X (--> l+ l- + X)
* Charged Drell-Yan p p(pbar) --> W+/W- + X (--> l+ nu_l / l- ~nu_l + X)
* Higgsstrahlung p p(pbar) --> W+/W- H + X 
* Higgs production via gluon fusion g g --> H + X [in progress]
* Bottom-quark fusion Higgs production b bbar --> H + X [in progress]

## Installation

Prerequisites:
* A C++11 compatible C++ compiler.
* The LHAPDF library installed on your computer.

The code requires GNU Scientific Library (gsl) version 2.6 or
higher. The library is shipped with the code.

Compilation of the code:

This is done in several steps. First untar the archive with
```shell
$ tar xzf n3loxs.tar.gz
```
and then enter the main directory with
```shell
$ cd ./n3loxs
```
Now compile the GSL library using the shell script `makegsl.sh`:
```shell
$ ./makegsl.sh
```
Once it is done, please make sure that the LHAPDF config program is
correctly assigned in the `Makefile`. If you have a custom
installation of LHAPDF on your computer, update the `Makefile`:
```shell
LHAPDFCONFIG = [absolute path to your LHAPDF installation]/bin/lhapdf-config
```
After this check, simply enter the following command to create all the
subprograms:
```shell
$ ./make all
```

The main executable is the shell script `n3loxs`. It accesses the
subprograms located in the directory `subprograms`. If you change the
location of the main executable, you have to update the fifth line of
the script so that it accesses correctly the subprograms:
```shell
maindir=[put here the absolute path of the directory n3loxs on your computer]
```

## Usage

The program calculates hadronic cross sections for various
processes up to N3LO in QCD. Drell-Yan (DY) processes are calculated
for an off-shell gauge boson (photon and W bosons) at a given
virtuality Q, that is chosen as the central scale by default. The
Higgs-strahlung processes (WH production) are inclusive calculations
and the central scale by default is the sum of the Higgs and W boson
masses. Only DY-type contributions are taken into acccount for the
Higgs-strahlung  processes. They are by far the dominant
contributions.

The program accepts up to 3 arguments on the command line;
* `lattice`: The lattice size that is used for the integration (integer)
* `process`: An integer to chose between the various processes
available. At the moment, 1 is for neutral Drell-Yan production
(offshell photon), 2 is for charged Drell-Yan production, 3 is for
inclusive WH Higgs-strahlung production
* `--scale`: An optional flag to calculate 15 different
predictions for the renormalization scale varied between 0.5 and 2
times the default central scale of the corresponding process. If
ommited the renormalization scale is set equal to the factorization
scale. The user can type `./n3loxs --help` or `./n3loxs -h` to display
the informations about these command-line arguments.

The program uses an input file for the physical parameters,
`n3loxs_parameters.in`. The following parameters can be modified:

---
`PDFset`

Name of the PDF set to be used in the calculation. Default:
`PDF4LHC15_nnlo_mc`.

---

`PDFnum`

Member of the PDF set to be used in the calculation. Default: `0`.

---

`order`

Up to which order in QCD the calculation shall be performed. `0`
stands for LO, `1` stands for NLO, `2` stands for NNLO, `3` stands for
N3LO. Default: `3`.

---

`collider`

The calculation can be performed for a p-p (`0`, LHC-type) or a p-pbar
(`1`, Tevatron-type) collider. Default: `0`.

---

`energy`

The hadronic center-of-mass energy of the collider, in TeV. Default:
`13.0`.

---

`Q`

The value, in GeV, for the virtuality of the gauge bosons in the
Drell-Yan processes (not relevant for Higgs-strahlung
processes). Default: `100.0`.

---

`xmuf`

Floating-point coefficient rescaling the factorization scale, so that
`muf = xmuf*mu0` where `mu0` stands for the default central scale of
the chosen process. Default: `1.0`.

---

`channel`

An integer to chose between W+ (`1`) or W- (`-1`) channels for the
charged Drell-Yan and Higgs-strahlung processes. Not relevant for
neutral Drell-Yan production.

---

Some parameters (e.g. masses, coupling constants) are hard-coded in the
include file `constants.h` located in the `include` directory. The
user can modify these parameters, however it requires a new
compilation of the code to be taken into account.

The program produces an output file containing the cross sections up
to the desired order in QCD, including the numerical error of the
integration. The factorization and renormalization scales used for the
calculation are also reported.

Please note that the PDF used in the
calculation is the same throughout the whole evaluation of the
program: When asking e.g. for an N3LO calculation, the LO, NLO, and
NNLO results are not using the corresponding PDFs, but simply the one
provided by the user in the input file `n3loxs_parameters.in`. The
evolution of the strong coupling constant is always done at the 4-loop
order.

## Citation policy

The program uses a quasi-Monte-Carlo (QMC) integration as implemented
by

[1] S. Borowka, G. Heinrich, S. Jahn, S. P. Jones, M. Kerner, and
J. Schenk, "A GPU compatible quasi-Monte Carlo integrator interfaced
to pySecDec". Comp. Phys. Commun. 240 (2019) 120. DOI:
[10.1016/j.cpc.2019.02.015](https://dx.doi.org/10.1016/j.cpc.2019.02.015),
arXiv:[1811.11720](https://arxiv.org/abs/1811.11720).

Their implementation can be found at
[this link](https://github.com/mppmu/qmc/).

In addition, the calculations underlying the Drell-Yan processes have
been described in the following references,

[2] C. Duhr, F. Dulat, and B. Mistlberger, "The Drell-Yan cross
section to third order in the strong coupling
constant". Phys. Rev. Lett. 125 (2020) 172001. DOI:
[10.1103/PhysRevLett.125.172001](https://dx.doi.org/10.1103/PhysRevLett.125.172001),
arXiv:[2001:07717](https://arxiv.org/abs/2001.07717).

[3] C. Duhr, F. Dulat, and B. Mistlberger, "Charged Current Drell-Yan
Production at N3LO". JHEP 11 (2020) 143. DOI:
[10.1007/JHEP11(2020)143](https://dx.doi.org/10.1007/JHEP11(2020)143),
arXiv:[2007:13313](https://arxiv.org/abs/2007.13313).

## Authors

* Julien Baglio (@jubaglio)
* Claude Duhr
* Bernhard Mistlberger
* Robert Szafron

<!---
✨ Note: ✨ bbH and ggH are yet to be implemented.
--->
