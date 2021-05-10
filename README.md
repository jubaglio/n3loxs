# n3loxs
N3LO cross sections calculator

A tool suite to calculate up to N3LO in QCD various cross sections at hadron colliders:

* Neutral Drell-Yan p p(pbar) --> gamma* + X (--> l+ l- + X)
* Charged Drell-Yan p p(pbar) --> W+/W- + X (--> l+ nu_l / l- ~nu_l + X)
* Higgsstrahlung p p(pbar) --> W+/W- H + X 
* Higgs production via gluon fusion g g --> H in the Born-improved approximation
* Bottom-quark fusion Higgs production b bbar --> H + X in the
  five-flavor scheme

## Installation

Prerequisites:
* A C++11 compatible C++ compiler.
* Python 3.5 or higher.
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

The main executable is the Python script `n3loxs`. It accesses the
subprograms located in the directory `subprograms`. If you change the
location of the main executable, you have to update the 44th line of
the script so that it accesses correctly the subprograms:
```shell
maindir=[put here the absolute path of the directory n3loxs on your computer in single quotes]
```

## Usage

The program calculates hadronic cross sections for various
processes up to N3LO in QCD:

* Drell-Yan (DY) processes are calculated for an off-shell gauge boson
(photon and W bosons) at a given virtuality Q, that is chosen as the
central scale by default.
* The Higgs-strahlung processes (WH production) are inclusive
calculations and the central scale by default is the sum of the Higgs
and W boson masses. Only DY-type contributions are taken into acccount
for the Higgs-strahlung  processes. They are by far the dominant
contributions.
* The gluon fusion process (ggH) is an inclusive calculation, in
the so-called Born-improved high top-mass limit (Born-improved HTL),
where the calculation is performed in an effective theory in which the
top-quark is decoupled, matched to the full Standard Model. The
predictions are rescaled to the exact leading order (one-loop) result
including the top-quark mass. The code allows the user to get results
either in the pure HTL (no rescaling) or in the Born-improved HTL, and
also to chose between the on-shell (OS) scheme or the MSbar scheme for
the top-quark mass. The central scale is by default muF = muR = MH/2.
* The bottom-quark fusion pocess (bbH) is an inclusive
calculation in the five-flavor scheme (5FS) for which the central
scale by default is muF = (MH+2*mbpole)/4, muR = MH. The bottom-quark
mass used in the Yukawa coupling is taken in the MSbar scheme.

The program accepts up to 3 arguments on the command line, as
well as an optional flag:
* `-lattice lattice`: The lattice size that is used for the
integration (integer), by default taken to be lattice=1000.
* `-seed seed`: The seed that is used to initialize the
pseudo-random-number generator (integer). By default taken to be
seed=1.
* `--filename filename`: The name of the input file (see below). By
  default this is `n3loxs_parameters.in`.
* `--scale` or `--7point`: An optional flag to calculate 15 different
predictions for the renormalization scale varied between 0.5 and 2
times the central scales of the process (`--scale`
flag); or an optional flag to calculate the seven-point scale
variation around the central scales of the process.

The user can type `./n3loxs --help` or `./n3loxs -h` to display
the informations about these command-line arguments.

The program uses an input file for the physical parameters, by default
this is `n3loxs_parameters.in`. Please note that it is possible to use
a custom input file, but this file needs to have the same structure as
the default input file. The following parameters can be modified:

---
`process`

An integer to select the process to be studied. At the
moment, 1 is for neutral Drell-Yan production (offshell photon), 2 is
for charged W+ Drell-Yan production, 3 is for charged W- Drell-Yan
production, 4 is for inclusive W+ H Higgs-strahlung production, 5 is
for inclusive W- H Higgs-strahlung production, 6 is for the inclusive
5FS bbH process, 7 is for the inclusive ggH process.

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
Drell-Yan processes (only relevan for these DY processes). Default: `100.0`.

---

`muf0`

The value, in GeV, of the user-defined central factorization scale
`muF0`. If `muf0=-1`, then `muF0` is internally fixed to the default
central scale of the chosen process. Default: `-1`.

---

`xmuf`

Floating-point coefficient rescaling the factorization scale `muF`, so
that `muF = xmuf*muF0` where `muF0` stands for the central
factorization scale of the chosen process. Default: `1.0`.

---

`mur0`

The value, in GeV, of the user-defined central renormalization scale
`muR0`. If `mur0=-1`, then `muR0` is internally fixed to the default
central scale of the chosen process. Default: `-1`.

---

`xmur`

Floating-point coefficient rescaling the factorization scale `muR`, so
that `muR = xmur*muR0` where `muR0` stands for the central
renormalization scale of the chosen process. This parameter is
discarded if the flags `--scale` or `--7point` are used. Default:
`1.0`.

---

`MW`

The value, in GeV, of the W-boson mass. Default: `80.398`.

---

`GammaW`

The value, in GeV, of the W-boson total decay width. Default: `2.085`.

---

`MZ`

The value, in GeV, of the Z-boson mass. Default: `91.1876`.

---

`MH`

The value, in GeV, of the Higgs mass. Default: `125.09`.

---

`Mt`

The value, in GeV, of the OS top-quark mass. Default: `172.5`.

---

`Mb`

The value, in GeV, of the OS bottom-quark mass. Default: `4.58`.

---

`mt(mt)`

The value, in GeV, of the MSbar top-quark mass at the scale of the
MSbar top-quark mass. Default: `162.7`.

---

`mb(mb)`

The value, in GeV, of the MSbar bottom-quark mass at the scale of the
MSbar bottom-quark mass. Default: `4.18`.

---

`vev`

The value, in GeV, of vacuum expectation value. Default: `246.221`.

---

`1/alpha`

The value of the inverse of the fine-structure constant. Default: `137.035999084`.

---

`mt_scheme`

The calculation of the ggH process can be performed in the OS
scheme (`0`) or in the MSbar scheme (`1`), as far as the top quark is
concerned. Default: `1`.

---

`htl_flag`

The results for the ggH process can be given in the pure HTL  (`0`) or
for Born-improved predictions (`1`). Default: `1`.

---

The CKM parameters are hard-coded in the
include file `constants.h` located in the `include` directory at the
moment. The user can modify these parameters, however it requires a
new compilation of the code to be taken into account.

The program produces an output file containing the cross sections up
to the desired order in QCD, including the numerical error of the
integration. The factorization and renormalization scales used for the
calculation are also reported.

Please note that the PDF used in the
calculation is the same throughout the whole evaluation of the
program: When asking e.g. for an N3LO calculation, the LO, NLO, and
NNLO results are not using the corresponding PDFs, but simply the one
provided by the user in the input file `n3loxs_parameters.in`. The
evolution of the strong coupling constant and of the MSbar masses are
consistently done at the given QCD order.

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

The calculations underlying the Drell-Yan processes have
been described in the following references,

[2] C. Duhr, F. Dulat, and B. Mistlberger, "The Drell-Yan cross
section to third order in the strong coupling
constant". Phys. Rev. Lett. 125 (2020) 172001. DOI:
[10.1103/PhysRevLett.125.172001](https://dx.doi.org/10.1103/PhysRevLett.125.172001),
arXiv:[2001.07717](https://arxiv.org/abs/2001.07717).

[3] C. Duhr, F. Dulat, and B. Mistlberger, "Charged Current Drell-Yan
Production at N3LO". JHEP 11 (2020) 143. DOI:
[10.1007/JHEP11(2020)143](https://dx.doi.org/10.1007/JHEP11(2020)143),
arXiv:[2007.13313](https://arxiv.org/abs/2007.13313).

The user should refer to the following papers when using the program
for the bottom-quark fusion process,

[4] C. Duhr, F. Dulat, and B. Mistlberger, "Higgs Boson Production in
Bottom-Quark Fusion to Third Order in the Strong
Coupling". Phys. Rev. Lett. 125 (2020) 051804. DOI:
[10.1103/PhysRevLett.125.051804](https://dx.doi.org/10.1103/PhysRevLett.125.051804),
arXiv:[1904.09990](https://arxiv.org/abs/1904.09990).

[5] C. Duhr, F. Dulat, V. Hirschi, and B. Mistlberger, "Higgs
production in bottom quark fusion: matching the 4- and 5-flavour
schemes to third order in the strong coupling". JHEP 08
(2020) 017. DOI:
[10.1007/JHEP08(2020)017](https://dx.doi.org/10.1007/JHEP08(2020)017),
arXiv:[2004.04752](https://arxiv.org/abs/2004.04752).

When using the program for the gluon fusion process, at least the following
references should be cited,

[6] C. Anastasiou, C. Duhr, F. Dulat, E. Furland, T. Gehrmann,
F. Herzog, and B. Mistlberger, "Higgs boson gluon-fusion production at
threshold in N3LO QCD". Phys. Lett. B 737 (2014) 325-328. DOI:
[10.1016/j.physletb.2014.08.067](https://dx.doi.org/10.1016/j.physletb.2014.08.067),
arXiv:[1403.4616](https://arxiv.org/abs/1403.4616).

[7] B. Mistlberger, "Higgs boson production at hadron colliders at
N3LO in QCD". JHEP 05 (2018) 028. DOI:
[10.1007/JHEP05(2018)028](https://dx.doi.org/10.1007/JHEP(2018)028),
arXiv:[1802.00833](https://arxiv.org/abs/1802.00833).

but the user is reminded that many other pappers relevant for the
desired process should also be referenced (e.g. LO, NLO, NNLO
calculations for example).

## Authors

* Julien Baglio (@jubaglio)
* Claude Duhr
* Bernhard Mistlberger
* Robert Szafron
