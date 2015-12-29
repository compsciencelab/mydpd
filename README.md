# MYDPD1.0.4 

30/08/05 (yes, this is very old, unmaintained, etc)
                                                                                                                           
Dissipative particle dynamics (DPD) consists of a set of N particles moving in continuous space. Each particle is defined by its position, momentum and mass m. The dynamics is specified by a set of Langevin equations very similar to the molecular dynamics equations, but where in addition to the conservative forces there are dissipative and fluctuating forces as well. This model is the natural extension of the Brownian dynamics (BD),  but while BD particles interact conserving only mass,  DPD conserves mass and momentum reproducing at larger scales not just the diffusive behavior but also hydrodynamics.

mydpd is a C++  code with periodic boundary conditions in three spatial dimensions. This code is realised under the GPL license, is simple and serial. It contains two integrators for the DPD stochastic equations, a simple DPD velocity Verlet and the stochastic Trotter integrator "DPD-Trotter" derived in the articles: 

G. De Fabritiis, M. Serrano, P. Espanol and P.V. Coveney, Efficient numerical integrators for stochastic models, Physica A 361, 429 (2005). pdf

M. Serrano, G. De Fabritiis, P. Espanol and P. V. Coveney, A stochastic Trotter integration scheme for dissipative particle dynamics, Math. Comput. Simul. 72, 190 (2006).


Scientific publications using this code should cite:

@ARTICLE {defabritiis05A,
    authors = " G. {De Fabritiis} and M. Serrano and P. Espa{\~{n}}ol  and P. V. Coveney",
    title = "Efficient numerical integrators for stochastic models",
    year = "2005",
    journal = "Physica A",
    note = "to appear" 
}

 
## COMPILING  MYDPD

Using g++, it should compile without problems by just typing
from the directory src
```
	make
```
This code  was developed and tested using g++ and Linux.
If you have problems try to use a recent (see date) Linux
distribution or fix it. It should work on SUSE 9.2, Fedora Core 2/3.


## TESTING

An input file to run a simple simulation is provided. 
The input file is in the directory test and consists in a 
list of three dimensional coordinates. 
From the test directory launch
```
	../src/mydpd 0.05 100 input.coo out > log
```

## NOTES ON RANDOM NUMBERS GENERATOR

The version furnished here is BSD implementation of the MT19937 algorithm coded by 
Takuji Nishimura and Makoto Matsumoto.  Check random.C for more information.



