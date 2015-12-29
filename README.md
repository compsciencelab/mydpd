# MYDPD1.0.3  - 30/08/05 (yes, this is very old, unmaintained, etc)

## INTRODUCTION                   
                                                                                                                           
This C++ code implements in three spatial dimensions and periodic boundary conditions
the mesoscopic model dissipative particle dynamics (DPD). 
For a short introduction on DPD see the reference [1] below.


Scientific publications using this code should cite:

[1] G. De Fabritiis, M. Serrano, P. Espanol and P. V. Coveney, 
    Efficient numerical integrators for stochastic models, to appear Physica A (2005).

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



