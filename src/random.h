/******************************************************************************

    Copyright (C) 2003, 2005  Gianni De Fabritiis  
    Email: g.defabritiis@ucl.ac.uk
    OPENMD web site: http://www.openmd.org

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation; either version 2 of the License, or (at your option)
    any later version.

    Any modified version of this code can be conveniently integrated into the
    main distribution (see website openmd.org). If you find more useful to
    redistribute, you have to retain this copyright notice unchanged and state
    that it is a modified version. You also have to avoid using the official
    name openmd.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
                                                                                                          
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc., 59
    Temple Place, Suite 330, Boston, MA 02111-1307 USA

*******************************************************************************/


#ifndef RANDOM_H
#define RANDOM_H


extern  void init_genrand(unsigned long s);
extern  void init_by_array(unsigned long init_key[], int key_length);
extern  double genrand_real3(); /* generates a random number on (0,1)-real-interval */

inline double uniform() {return genrand_real3();}
inline double uniform(double a, double b) {return a + (b-a)*uniform();}

inline double normal(void) {
	double V1, V2, W, Y;  
	double X1, X2;       

	do {
		V1 = 2.0 * uniform() - 1.0;
		V2 = 2.0 * uniform() - 1.0;
		W = (V1 * V1) + (V2 * V2);
	} while (W > 1.0);

	Y = sqrt((-2.0 * log(W)) / W);
	X1 = V1 * Y;
	X2 = V2 * Y;

	return(X1);
}



#endif


