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


#ifndef TIMER_H
#define TIMER_H

#include <time.h>	

class Timer {
	clock_t t0;
	clock_t final;
public:
	void tic() {t0=clock();}
	void toc() {final+=clock()-t0;}
	double time() {return double(final)/double(CLOCKS_PER_SEC);}
	Timer() : t0(0), final(0) {}
};
inline	void print_time() {
		time_t rawtime;
		struct tm * timeinfo;
	 	time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		cout<<endl<<asctime(timeinfo)<<endl;
	}
	
#endif

