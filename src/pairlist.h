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


#ifndef PAIRLIST_H
#define PAIRLIST_H
#include <vector>
#include "Vector.h"
#include "distance.h"

using namespace std;

typedef vector<int> IntList;
//typedef vector<Vector> PointList;
typedef vector<IntList> BoxList;

class Pairlist {
public:
	Pairlist(double _mindist, int _steps=1) : mindist(_mindist), steps(_steps), n(0), margin(_mindist,_mindist,_mindist) {}
	void compute( const  Distance& dist, const Vector* pos, int  np, int* aindex = 0);
	int findbox(Vector p, const Distance& dist ) const;
	int num_boxes() const {return xb*yb*zb;}
//	vector<PointList> pointlist;
	BoxList boxlist;
	BoxList nbrlist;	
	BoxList fullnbrlist;	
private:
	int n;
	double mindist;
	Vector margin;
	Vector lmin;
	Vector pairdist;
	int steps;
	int xb,yb,zb, xytotb;

	void set_neighbors(const Distance& dist); 

};
#endif
