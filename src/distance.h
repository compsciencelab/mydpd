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




#ifndef DISTANCE_H
#define DISTANCE_H

#include "Vector.h"

class Distance {
public:
	Distance(void) : lx(0.0), ly(0.0), lz(0.0) {}
  	Distance(double _lx, double _ly, double _lz, Vector _O) {set(_lx,_ly,_lz,_O);}

  	void set(double _lx, double _ly, double _lz, Vector _O) {
	lx = _lx;
	ly = _ly;  
	lz = _lz;  
	ilx = (lx? 1.0/lx : 0.0 );
	ily = (ly? 1.0/ly : 0.0 );
	ilz = (lz? 1.0/lz : 0.0 );
	O = _O;
	}

	inline Vector delta(const Vector& p1,const  Vector& p2) const {
	Vector delta = p1 - p2;
	delta.x -= lx * rint( ilx * delta.x );
	delta.y -= ly * rint( ily * delta.y );
	delta.z -= lz * rint( ilz * delta.z );
	return delta;
	}

	inline Vector delta(const Vector& p1) const {
	Vector delta = p1 - O;
	delta.x -= lx * rint( ilx * delta.x );
	delta.y -= ly * rint( ily * delta.y );
	delta.z -= lz * rint( ilz * delta.z );
	return delta;
	}
	
	inline Vector image(const Vector& p1) const {
	Vector delta = p1 - O;
	delta.x -= lx * rint( ilx * delta.x );
	delta.y -= ly * rint( ily * delta.y );
	delta.z -= lz * rint( ilz * delta.z );
	return delta + O;
	}

	inline double a() const { return lx; }
	inline double b() const { return ly; }
	inline double c() const { return lz; }
	inline Vector origin() const { return O; }
	inline double volume() const { return lx*ly*lz; } 
private:
	double lx,ly,lz; 
	double ilx,ily,ilz;
	Vector O;
};


#endif


