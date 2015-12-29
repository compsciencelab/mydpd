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


#ifndef VECTOR_H
#define VECTOR_H

#include <iostream.h>
#include <fstream.h>
#include <math.h>


class vector3 {
public:
      typedef double DataType;
      vector3(DataType ix, DataType iy, DataType iz) { x=ix; y=iy; z=iz;} 
      vector3()  {x=0;y=0;z=0; }

      const vector3& operator+=(const vector3& v)  {x+=v.x;y+=v.y;z+=v.z;return *this;}
      const vector3& operator-=(const vector3& v)  {x-=v.x;y-=v.y;z-=v.z;return *this;}
      const vector3& operator*=(const DataType a) {x*=a;y*=a;z*=a;return *this;}
      const vector3& operator/=(const DataType a) {return *this *= 1.0/a;}

      inline double length() const {return sqrt(x*x+y*y+z*z);}
      double length2() const {return x*x + y*y + z*z;}

      DataType	x,y,z;
};


typedef vector3 Vector;

template <class InputStreamT>
inline InputStreamT& operator>>(InputStreamT& s,vector3& v) {
   s >> v.x >> v.y >> v.z;
   return s;
}

template <class OutputStreamT>
inline OutputStreamT& operator<<(OutputStreamT& s, const vector3& v) {
  s <<  v.x << ' ' << v.y << ' '<<v.z<<' ';
  return s;
}

inline vector3 operator+(const vector3& u, const vector3& v) 
	{return vector3(u.x + v.x, u.y + v.y,u.z + v.z);}
inline vector3 operator-(const vector3& v) 
	{return vector3(-v.x, -v.y, -v.z);}
inline vector3::DataType operator*(const vector3& u, const vector3& v) 
	{return u.x * v.x + u.y * v.y + u.z*v.z;}
inline vector3 operator*(const vector3& v, const vector3::DataType a) 
	{return vector3(v.x * a, v.y * a, v.z*a);}
inline bool operator==(const vector3& u, const vector3& v) 
	{return (u.x == v.x && u.y == v.y && u.z==v.z);}
inline bool operator!=(const vector3& u, const vector3& v) 
	{return (u.x != v.x || u.y != v.y || u.z!=v.z);}
inline vector3 abs(const vector3& v) 
	{return vector3(fabs(v.x), fabs(v.y), fabs(v.z));}
inline vector3 cross( const vector3 &v1, const vector3 &v2) {
       return vector3( (v1.y*v2.z-v2.y*v1.z),
                      (v2.x*v1.z-v1.x*v2.z),
                      (v1.x*v2.y-v2.x*v1.y) );
}

/* Secondary operations */
inline vector3 operator-(const vector3& u, const vector3& v) {return u + -v;}
inline vector3 operator*(const vector3::DataType a, const vector3& v) {return v * a;}
inline vector3 operator/(const vector3& v, const vector3::DataType a) {return v * (1.0/a);}

#endif 

