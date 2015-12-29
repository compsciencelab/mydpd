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



#include <iostream.h>
#include "pairlist.h"

using namespace std;

inline void find_minmax(const Vector* r, int n, Vector& min, Vector& max) {
	min=max=r[0];
        for(int i=1; i < n; i++)   {
        const Vector* tmp = r+i;
        if(tmp->x < min.x) min.x = tmp->x; else if(tmp->x > max.x) max.x = tmp->x;
        if(tmp->y < min.y) min.y = tmp->y; else if(tmp->y > max.y) max.y = tmp->y;
        if(tmp->z < min.z) min.z = tmp->z; else if(tmp->z > max.z) max.z = tmp->z;
        }
}


int Pairlist::findbox(Vector p, const Distance& dist) const {
	p = dist.delta(p) - lmin  ;
	int axb = (int)( p.x / pairdist.x) ;
	int ayb = (int)( p.y / pairdist.y) ;
	int azb = (int)( p.z / pairdist.z) ;
	return  azb * xytotb + ayb * xb + axb; 
}


void Pairlist::set_neighbors(const Distance& domain) {
static  int oxb=-1,oyb=-1,ozb=-1;

	if (xb!=oxb || yb!=oyb || zb!=ozb) {  //boxes changed
#ifdef DEBUG
	cout<<"PairList: setting neighbors: xb "<<xb<<" yb "<<yb<<" zb "<<zb<<" pairdist = "<<pairdist<<endl;
#endif

	//IMPORTANT NOTE: 
	//Resizing does not require to reallocate, but it means that
	//what is left has to be cleared with the loop

	nbrlist.resize(xb*yb*zb);	
	fullnbrlist.resize(xb*yb*zb);

	for (int i=0;i<xb*yb*zb;i++) {
		fullnbrlist[i].clear();
		nbrlist[i].clear();
	}

	xytotb = xb*yb;
	oxb=xb; oyb=yb; ozb=zb;
	
	for (int zi=0; zi<zb; zi++) 
        for (int yi=0; yi<yb; yi++) 
        for (int xi=0; xi<xb; xi++) 
        for (int zz=-1; zz<=1; zz++)   
        for (int yy=-1; yy<=1; yy++)  
        for (int xx=-1; xx<=1; xx++)  {
		bool cond =  (xx  > 0) || (xx == 0 && yy  == 1) || ( zz > -1 && xx == 0 && yy == 0);
		int axx = xx + xi; 
	       	int ayy = yy + yi; 
		int azz = zz + zi;
		bool flag = true;    
		if (domain.a()) axx = axx - int (xb * floor (float (axx) / float (xb)));
		else if (axx < 0 || axx >= xb) flag = false;
		if (domain.b()) ayy = ayy - int (yb * floor (float (ayy) / float (yb)));
		else if (ayy < 0 || ayy >= yb) flag = false;
		if (domain.c()) azz = azz - int (zb * floor (float (azz) / float (zb)));
		else if (azz < 0 || azz >= zb) flag = false;
		if (cond && flag) nbrlist[ zi * xytotb + yi * xb + xi ].push_back (  azz * xytotb + ayy * xb + axx );
       		if (flag) fullnbrlist[ zi * xytotb + yi * xb + xi ].push_back (  azz * xytotb + ayy * xb + axx );
	}
	}
}


void Pairlist::compute(const Distance& domain, const Vector* r, int np, int* aindex) {
//if (n%steps==0) {
#ifdef DEBUG1
	cout<<"Pairlist::compute"<<endl;
#endif

	double lx,ly,lz;
	lx = domain.a(); //set box size
	ly = domain.b();
	lz = domain.c();
	lmin = domain.origin() - Vector(0.5*lx, 0.5*ly, 0.5*lz);//in this case domain is centred 
	Vector min,max,side;
	if (!domain.a() || !domain.b() || !domain.c())  { //at least one direction is NOT periodic
		find_minmax(r,np,min,max);
		max+=margin;// 0.1A margin
		min-=margin;
		side = max - min ;
		if (!domain.a()) {lx = side.x; lmin.x = min.x;}
		if (!domain.b()) {ly = side.y; lmin.y = min.y;} 
		if (!domain.c()) {lz = side.z; lmin.z = min.z;} 
	}
#ifdef DEBUG1
        cout<<"Pairlist::compute set lx,ly,lz "<<lx<<" "<<ly<<" "<<lz<<" lmin "<<lmin<<" min "<<min<<" max "<<max<<endl;
#endif
	xb = (int)(lx / mindist);if (xb==0) xb=1; pairdist.x = lx/double(xb);
	yb = (int)(ly / mindist);if (yb==0) yb=1; pairdist.y = ly/double(yb);
	zb = (int)(lz / mindist);if (zb==0) zb=1; pairdist.z = lz/double(zb);
#ifdef DEBUG1
        cout<<"Pairlist:: neighbors xb,yb,zb "<<xb<<" "<<yb<<" "<<zb<<" pairdist "<<pairdist<<endl;
#endif
	set_neighbors(domain);

	boxlist.resize(xb*yb*zb);
//	pointlist.resize(xb*yb*zb);
	for (int i=0;i<xb*yb*zb;i++) { 
		boxlist[i].clear();
//		pointlist[i].clear();
	}
#ifdef DEBUG1
        cout<<"Pairlist::compute push_back atoms np  "<<np<<endl;
#endif

	for (int k=0;k<np;k++) {  
		int box = findbox(r[k],domain);
		if (box >= boxlist.size()) {
			cout<<" BOX ERROR "<<box<<" pos "<<r[k]<<" delta "<<domain.delta(r[k])<<endl;
		}
		int anum;
		if (aindex != 0) anum = aindex[k]; else anum = k;
		boxlist[box].push_back( anum );
//		pointlist[box].push_back( r[k] );
	}
#ifdef DEBUG1
        cout<<"Pairlist::End  "<<endl;
#endif

//}	
n++;
}

