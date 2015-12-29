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


#ifndef MESOFLUID_H
#define MESOFLUID_H
#include <vector>
#include "Vector.h"
#include "pairlist.h"
#include "constants.h"

using namespace std;

struct Parameters {
	int Nparticles;
	double L;
	int Niter;
	double rcut;
	double pairdist;
	double tau;
	char* infname;
	double a;
	double gamma;
	double sigma;
	double mass;
	};


class Fields { 
public:
	Vector P;
	Fields() :  P(0.0,0.0,0.0) {}
};

class Pair {
public:
	int k,l;
	Pair(int _k,int _l) :k(_k),l(_l) {}
};

class Mesofluid  {
public:
	Mesofluid(Vector* pp,Fields* fp,int _np,const Distance& _domain, const Parameters& par);
	Mesofluid(const Parameters& par, const Distance& _domain);
	void integrate_euler();
	void integrate_DPD_VV();
	void integrate_trotter();
	void state(double& T, Vector& P) const;
	int read_pos(const char* filename);
	void write_pos(const char* filename, bool wrap=1);
	Parameters Par;
private:
	Distance domain;
	Pairlist pl;
	vector<Pair> pairs;
	vector<Fields> fp;
	vector<Vector> r;
	int np;
	double compute_forces(double dt, vector<Fields>& dfp, int NC=1);
	double compute_forces_trotter(double dt, int order=0);

};

#endif
