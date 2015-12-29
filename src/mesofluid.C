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



#include <fstream.h>
#include "mesofluid.h"
#include "random.h"

using namespace std;

	
inline Vector Maxwell_Boltzmann(double m, const Vector& mean_v, double kT) {
        return sqrt(kT/m)*Vector(normal(), normal(), normal()) + mean_v;
}


void compute_pairs(Pairlist& pl, Distance& domain, vector<Vector>& r, double RCUT, vector<Pair>& pairs) {
	Vector pos[r.size()];
	for (int i=0;i<r.size();i++) pos[i]=r[i];
	pl.compute(domain,pos,r.size());
	pairs.clear();
	for (int k = 0; k < pl.boxlist.size(); k++)//iterate over boxes 
	for (int nbr=0; nbr<pl.nbrlist[k].size(); nbr++) {//iterate over nbr boxes
		int k1 =  pl.nbrlist[k][nbr];
		for (int l = 0; l < pl.boxlist[ k ].size(); l++) {//iterate over atoms
			int i =  pl.boxlist[k][l];
			int startl1=0;
			if (k1 == k) startl1 = l + 1;
			for (int l1=startl1;l1<pl.boxlist[ k1 ].size();l1++) {//iterate over atoms
				int j = pl.boxlist[k1][l1];
				Vector r_kl = domain.delta(r[i],r[j]);
				double R_kl=r_kl.length();
                		if (R_kl<RCUT && i!=j) 
                		 pairs.push_back(Pair(i,j));
			}
		}
	}
//	cout<<" pairlist = "<<pairs.size()<<endl;
}
/*
void compute_pairs(Pairlist& pl, Distance& domain, vector<Vector>& r, double RCUT, vector<Pair>& pairs) {
	pairs.clear();
	for (int i = 0; i < r.size(); i++)
	for (int j = 0; j < r.size(); j++) {
				Vector r_kl = domain.delta(r[i],r[j]);
				double R_kl=r_kl.length();
                		if (R_kl<RCUT && i<j) pairs.push_back(Pair(i,j));
	}
	cout<<" pairlist = "<<pairs.size()<<endl;
}
*/

Mesofluid::Mesofluid(const Parameters& par, const Distance& _domain)  
: Par(par),  domain(_domain) , pl(Pairlist(par.pairdist,1))
{
 	np = read_pos(par.infname);
	fp.resize(np,Fields());
	double  T;
	if (Par.gamma>0) T=0.5*Par.sigma*Par.sigma/Par.gamma;	
	else T=0;
	for (int i=0;i<np;i++) {
		fp[i].P=Par.mass*Maxwell_Boltzmann(Par.mass, Vector(), T);
	}
	Vector P;
	state(T,P);
	for (int i=0;i<np;i++) fp[i].P-=P;//remove average momentum

	compute_pairs(pl,domain,r,par.rcut,pairs);

	cout<<" np = "<<np<<" a = "<<Par.a<<" sigma ="<<Par.sigma<<" gamma = "<<Par.gamma<<" dt ="
		<<Par.tau<<" mass ="<<Par.mass<<" rcut = "<<Par.rcut<<endl;
}

int    Mesofluid::read_pos(const char* filename) {
	ifstream is(filename);
	char buf[100];
	while ( is.getline(buf,100) ) {
		float x,y,z;
		sscanf(buf,"%f %f %f",&x,&y,&z);
		r.push_back(Vector(x,y,z));
	}
	is.close();
	return r.size();
}

void Mesofluid::write_pos(const char* filename, bool wrap) {
	ofstream os(filename);
	for (int i=0;i<np;i++) {
		Vector pos; 
		if (wrap) pos = domain.image(r[i]);
		else pos = r[i];
		os<<pos.x<<" "<<pos.y<<" "<<pos.z<<"\n";
	}
	os.close();
}
 
void   Mesofluid::state(double& T, Vector& P) const
{
	double tot=0.0;
	Vector tmp;
	for(int i=0;i<np;i++)
	{ 
		tot+=0.5*fp[i].P*fp[i].P/Par.mass;
		tmp+=fp[i].P;
	}
	T = 2.0/3.0*tot/double(np-1);
	P=tmp/double(np);
}

		
inline double omega_R(double r, double rcut) {
	return 1.0 -  r/rcut;
}

inline double cforce(double r, double rcut) {
	return (1.0 -  r/rcut);
}

const double SQRT3 = sqrt(3.0);

inline double weakrnd() {
	double tmp = uniform(0.0,1.0);
	if (tmp>=1.0/3.0) return 0.0;
	else if (tmp<1.0/6.0) return SQRT3;
	else return -SQRT3;
}

double    Mesofluid::compute_forces(double dt, vector<Fields>& F, int NC)
{
	double Ec=0.0;
	Vector Fr_kl,Fd_kl,Fc_kl;
	for(int i=0;i<np;i++) F[i]=Fields();

	double sqrtdt =  sqrt(dt) ;
	for (int i=0;i<pairs.size();i++) {
		int k=pairs[i].k;
		int l=pairs[i].l;
		Vector r_kl = domain.delta(r[k],r[l]);
		double R_kl=r_kl.length();
		if (R_kl==0.0) cout<<"R_KL=0"<<endl;
                if (R_kl<Par.rcut) {
			//k
			Fields fk=fp[k];
			Vector U_k=fk.P/Par.mass;
			//l
			Fields fl=fp[l];
			Vector U_l=fl.P/Par.mass;
			//geometry
			Vector e_kl(r_kl/R_kl);
			Vector U_kl(U_k-U_l);
			double w_R=omega_R(R_kl,Par.rcut);
			double w_D=w_R*w_R;
			

			Fc_kl = Par.a*cforce(R_kl,Par.rcut)*e_kl;
			Ec +=   0.5*Par.a*cforce(R_kl,Par.rcut)*cforce(R_kl,Par.rcut);
			Fr_kl = double(NC)*Par.sigma* w_R*normal()*e_kl/sqrtdt;
			Fd_kl = - double(NC)* Par.gamma*w_D*(U_kl*e_kl)*e_kl;
		
			F[k].P +=  Fc_kl  + Fd_kl  + Fr_kl ;
			F[l].P += -Fc_kl  - Fd_kl  - Fr_kl ;
		}
	}
return Ec;
}

inline double rho_p(int p) {
	double r=0.0;
	for(int i=1;i<=p;i++) r+=1.0/double(i*i);
	return 1.0/12.0 - 1.0/(2.0*Pi*Pi)*r;
}

double    Mesofluid::compute_forces_trotter(double dt, int order)
{
//	static double rhop =  rho_p(20);
	double Ec = 0.0;
	Vector Fr_kl,Fd_kl,Fc_kl;

	double sqrtdt =  sqrt(dt) ;

	int start,end,step;
	if (order==0) {start=0;end=pairs.size()-1;}
	else if (order==1) {start=-pairs.size()+1;end=0;}
	else abort();

	for (int ii=start;ii<=end;ii++) {
		int i = abs(ii);
		int k=pairs[i].k;
		int l=pairs[i].l;
//		cout<<" k = "<<k<<" l = "<<l<<endl;
		Vector r_kl = domain.delta(r[k],r[l]);
		double R_kl=r_kl.length();
		if (R_kl==0.0) cout<<"R_KL=0"<<endl;
                if (R_kl<Par.rcut) {
			//k
			Fields fk=fp[k];
			Vector P_k=fk.P;
			//l
			Fields fl=fp[l];
			Vector P_l=fl.P;
			//geometry
			Vector e_kl(r_kl/R_kl);
			Vector P_kl(P_k-P_l);
			double w_R=omega_R(R_kl,Par.rcut);
			double w_D=w_R*w_R;

			if (Par.gamma>0 && Par.sigma>0) {
			double gamma1 = 2.0*Par.gamma*w_D;
			Fr_kl = Par.sigma*w_R*e_kl*sqrt( (1.0-exp(-2.0*gamma1*dt))/(2.0*gamma1) )*normal();
			Fd_kl = 0.5*( P_kl*e_kl - Par.a*cforce(R_kl,Par.rcut)/(Par.gamma*w_D) )
				*e_kl*( exp(-2.0*Par.gamma*w_D*dt)-1.0);
			} else {
				Fr_kl*=0.0 ;
				Fd_kl = Par.a*cforce(R_kl,Par.rcut)*e_kl*dt;
			}
			Ec += 0.5*Par.a*cforce(R_kl,Par.rcut)*cforce(R_kl,Par.rcut);

			fp[k].P +=     Fd_kl  + Fr_kl ;
			fp[l].P +=    - Fd_kl  - Fr_kl ;
		}
	}
return Ec;
}

void Mesofluid::integrate_euler()
{
	vector<Fields> F(np,Fields());
	double dt = Par.tau;
	compute_pairs(pl,domain,r,Par.rcut,pairs);
	compute_forces(dt,F);

	for(int k=0;k<np;k++) r[k]+= fp[k].P/Par.mass*dt;
	for(int k=0;k<np;k++) fp[k].P += F[k].P*dt; 
	
}

void Mesofluid::integrate_DPD_VV()
{
static	vector<Fields> F(np,Fields());
	vector<Vector> tmp_P(np,Vector());
	vector<Vector> tmp_FP(np,Vector());

	double dt = Par.tau;
	double m = Par.mass; 
	for(int k=0;k<np;k++) r[k] += fp[k].P/m*dt + 0.5*F[k].P/m*dt*dt;
	for(int k=0;k<np;k++) {tmp_P[k] = fp[k].P;tmp_FP[k]=F[k].P;}
	for(int k=0;k<np;k++) fp[k].P += 0.5*F[k].P*dt; 
	compute_pairs(pl,domain,r,Par.rcut,pairs);
	double Ec = compute_forces(dt,F);        	
       	for(int k=0;k<np;k++) fp[k].P = tmp_P[k] + 0.5*(F[k].P + tmp_FP[k])*dt; 
        //cout<< "Conservative energy [DPD-VV] "<<Ec<<endl;

}

void Mesofluid::integrate_trotter()
{
	double dt2 = 0.5*Par.tau;
	double dt =  Par.tau;
	double m = Par.mass; 

	double Ec = compute_forces_trotter(dt2,0);//2
	for(int k=0;k<np;k++) r[k] += fp[k].P/m*dt;//1
	compute_pairs(pl,domain,r,Par.rcut,pairs);

	Ec = compute_forces_trotter(dt2,1);//2
	//cout<< "Conservative energy 2 [DPD-Trotter] "<<Ec<<endl;
		
}


