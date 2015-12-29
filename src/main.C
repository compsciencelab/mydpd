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
#include <stdio.h>
#include "mesofluid.h"
#include <fstream.h>
#include "random.h"
#include "timer.h"

struct Controls {
		int transient;
		int update_step;
		char* fname;

};

	
Timer TT[10];


int main(int argc,char* argv[]) 
{	
	if (argc!=5) {
		cout<<endl<<"USAGE: "<<argv[0]<<" dt N input_coo dat_file_name"<<endl;
		exit(1);
	}
TT[0].tic();
	print_time(); 
	Parameters Par;
	Controls C;

	Par.Nparticles=4000;
	Par.L= 10;
	Par.Niter=atoi(argv[2]);
	Par.rcut=1.0;
	Par.pairdist=1.0;
	Par.tau=atof(argv[1]);
	Par.infname =  argv[3];
	Par.a=25;
	Par.gamma=4.5;
	Par.sigma=3.0; 
	Par.mass=1.0;
	C.transient=0;
	C.update_step=1;
	C.fname = argv[4];
 
        Distance domain(Par.L,Par.L,Par.L,Vector(0,0,0));
	
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	init_by_array(init, length);

	Mesofluid model(Par, domain );
	
	ofstream os(C.fname);os.close();//reset file

	double t=0.0;
	
	for (int n=0;n<Par.Niter;n++)
	{
		if((n%C.update_step==0)&&(n>=C.transient))
		{
		double T;Vector tot_P;
	        model.state(T,tot_P);		
		cout<<"n="<<n<<" t="<<t<<" temp="<<T<<" total_momentum = "<<tot_P.length()<<endl;
		ofstream os(C.fname,ios::app);
		os<<t<<' '<<T<<'\n';
		os.close();
		}
TT[2].tic();
		model.integrate_trotter();
TT[2].toc();
		t+=model.Par.tau;
	}
TT[0].toc();
	int tt;
	for(tt=0;tt<3;tt++) cout<<" "<<TT[tt].time();
	cout<<'\n';
	model.write_pos("pos.coo");
	print_time();
}

