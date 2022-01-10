/*
 *  The aim of this code is to compute the equilibrium configurations
 *  for a Fermion-Axion star system. In order to use this code please
 *  read the README file because it is all automated with the code
 *  wrapper.py.
 *
 *  OUTPUT: txt file with 4 columns
 *      1. phi0: central value of the scalar axion field
 *      2. P=k*(rho0)^gamma : value of the central pressure of the 
 *                            fermion part of the star (polytrope EOS)
 *      3. omega_shooting: frequency result from the shooting method
 *      4. R: radius of the fermion-axion star
 *		 
 */
 
 
/*
 *   Copyright (c) 2021-2022 Saeed Fakhry (the main author of this code),
 *   Simone Albanesi, Fabrizio di Giovanni, Davide Guerra, Miquel Miravet
 *
 *   This file is part of FAS (Fermion-Axion Star program).
 *
 *   FAS is free software (we accept a coffee or a beer if you want); 
 *   you can redistribute it and/or modify it under the terms of the 
 *   people that develope it, otherwise ciao ciao.
 *
 *   Here there are the funcions used by Shoot.cpp and Solve.cpp in order
 *   to build the 2-D existence plot of a fermion-axion star. 
 *
 */


/*
 *
 *   Revision 1.1  2021/12/02  09:27:00  Davide Guerra
 *       *** first patch ***
 *
 *
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <string>
#include "shooting.h"

using namespace std;


/****************/
/*  PARAMETERS  */
/****************/

long double* riso, *riso2, *beta;
long double dpsidr, ADM_mass;
long double phi0c, P0c, omega, l_1, l_2;
long double rtmp;
long double numbos, numfer, numbos2,numfer2, r99bos, r99fer, r99tot, rbos_iso, rfer_iso, r_isotropic_max, factor, rmax;
int nbos,nfer, nmax;
long double massbos, massfer, massT, massT2;
long double* conformal, *conformal2;
bool part, isotropic;
long double* bospotential;

long double epsilon = 1e-8;
int flag = 0;


/**********/
/*  MAIN  */
/**********/

int main(int argc, char *argv[]){

	
	if (argc<2){
	  cout<<"You should run './Solve data' "<<'\n';
	  cout<<"where data is the folder where you want to put your output"<<endl;
	  return 0;
	}

	int n, ntmp, input_num=1;
	
  // setting the input/output files and directories
  ifstream  input;
	string folder;
	folder = argv[1];
	string tmp;
	string fa_str;
	string phi_str;
	string rho_str;
	input.open(folder + "/Output.txt", ios::in );
	input >> input_num >> num;

	//dynamical allocation
	var = new long double*[5];
	for (int i = 0;i<5;i++)
		var[i] = new long double[num];

	riso = new long double[num];
	riso2 = new long double[num];
	beta = new long double[num];
	conformal = new long double[num];
	conformal2 = new long double[num];
	bospotential = new long double[num];

	ofstream rfile;
	ofstream solfile;
	ofstream solfile1;
	ofstream solfile2;
	ofstream solfile3;
	ofstream solfile4;
	ofstream solfile5;
	ofstream alpsolfile;
	ofstream asolfile;
  ofstream evolutionfile;
	ofstream massfile;
	ofstream conffile;
	ofstream infofile;
	ofstream bospotential_file;
	ofstream existence_plot;
  ofstream blackline;

	fa_str = to_string(fa_axion);

  // reading the input file
	for(int j=1; j<=input_num; j++) {
		y[0]=1;
		y[1]=1;
		input >> phi0c >> P0c >> omega >> R;
		cout<< phi0c<<" "<<P0c<< " "<<omega<<endl;
		y[2]=phi0c;
		y[3]=0;
		y[4]=P0c;

		n=0;
		omega2=pow(omega,2);
		delta_r=R/(num);
		cout << setprecision(9);
		r=delta_r;
		part = 0;
		numbos = 0.0;
		numfer = 0.0;
   
    while (r<=R  +delta_r) {
			for (int i=0; i<=4; i++) {
				var[i][n]=y[i];
			}
	    n++;
	    Runge_Kutta_4(r, omega2, delta_r);
	    for (int i=0; i<=4; i++) {
	      y[i]+=diff[i];
    	}
    
      //control for P not to become negative
			if(y[4]<=0)	y[4]=0;

			if ((y[2] < 0) || diff[2] >0) {
				if(flag==0){
					cout<<r<<endl;
					flag=1;
					rmax = r;
					nmax = n;
				}

        // Define a value for phi, because of the derivative of axion potential has phi at the denominator! 		
        y[2]=epsilon * phi0c;
				y[3]=0;    		
			}
			
      if (y[4] == var[4][n-1] || var[4][n-1]==0) 	y[4]=0;

			r+=delta_r;
		}

    //number of bosons and fermions
    tmp = to_string(j);
		r=delta_r;
		n=0;
		numbos=0;
		numfer=0;
		
    while (r<=R) {
			numbos = numbos + 4.0*M_PI*var[0][n]*omega*pow(var[2][n]*r,2.0)/var[1][n]*delta_r;
			numfer = numfer + 4.0*M_PI*var[0][n]*pow(var[4][n]/KK,1/polind)*pow(r,2.0)*delta_r;
	    n++;
			r+=delta_r;
		}

		r=delta_r;
		n=0;
		numbos2=0;
		numfer2=0;

    //R99 or R95 for boson and fermions and total

		while (r<R) {
			numbos2 = numbos2 + 4.0*M_PI*var[0][n]*omega*pow(var[2][n]*r,2.0)/var[1][n]*delta_r;
			numfer2 = numfer2 + 4.0*M_PI*var[0][n]*pow(var[4][n]/KK,1/polind)*pow(r,2.0)*delta_r;
	    		n++;
			if(numbos2/numbos<0.99){
				r99bos=r;
				nbos = n;
				}
			if(numfer2/numfer<0.99){
				r99fer=r;
				nfer = n;
				}

			r+=delta_r;
		}

    //Change of R99 to isotropic coordinates
		r_isotropic_max = pow( (1+sqrt(var[0][num-1]))/2,2)*R/var[0][num-1];


    //Isotropic radius
	  isotropic=0;
	  if(isotropic){	
	    r=delta_r;
	    n=0;
	    while(r<R){
		    rtmp=r;
  		  factor=0;
  		  ntmp=n;
  		  while(rtmp<R){

	  		  factor = factor + var[0][ntmp]*delta_r/rtmp;
	  		  ntmp++;
	  		  rtmp+=delta_r;		
		
	  	  }
	  	  riso[n]=r_isotropic_max*exp(-factor);
	  	  r+=100*delta_r;
	  	  n+=100;
	    }

    //using beta-factor improves the change to isotropic coordinates (see arXiv:gr-qc/9707045 Sec. IIIA)
    	r=delta_r;
    	n=0;
    	while(r<R){
    		rtmp=r;
    		factor=0;
    		ntmp=n;
    		while(rtmp<R){
  
    			factor = factor + (var[0][ntmp]-1)*delta_r/rtmp;
    			ntmp++;
    			rtmp+=delta_r;		
	  	
    		}
	    	beta[n] =r_isotropic_max/R * exp(-factor);
	    	riso2[n]=r*beta[n];
	    	r+=100*delta_r;
	    	n+=100;  

	    }

	    r=delta_r;
	    n=0;

	    while(r<R){
	
		    conformal[n] = sqrt(r/riso[n]);
		    conformal2[n] = sqrt(1.0/beta[n]);
	  	  r+=100*delta_r;
  	  	n+=100;
  
    	}

	    dpsidr = (conformal2[n-400]-8.0*conformal2[n-300]+ 8.0*conformal2[n-100]-conformal2[n])* 1.0 / (12.0*delta_r*100);
    	ADM_mass = -2.0*riso2[n-200]*riso2[n-200]*dpsidr;
    	cout<<"dpsidr= "<<dpsidr<<" ADM mass = "<<ADM_mass<<endl;
    }

    //Boson and fermion mass as spatial volume integrals
		r=delta_r;
		n=0;
		massbos=0;
		massfer=0;
		
		while (r<=R) {
			massfer = massfer + 4.0*M_PI*var[0][n]*pow(r,2)*( pow(var[4][n]/KK,1/polind) + var[4][n]/(polind-1) ) *delta_r;
			massbos = massbos + 4.0*M_PI*var[0][n]*pow(r,2.0)*(0.5*( omega2/var[1][n]/var[1][n]*var[2][n]*var[2][n] + 
				1/var[0][n]*var[3][n]*var[3][n] + V(var[2][n]) ))*delta_r;
			bospotential[n] =  4.0*M_PI*var[0][n]*pow(r,2.0)*(0.5*( omega2/var[1][n]/var[1][n]*var[2][n]*var[2][n] + 
				1/var[0][n]*var[3][n]*var[3][n] + V(var[2][n]) ))*delta_r;
  		
      n++;
			r+=delta_r;
		}	
	

		massT=((R)/2)* (1-1/pow(var[0][num-1],2));
		massT2=0;
		r=delta_r;
		n=0;

		while (r<R) {
			massT2 = ((r)/2)* (1-1/pow(var[0][n],2));

			if(massT2>=0.99*massT) {
				r99tot=r;
				break;
				}
 		  
      n++;
			r+=delta_r;
		}


		//Normalization of the lapse
		r=delta_r;
		n=0;
		omega = omega/(var[0][num-1]*var[1][num-1]);
		while (r<R) {
			var[1][n] = var[1][n]/(var[0][num-1]*var[1][num-1]);
			n++;
			r+=delta_r;
		}

    tmp = to_string(j);

		phi_str = to_string(phi0c);
		rho_str = to_string(pow(P0c/KK,1/polind));
		
    rfile.open(folder + "/solr"+tmp+".txt");
    solfile.open(folder + "/solphi"+tmp+".txt");
    solfile1.open(folder + "/solP"+tmp+".txt");
    solfile2.open(folder + "/solrho"+tmp+".txt");
    solfile3.open(folder + "/solrhoiso"+tmp+".txt");
    solfile4.open(folder + "/solphiiso"+tmp+".txt");
		solfile5.open(folder + "/solrhobosoniso"+tmp+".txt");
    alpsolfile.open(folder + "/solalp"+tmp+".txt");
    asolfile.open(folder + "/sola"+tmp+".txt");
    massfile.open(folder + "/solmass"+tmp+".txt");
    conffile.open(folder + "/conformal"+tmp+".txt");
    infofile.open(folder + "/info"+tmp+".txt");
    evolutionfile.open(folder + "/Mixed_"+phi_str +"_"+ rho_str +"_"+ fa_str+".dat");
    bospotential_file.open(folder + "/bospotential"+tmp+".txt");
		evolutionfile<<delta_r*10<<" "<<omega<<endl;
    blackline.open(folder + "/blackline"+tmp+".txt");

		//change number here, depending on how many models for each rho 
		if(j%input_num==1){
	           existence_plot.open(folder + "/existence_plot"+rho_str+".txt");
		  }

		n=0;
		r=delta_r;
		while (r<=R) {
      if (n%10 == 0) {
  			evolutionfile << r <<" "<<var[2][n]<<" "<<pow(var[1][n],2.)<<" "<<pow(var[0][n],2.)
					<<" "<<var[4][n]<<" "<<	pow(var[4][n]/KK,1/polind)<<"\n";	
      }


			if (n%10 == 0) {
				rfile << r  << "\n";

				solfile << r<<" "<<var[2][n]  << "\n";

				solfile1 << r<<" "<<var[4][n]  << "\n";

				solfile2 << r <<" "<< pow(var[4][n]/KK,1/polind)  << "\n";
	
				alpsolfile << r <<" "<< var[1][n] << "\n";

				asolfile << (r)<<" "<< var[0][n] <<"\n";

				massfile << r <<" "<< ((r)/2)* \
				(1-1/pow(var[0][n],2)) << "\n";	

				solfile3 <<riso[n]<<" "<< pow(var[4][n]/KK,1/polind)  << "\n";

				solfile4 <<riso[n]<<" "<< var[2][n]  << "\n";

				solfile5 <<riso[n]<<" "<<0.5*( omega2/var[1][n]/var[1][n]*var[2][n]*var[2][n] + 
					1/var[0][n]/var[0][n]*var[3][n]*var[3][n] + V(var[2][n])) <<"\n";

				conffile <<riso[n]<<" "<< conformal[n]  << " "<<riso2[n]<<" "<<conformal2[n]<< "\n";

				bospotential_file << r <<" "<< bospotential[n] << "\n";	
	

			
			}
	
			n++;
			r+=delta_r;
		}

		existence_plot<<var[2][0]<<" "<<omega<<" "<<((R)/2)* (1-1/pow(var[0][num-1],2))<<" "<< pow(var[4][0]/KK,1/polind)<<endl;  
    blackline << numbos << " " << numfer << " " << endl;
		rfile << endl;
		solfile << endl;
		solfile1 << endl;
		solfile2 << endl;
		solfile3 << endl;
		solfile4 << endl;
		solfile5 <<endl;
		alpsolfile << endl;
		massfile << endl;
		conffile << endl;
		asolfile << endl;
		evolutionfile<<endl;
		bospotential_file<<endl;


    cout<< "mass= "<<((R)/2)* (1-1/pow(var[0][num-1],2))<<" numbos: "<<numbos<< " numfer: "
      <<numfer<<" omega: "<<omega<<" r99bos: "<<r99bos<<" r99fer: "<<r99fer<<" r99tot: "<<r99tot<<endl;
		cout<<" Nb/Nf: "<<numbos/numfer<< " Eb/Ef: "<<massbos/massfer<<endl;
    cout<< "mass_bos= "<<massbos<<" mass_fer= "<<massfer<<endl;
		cout<<"compactness= "<<((R)/2)* (1-1/pow(var[0][num-1],2))/r99tot<<endl;
		cout<<"compactness_boson= "<<massbos/r99bos<<endl;
	  //cout<< "mass_corr= "<<((R)/2)* (1-1/pow(var[0][num-1],2))<<endl;
		cout<< "E_binding= "<< ((R)/2)* (1-1/pow(var[0][num-1],2)) - numbos - numfer <<endl;

  	if(isotropic)
		 cout<<"r99bos_iso: "<<riso[nbos]<<" r99fer_iso: "<<riso[nfer]<<endl;


    infofile<< "mass= "<<((R)/2)* (1-1/pow(var[0][num-1],2))<< "numbos: "<<numbos<< " numfer: "
      <<numfer<<" omega: "<<omega<<" r99bos: "<<r99bos<<" r99fer: "<<r99fer<<" r99tot: "<<r99tot<<endl;
		infofile<<" Nb/Nf: "<<numbos/numfer<< " Eb/Ef: "<<massbos/massfer<<endl;
    infofile<< "mass_bos= "<<massbos<<" mass_fer= "<<massfer<<endl;
		infofile<<"compactness= "<<((R)/2)* (1-1/pow(var[0][num-1],2))/r99tot<<endl;
		infofile<<"compactness_boson= "<<massbos/r99bos<<endl;
	  //infofile<< "mass_corr= "<<((rmax)/2)* (1-1/pow(var[0][nmax],2))<<endl;
    infofile<< "E_binding= "<< ((R)/2)* (1-1/pow(var[0][num-1],2)) - numbos - numfer <<endl;

	  if(isotropic)
      infofile<<"r99bos_iso: "<<riso[nbos]<<" r99fer_iso: "<<riso[nfer]<<endl;

    rfile.close();
	  solfile.close();
	  solfile1.close();
	  solfile2.close();
	  solfile3.close();
	  solfile4.close();
	  solfile5.close();
	  alpsolfile.close();
	  massfile.close();
	  conffile.close();
	  asolfile.close();
	  evolutionfile.close();
	  infofile.close();
	  bospotential_file.close();
    blackline.close();
		if(j%input_num==0){
      existence_plot.close();
    }
  }

	for(int i = 0; i < 5; i++) {
    delete [] var[i];
	}
  delete [] var;
  delete [] riso;
	delete [] riso2;
	delete [] beta;
	delete [] conformal;
	delete [] conformal2;
	delete [] bospotential;

	input.close();

	return 0;
}
