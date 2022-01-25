/*
 *  This code contains all the functions used in Shoot.cpp and Solve.cpp
 *  including the shooting method to evaluate the field frequency		
 *		 
 */
 
 
/*
 *   Copyright (c) 2021-2022 Simone Albanesi, Fabrizio di Giovanni,
 *   Davide Guerra, Miquel Miravet
 *
 *   This file is part of FAS (Fermion-Axion Star program).
 *
 *   FAS is free software (we accept a coffee or a beer if you want); 
 *   you can redistribute it and/or modify it under the terms of the 
 *   people that develope it.
 *
 *   Here there are the funcions used by Shoot.cpp and Solve.cpp in
 *   to build the 2-D existence plot of a fermion-axion star. 
 *
 */


/*
 *
 *   Revision 1.0  2021/12/02  20:30:00  Davide Guerra
 *       *** first version ***
 *
 *
*/


#include "shooting.h"

/***********************/
/*  GLOBAL PARAMETERS  */
/***********************/

long double y[5], y_axi[5], diff[5], F[5];              // variables for the ODE solver and the Shooting solver
long double** var;
long double k1[5], k2[5], k3[5], k4[5];


//const long double |KK=VALUE|, |m=VALUE|, 
//                  |fa_axion=VALUE|, |polind=VALUE|;     // constants related with the politropic eos for the
                                                        // fermion part (KK and polind) and the mass (m) of
                                                        // the scalar field and the decay constant (fa) for
                                                        // the axion potential
                                                                               

long double eps, grr, mass;                             // eps=pow(KK,1/polind)*pow(P,(polind-1)/polind)/(polind-1)
                                                        // grr is the r,r component of the metric tensor gmunu
                                                        // mass is the final mass of the fermion-axion model
                                                        

int num;                        
long double R, Rmax, omega2, Cutphi, CutP, delta_r, r;  // variable for the ODE and the RK4 solver


/**************************/
/*  FUNCTIONS DEFINITION  */
/**************************/
                      
// function that evaluate the sign of a variable (1 if positive, 0 if null, -1 if negative)                               
int sign (long double z) {
  if (z>0) return 1;
  if (z==0) return  0;
  if (z<0) return -1;
}

// function that return the minimun between two input variables
long double min(long double z_1, long double z_2){
	if (z_1<z_2) return z_1;
	else return z_2;
}

// definition of the axion potential according from "D.Guerra, C.Macedo, P.Pani, Axion boson stars, 2019"
long double V(long double phi) {
  return 2.0 * pow(m,2) * pow(fa_axion,2) / 0.22 * (1.0 - sqrt(1.0 - 4 * 0.22 * pow(sin( phi/(2.0*fa_axion)),2)));
}
long double dV(long double phi) {
  if (phi!=0.0){
    return 2.0 * pow(m,2)  * fa_axion * sin(phi/(2.0*fa_axion)) * cos(phi/(2.0*fa_axion)) /
      (phi * sqrt(1.0 - 4.0 * 0.22 * pow(sin(phi/(2.0*fa_axion)),2) ));
  } else {
    return 0;
  }
}


// function used for the Runge Kutta calcolus 
int func(long double r,  long double omega2, long double h) {

  long double a=y_axi[0];
  long double alp=y_axi[1];
  long double phi=y_axi[2];
  long double dphi=y_axi[3];
  long double P=y_axi[4];
  long double rho_tmp;
    
  //OSS: eps=K*rho --> eps=10*sqrt(P) eps = pow(KK,1/polind)*pow(P,(polind-1)/polind)/(polind-1)
  eps = pow(KK,1/polind)*pow(P,(polind-1)/polind)/(polind-1);
  rho_tmp = pow(P/KK,1/polind);

  F[0] = 0.5*a*((1-a*a)/r + 8*Pi*r*( pow(dphi,2) + omega2 * pow(a*phi/alp, 2) + a*a * V(phi)+a*a*rho_tmp*(1+eps)));
  F[1] = 0.5*alp*((a*a-1)/r + 8*Pi*r*( pow(dphi,2) + omega2 * pow(a*phi/alp, 2) - a*a * V(phi)+a*a*P));
  F[2] = dphi;
  F[3] = -(1+a*a-8*Pi*pow(r*a,2)*(V(phi)+0.5*(rho_tmp*(1+eps)-P)))*dphi/r+(dV(phi)-omega2/pow(alp,2))* a*a * phi;
  F[4] =-(rho_tmp*(1+eps)+P)*0.5*((a*a-1)/r + 8*Pi*r*( pow(dphi,2) + omega2 * pow(a*phi/alp, 2) - a*a * V(phi)+a*a*P));

  return 0;
}


// Runge-Kutta to the 4th order algorithm
int Runge_Kutta_4(long double r,  long double omega2, long double h)  {
    
  long double half_h=0.5*h;

  for (int i=0; i<=4; i++) y_axi[i]=y[i];
    
  func(r, omega2, h);
  for (int i=0; i<=4; i++) k1[i]=h*F[i];
  for (int i=0; i<=4; i++) y_axi[i]=y[i]+0.5*k1[i];

  //control for P not to become negative
  if (y_axi[4]<=0) {
	  y_axi[4]=0;
		k1[4]=0;
  }
    
  func(r+half_h, omega2, h);
  for (int i=0; i<=4; i++) k2[i]=h*F[i];
  for (int i=0; i<=4; i++) y_axi[i]=y[i]+0.5*k2[i];
    
  //control for P not to become negative
  if (y_axi[4]<=0) {
	  y_axi[4]=0;
		k2[4]=0;
  }
    
  func(r+half_h, omega2, h);
  for (int i=0; i<=4; i++) k3[i]=h*F[i];
  for (int i=0; i<=4; i++) y_axi[i]=y[i]+k3[i];
   
  //control for P not to become negative
  if (y_axi[4]<=0) {
	  y_axi[4]=0;
		k3[4]=0;
  }
    
  func(r+h, omega2, h);
  for (int i=0; i<=4; i++) k4[i]=h*F[i];
  for (int i=0; i<=4; i++) diff[i]=(k1[i]+k4[i])/6 + (k2[i]+k3[i])/3;
    
  return 0;
}


// solver for ODE equations
int ode_int(long double phi0c, long double rho0c, long double omega, long double* a){
    
  int n=0;

  y[0]=1;
  y[1]=1;
  y[2]=phi0c;
  y[3]=0;
  y[4]=KK*pow(rho0c,polind);
  delta_r=R/(num);
  omega2=pow(omega,2);
  r=delta_r;
  
    
  while (r<R) {
    for (int i=0; i<=4; i++) {
      var[i][n]=y[i];
    }
    n++;
    Runge_Kutta_4(r, omega2, delta_r);
	
    for (int i=0; i<=4; i++) {
      y[i]+=diff[i];
    }
        
    //control for P not to become negative
	if(y[4]<=0) y[4]=0;

	//Setting to 0 the pressure when it starts to do steps (for low values of P)
    if (y[4] == var[4][n-1] || var[4][n-1]==0) y[4]=0;


    if (y[2] > Cutphi || y[2] < 0 || y[4]> CutP || diff[2]>0) {
		  Rmax = r;
      break;
    }
        
    r+=delta_r;
  }
    
  *a = y[0];

  return sign(y[2]);

}



// Shooting method that returns the mass of the fermion-axion star and
// evalutate the value of the field frequency saving in the pointer 
long double shooting_omega(long double phi0c, long double rho0c, long double ord, long double* mid_omega, long double* a){

  //Cuts for phi and P: if the solution goes over these values it breaks. To avoid having nans. 
  Cutphi = 1.1*phi0c;
  CutP = 1.1* KK*pow(rho0c,polind);
     		
  long double omega_1=1.0*m,omega_2=2.0*m;
  long double f_1=ode_int(phi0c, rho0c, omega_1, a);
  long double f_2=ode_int(phi0c, rho0c, omega_2, a);
  long double f_mid;
  
    
  while (f_1*f_2>0) {
    omega_2=omega_2*1.1;
    omega_1=omega_1/1.1;
    f_1=ode_int(phi0c, rho0c, omega_1, a);
    f_2=ode_int(phi0c, rho0c, omega_2, a);
			
    //control for nans
	  if (f_1*f_2 != 1 && f_1*f_2 != -1 && f_1*f_2!=0){
			cout<<"nan detected;"<<'\n';
			cout<<"omega_1= "<<omega_1 <<" omega_2= "<<omega_2<<endl; 
			cout<<"f_1= "<<f_1 <<" f_2= "<<f_2<<endl; 
		}
		// end control
     	
  }
    
  while((omega_2-omega_1)/omega_1>2*ord) {
    *mid_omega=(omega_1+omega_2)/2;
    cout << *mid_omega << "  ";
    f_mid=ode_int(phi0c, rho0c, *mid_omega, a);
			
		//control for nans
		if (f_mid != 1 && f_mid != -1 && f_1*f_2!=0){
			cout << "nan detected;" <<'\n';
			cout << "omega_mid= " << *mid_omega << endl; 
			cout << "f_mid= " << f_mid << endl; 
		}
		// end control
      
    if (f_mid>0){
      if(f_1<0){
			  omega_2=*mid_omega;
				f_2=f_mid;
      }
		  else {
		    omega_1=*mid_omega;
			  f_1=f_mid;
		  }
    }
		else{
		  if(f_1<0){
			  omega_1=*mid_omega;
				f_1=f_mid;
      }
			else {
			  omega_2=*mid_omega;
				f_2=f_mid;
			}
    }

  }
		
  cout<<endl;
  
  return ((Rmax/2.) * (1-1/pow(*a,2)));
  
  
}
