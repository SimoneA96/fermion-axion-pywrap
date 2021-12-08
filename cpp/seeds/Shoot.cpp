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
 *   people that develope it.
 *
 *   Here there are the funcions used by Shoot.cpp and Solve.cpp in order
 *   to build the 2-D existence plot of a fermion-axion star. 
 *
 */


/*
 *
 *   Revision 1.1  2021/12/02  20:30:00  Davide Guerra
 *       *** first patch ***
 *
 *
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include "shooting.h"

using namespace std;


/********** MAIN ***********/

int main(int argc, char *argv[]){
	int input_num;
	long double order=-19, phi0c, rho0c, omega_mid;
	
 
	if (argc<2){
	  cout<<"You should run './Shoot data' "<<'\n';
	  cout<<"where data is the folder where you want to put your output"<<endl;
	  return 0;
	}

// defining input and output files/directories

  ifstream  input;
 	ofstream  output;
	string folder;
	folder = argv[1];
	char* folder2;
	folder2 = argv[1];
	
  mkdir( folder2, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0;
   	input.open("Input.txt", ios::in);
   	output.open(folder + "/Output.txt", ios::out);
   	input >> order >> num >> input_num;

	var = new long double*[5];
	for (int i = 0;i<5;i++)
		var[i] = new long double[num];

  output << input_num <<" "<< num <<endl;
  output << setprecision(2-order);
 	cout << setprecision(2-order); 


  if (order<0) {int j=(-1)*order;
    order=1;
		for (int i=0; i<j; i++) {
		  order=order*0.1;
		}
    cout << order << "\n";
 	}

  // reading from the input file 	
  for (int j=0; j<input_num; j++) {
    input >> phi0c >> rho0c >> R;
    Rmax=R;
		cout << j <<" "<<phi0c<<" "<<rho0c<< endl;
		
    // call the shooting method function
		mass = shooting_omega(phi0c, rho0c, order, &omega_mid, &grr);
   
    // printing the output
    output <<  phi0c << " " << KK*pow(rho0c,polind) << " " << omega_mid <<" "<<R<< endl;
    cout<<"grr= "<<grr<<" mass= "<<mass<<endl;
	}

	for (int i=0;i<5;i++)
	  delete [] var[i];
 
  delete [] var;
  input.close();
	output.close();
	return 0;
}
