/*
 *  Header with all the functions used in Shoot.cpp and Solve.cpp
 *  including the shooting method to evaluate the field frequency		
 *		 
 */
 
 
/*
 *   Copyright (c) 2021-2022 Simone Albanesi, Fabrizio di Giovanni,
 *   Davide Guerra
 *
 *   This file is part of FAS (Fermion-Axion Star program).
 *
 *   FAS is free software (we accept a coffee or a beer if you want); 
 *   you can redistribute it and/or modify it under the terms of the 
 *   3 people that develope it.
 *
 *   Here there are the funcions used by Shoot.cpp and Solve.cpp in order
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

#include <iostream>
#include <cmath> 
#include <fstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>

#define DEBUG 0

using namespace std;

#ifndef HEADER_H
#define HEADER_H

/***********************/
/*  GLOBAL PARAMETERS  */
/***********************/

#define Pi 3.14159265358979323846            // value of pi

#define KK 100                           // constants related with the politropic eos for the
#define m 1.0                            // fermion part (KK and polind) and the mass (m) of
#define fa_axion 0.02                     // the scalar field and the decay constant (fa) for
#define polind 2.                       // the axion potential


extern long double y[5], y_axi[5], diff[5], F[5];               // variables for the ODE solver and the Shooting solver
extern long double** var;
extern long double k1[5], k2[5], k3[5], k4[5];                                                                             

extern long double eps, grr, mass;                             // eps=pow(KK,1/polind)*pow(P,(polind-1)/polind)/(polind-1)
                                                               // grr is the r,r component of the metric tensor gmunu
                                                               // mass is the final mass of the fermion-axion model

extern int num;                        
extern long double R, Rmax, Cutphi, CutP, omega2, delta_r, r;  // variable for the ODE and the RK4 solver

/**************************/
/*  FUNCTIONS DEFINITION  */
/**************************/
                      
// function that evaluate the sign of a variable (1 if positive, 0 if null, -1 if negative)                               
int sign (long double z);


// function that return the minimun between two input variables
long double min(long double z_1, long double z_2);


// definition of the axion potential according from "D.Guerra, C.Macedo, P.Pani, Axion boson stars, 2019"
long double V(long double phi);
long double dV(long double phi);


// function used for the Runge Kutta calcolus 
int func(long double r,  long double omega2, long double h);


// Runge-Kutta to the 4th order algorithm
int Runge_Kutta_4(long double r,  long double omega2, long double h);


// solver for ODE equations
int ode_int(long double phi0c, long double rho0c, long double omega, long double* a);


// Shooting method that returns the mass of the fermion-axion star and
// evalutate the value of the field frequency saving in the pointer 
long double shooting_omega(long double phi0c, long double rho0c, long double ord, long double* mid_omega, long double* a);

#endif