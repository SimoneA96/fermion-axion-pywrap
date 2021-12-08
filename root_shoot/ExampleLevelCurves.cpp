/*
 *   Name convention for vertices and sides. The starting direction is 
 *   perpendicular to the sides and pointing outside the square (e.g.
 *   start_direction=3 means --->, start_direction=1 means <---).
 *
 *                2
 *        2---------------3
 *        |               | 
 *        |               |
 *      1 |               | 3
 *        |               |
 *        |               |
 *        1---------------4
 *                4
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
 *   DESCRIPTION NEEDED
 *
 */


/*
 *
 *   Revision 1.0  2021/12/02  21:30:00  Davide Guerra
 *       *** first version ***
 *
 *
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include "LevelCurves.h"
#include "shooting.h"


using namespace std;


int main(int argc, char *argv[]){
    double xmin, xmax, ymin, ymax;
    double x_init, y_init, ds,  start_direction, L, tol;
    int maxpoints, verbose;
    char fname[]="out_example.txt";

    x_init          = 0.001;
    y_init          = 0.0;
    start_direction = 2;
    ds              = 1e-2;  // used for initial step and for squares
    L               = 1e-4;  // distance(C,D) = 2*L (C and D points constructed for the bisection wo/ squares)
    xmin            = 0;
    xmax            = 0.008;
    ymin            = 0;
    ymax            = 0.42;
    tol             = 1e-4; // tolerance in bisection
    maxpoints       = 50;   // maximum number of points for the level curve 

    // Find the level curve and prin the ouput on "output.txt"
    FindLevelCurve(x_init, y_init, ds, start_direction, L, tol, shooting_omega, xmin, xmax, ymin, ymax, fname, maxpoints);
    
    return 0;
}
