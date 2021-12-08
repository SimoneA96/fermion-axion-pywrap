/*
 *      move to the previous directory to lunch the code
 *
 */


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
    double x_init, y_init, ds, L, tol;
    int maxpoints, verbose;
    char fname[]="out_example.txt";
    char start_direction[]="up";
    
    x_init          = 0.002;
    y_init          = 0.01;
    ds              = 1e-2;  // used for initial step and for squares
    xmin            = 0;
    xmax            = 0.008;
    ymin            = 0;
    ymax            = 0.42;
    tol             = 1e-4; // tolerance in bisection
    maxpoints       = 20;   // maximum number of points for the level curve 

    // Find the level curve and prin the ouput on "output.txt"
    FindLevelCurve(x_init, y_init, ds, start_direction, tol, shooting_omega, xmin, xmax, ymin, ymax, fname, maxpoints);
    
    return 0;
}
