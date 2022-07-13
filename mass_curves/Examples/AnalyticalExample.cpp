/*
 *	move to the previous directory to run the code
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "LevelCurves.h"

/*
 Name convention for vertices and sides. The starting direction is 
 perpendicular to the sides and pointing outside the square (e.g.
 start_direction=3 means --->, start_direction=1 means <---).
 
                 2
         2---------------3
         |               | 
         |               |
       1 |               | 3
         |               |
         |               |
         1---------------4
                 4
*/

using namespace std;

long double MyFunction(long double, long double);

int main(int argc, char *argv[]){
    double xmin, xmax, ymin, ymax;
    double x_init, y_init, ds, tol;
    int maxpoints, verbose;
    char fname[]="out_example.txt";
    char start_direction[] = "up";
    
    x_init          = -0.6; // closed curve
    y_init          = -0.30103;
    ds              = 0.01;  // used for initial step and for squares
    xmin            = -1;
    xmax            = 0.2;
    ymin            = -0.6;
    ymax            = 0.6;
    tol             = 1e-10; // tolerance in bisection
    maxpoints       = 500;   // maximum number of points for the level curve 

    // Find the level curve and prin the ouput on "out_example.txt"
    FindLevelCurve(x_init, y_init, ds, start_direction, tol, MyFunction, xmin, xmax, ymin, ymax, fname, maxpoints);
    
    return 0;
}

long double MyFunction(long double x, long double y){
    return sin(4*x)*cos(4*y);
}
