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
    char start_direction[] = "left";
    
    //x_init          = 0.5; // curve that ends on boundary
    //y_init          = 0.5;
    //start_direction = 1;
    x_init          = 0.3; // closed curve
    y_init          = -1;
    ds              = 0.02;  // used for initial step and for squares
    xmin            = -1;
    xmax            = 1;
    ymin            = -1;
    ymax            = 1;
    tol             = 1e-10; // tolerance in bisection
    maxpoints       = 100;   // maximum number of points for the level curve 

    // Find the level curve and prin the ouput on "output.txt"
    FindLevelCurve(x_init, y_init, ds, start_direction, tol, MyFunction, xmin, xmax, ymin, ymax, fname, maxpoints);
    
    return 0;
}

long double MyFunction(long double x, long double y){
    return sin(4*x)*cos(4*y);
}
