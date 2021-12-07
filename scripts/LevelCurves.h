#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double L2SquareDistance(double,double,double,double);
int BisectionForLevelCurves(double x1, double y1, double x2, double y2, double (*func)(double,double), double K, double tol, 
              double &xroot, double &yroot);
void CreateSquare(double x0, double y0, double ds, double xmin, double xmax, double ymin, double ymax, 
                  double square[][2]);
int FindSquareSide(double K, double x0, double y0, double square[][2], double x_previous, double y_previous, 
                   double (*func)(double, double));
void FindPointsForRootFinder(double xp,  double yp,  double x0,  double y0, double L, 
                             double &x1, double &y1, double &x2, double &y2);
bool FindSecondPoint(double xp, double yp, double ds, int start_direction, 
                     double (*func)(double,double), double K, double tol, double &x0, double &y0);
int FindPoints(double x_init, double y_init, double ds, double x0, double y0, 
                    double xmin, double xmax, double ymin, double ymax, double L, double tol, double **points, int N,
                    double (*func)(double,double), double K);
void FindLevelCurve(double x_init, double y_init, double ds, int start_direction, double L, double tol, double
                   (*func)(double, double), double xmin, double xmax, double ymin, double ymax, 
                   char *fname, int maxpoints);

