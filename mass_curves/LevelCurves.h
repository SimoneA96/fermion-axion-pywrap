#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace std;

double L2SquareDistance(double,double,double,double);
int BisectionForLevelCurves(double x1, double y1, double x2, double y2, long double (*func)(long double, long double), double K, double tol, 
              double &xroot, double &yroot);
void CreateSquare(double x0, double y0, double dx, double dy, double xmin, double xmax, double ymin, double ymax, 
                  double square[][2]);
int FindSquareSide(double K, double x0, double y0, double square[][2], double x_previous, double y_previous, 
                   double (*func)(double, double));
void FindPointsForRootFinder(double xp,  double yp,  double x0,  double y0, double L, 
                             double &x1, double &y1, double &x2, double &y2);
int FindSecondPoint(double xp, double yp, double dx, double dy, int start_direction, 
                     double xmin, double xmax, double ymin, double ymax,
                     long double (*func)(long double, long double), double K, double tol, double &x0, double &y0, double &L);
int Iterate(double x_init, double y_init, double dx, double dy, double x0, double y0, 
                    double xmin, double xmax, double ymin, double ymax, double L, double tol, double **points, int N,
                    long double (*func)(long double, long double), double K);
void FindLevelCurve(double x_init, double y_init, double ds, char *start_direction, double tol, 
                    long double (*func)(long double, long double), double xmin, double xmax, 
                    double ymin, double ymax, char *fname, int maxpoints);

