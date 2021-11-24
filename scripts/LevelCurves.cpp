#include "LevelCurves.h"


double L2SquareDistance(double x1, double y1, double x2, double y2){
   return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2); 
}
bool TooClose(double a, double b){
    double eps = 1e-10;
    if (fabs(b-a)<eps){
        return 1;
    } else {
        return 0;
    }
}
int BisectionForLevelCurves(double x1, double y1, double x2, double y2, double (*func)(double,double), double K, double tol, 
              double &xroot, double &yroot){
    double tol2;
    int itermax, iter;
    double f1, f2;
    double xm, ym, fm;

    tol2    = tol*tol;
    itermax = 100;
    iter    = -1;

    f1 = func(x1,y1)-K;
    f2 = func(x2,y2)-K;

    if (f1*f2>0){
        return iter;
    }
    
    iter++;
    
    while (1) {
        xm = (x1+x2)/2;
        ym = (y1+y2)/2;

        if (L2SquareDistance(x1,y1,x2,y2)<tol2 || iter>itermax){
            if (iter>itermax) {
                printf("reached itermax=%d!\n", itermax);
            }
            xroot = xm;
            yroot = ym;
            return iter;
        }
        
        fm = func(xm,ym)-K;
        iter++;
        
        if (f1*fm<0){
            x2 = xm;
            y2 = ym;
        } else if (f1*fm>0){
            x1 = xm;
            y1 = ym;
            f1 = fm;
        } else {
            xroot = xm;
            yroot = ym;
            return iter;
        }
    }
}

void CreateSquare(double x0, double y0, double dx, double dy, double xmin, double xmax, double ymin, double ymax, 
                  double square[][2]){
    double eps = 1e-10;
    double x1, x2, y1, y2;
    x1  = x0-dx;
    x2  = x0+dx;
    y1  = y0-dy;
    y2  = y0+dy;
    if (x1<xmin){
        x1 = xmin+eps;
    }
    if (x2>xmax){
        x2 = xmax-eps;
    }
    if (y1<ymin){
        y1 = ymin+eps;
    }
    if (y2>ymax){
        y2 = ymax-eps;
    }
    square[0][0] = x1;
    square[0][1] = y1;
    square[1][0] = x1;
    square[1][1] = y2;
    square[2][0] = x2;
    square[2][1] = y2;
    square[3][0] = x2;
    square[3][1] = y1;
    return;
}

int FindSquareSide(double K, double x0, double y0, double square[][2], double x_previous, double y_previous, 
                   double (*func)(double, double)){
    double k1, k2, k3, k4;
    double intervals[4][2]={{0}};
    double logic_vec[4]={0};

    k1 = func(square[0][0], square[0][1]);
    k2 = func(square[1][0], square[1][1]);
    k3 = func(square[2][0], square[2][1]);
    k4 = func(square[3][0], square[3][1]);

    if (k1<k2){
        intervals[0][0] = k1;
        intervals[0][1] = k2;
    } else {
        intervals[0][0] = k2;
        intervals[0][1] = k1;
    }

    if (k2<k3){
        intervals[1][0] = k2;
        intervals[1][1] = k3;
    } else {
        intervals[1][0] = k3;
        intervals[1][1] = k2;
    }

    if (k3<k4){
        intervals[2][0] = k3;
        intervals[2][1] = k4;
    } else {
        intervals[2][0] = k4;
        intervals[2][1] = k3;
    }

    if (k4<k1){
        intervals[3][0] = k4;
        intervals[3][1] = k1;
    } else {
        intervals[3][0] = k1;
        intervals[3][1] = k4;
    }


    for(int i=0; i<4; i++){
        if (K>intervals[i][0] && K<intervals[i][1])
            logic_vec[i] = 1;
    }

    int side_number;

    if ((logic_vec[0]+logic_vec[1]+logic_vec[2]+logic_vec[3])<1e-12)
        side_number = 0;

    if (logic_vec[0] && logic_vec[2]){
        if (x0>x_previous){
            side_number = 3;
        } else {
            side_number = 1;
        }
    } else if (logic_vec[1] && logic_vec[3]) {
        if (y0>y_previous){
            side_number = 2;
        } else {
            side_number = 4;
        }
    } else if (logic_vec[0] && logic_vec[1]){
        if (y0>y_previous){ 
            side_number = 2;
        } else {
            side_number = 1;
        }
    } else if (logic_vec[1] && logic_vec[2]){
        if (y0>y_previous){
            side_number = 2;
        } else {
            side_number = 3;
        }
    } else if (logic_vec[2] && logic_vec[3]){
        if (y0>y_previous){
            side_number = 3;
        } else {
            side_number = 4;
        }
    } else if (logic_vec[0] && logic_vec[3]){
        if (y0>y_previous){
            side_number = 1;
        } else {
            side_number = 4;
        }
    }
    return side_number;
}

void FindPointsForRootFinder(double xp,  double yp,  double x0,  double y0, double L, 
                             double &x1, double &y1, double &x2, double &y2){
    if (xp==x0){
        x1 = x0;
        y1 = y0-L;
        x2 = x0;
        y2 = y0+L;
    } else {
        double m, xm, ym, M, Q, a, b, c;
        // angular coeff of the line passing for (xp,yp) and (x0,y0)
        m = (y0-yp)/(x0-xp);
        // middle point of (x1,y1) and (x2,y2)
        xm = 2*x0-xp;
        ym = 2*y0-yp;
        // line perpendicular to y=m*x and passing for the middle point
        M = -1/m;
        Q = ym - M*xm;
        // find (x1,y1) and (x2,y2)
        a  = (1+M*M);
        b  = 2*M*Q-2*xm-2*ym*M;
        c  = xm*xm+ym*ym+Q*Q-2*ym*Q-L*L;
        x1 = (-b-sqrt(b*b-4*a*c))/(2*a);
        y1 = M*x1+Q;
        x2 = (-b+sqrt(b*b-4*a*c))/(2*a);
        y2 = M*x2+Q;
    }

    return;
}

bool FindSecondPoint(double xp, double yp, double dx, double dy, int start_direction, 
                     double (*func)(double,double), double K, double tol, double &x0, double &y0){
    printf("Searching level-curve for M=%.3f starting from (x,y)=(%.3f,%.3f)\n", K, xp, yp);
    double a,b,c;
    bool isyaxis;
    // find second point 
    if (start_direction==1) {
        a       = yp-dy;
        b       = yp+dy;
        c       = xp-dx;
        isyaxis = 1;
    } else if (start_direction==2){
        a       = xp-dx;
        b       = xp+dx;
        c       = yp+dy;
        isyaxis = 0;
    } else if (start_direction==3){
        a       = yp-dy;
        b       = yp+dy;
        c       = xp+dx;
        isyaxis = 1;
    } else if (start_direction==4){
        a       = xp-dx;
        b       = xp+dx;
        c       = yp-dy;
        isyaxis = 0;
    } else {
        printf("start_direction must be 1,2,3 or 4");
        return 0;
    }

    if (isyaxis) {
        return BisectionForLevelCurves(c, a, c, b, func, K, tol, x0, y0);
    } else {
        return BisectionForLevelCurves(a, c, b, c, func, K, tol, x0, y0);
    }
}

int FindPoints(double x_init, double y_init, double dx, double dy, double x0, double y0,
                    double xmin, double xmax, double ymin, double ymax, double L, double tol, double **points, int N,
                    double (*func)(double,double), double K){
    double a,b,c;
    bool isyaxis;

    double x1_bis, y1_bis, x2_bis, y2_bis;
    double x0_tmp, y0_tmp;

    int found_root;

    double square[4][2];

    double x1,y1,x2,y2;
    

    double xp, yp;

    xp = x_init;
    yp = y_init;
    
    int side_number;
    int iter_bisec;

    // find second point 
    int points_counter = 2;;
    for(int i=2; i<N; i++){
        FindPointsForRootFinder(xp,yp,x0,y0,L,x1_bis,y1_bis,x2_bis,y2_bis);
        iter_bisec = BisectionForLevelCurves(x1_bis, y1_bis, x2_bis, y2_bis, func, K, tol, x0_tmp, y0_tmp);
        if (iter_bisec>0){
            xp  = x0;
            yp  = y0;
            x0 = x0_tmp;
            y0 = y0_tmp;
        } else {
            CreateSquare(x0,y0,dx,dy,xmin,xmax,ymin,ymax,square);
            x1 = square[0][0];
            x2 = square[3][0]; 
            y1 = square[0][1];
            y2 = square[1][1];
            side_number = FindSquareSide(K, x0, y0, square, xp, yp, func);
            if (side_number==1){
                a       = y1;
                b       = y2;
                c       = x1;
                isyaxis = 1;
            } else if (side_number==2) {
                a       = x1;
                b       = x2;
                c       = y2;
                isyaxis = 0;
            } else if (side_number==3) {
                a       = y1;
                b       = y2;
                c       = x2;
                isyaxis = 1;
            } else if (side_number==4) {
                a       = x1;
                b       = x2;
                c       = y1;
                isyaxis = 0;
            } else {
                printf("Side of the square not found! Try to change the initial direction\n");
            }

            xp = x0;
            yp = y0;

            if (isyaxis) {
                BisectionForLevelCurves(c, a, c, b, func, K, tol, x0, y0);
            } else {
                BisectionForLevelCurves(a, c, b, c, func, K, tol, x0, y0);
            }
        }
        
        points[i][0] = x0;
        points[i][1] = y0;
        points_counter++;

        if (TooClose(x0, xmin)   || TooClose(x0, xmax) || TooClose(y0, ymax)   || TooClose(y0, ymin) ||
           (fabs(x0-x_init)<dx && fabs(y0-y_init)<dy)){
            printf("Stopped at %d iteration\n", i);
            break;
        }
    }
    return points_counter; 
}

void FindLevelCurve(double x_init, double y_init, double dx, double dy, int start_direction, double L, double tol, double
                   (*func)(double, double), double xmin, double xmax, double ymin, double ymax, 
                   char *fname, int maxpoints){
    
    double x0, y0,K;
    K = func(x_init,y_init);
    FindSecondPoint(x_init, y_init, dx, dy, start_direction, func, K, tol, x0, y0);
    
    // allocate points
    double **points; 
    points = new double*[maxpoints]; 
    for (int i=0; i<maxpoints; i++)
        points[i] = new double[2];
    
    points[0][0] = x_init;
    points[0][1] = y_init;
    points[1][0] = x0;
    points[1][1] = y0;
    
    int iter;
    iter = FindPoints(x_init, y_init, dx, dy, x0, y0, xmin, xmax, ymin, ymax, L, tol, points, maxpoints, func, K);
    
    FILE *fp;
    fp=fopen(fname,"w");
    for(int i=0; i<iter; i++)
        fprintf(fp, "%23.15e %23.15e\n", points[i][0], points[i][1]);
    fclose(fp);

    cout<<"Ouput printed on '"<<fname<<"'"<<endl;
    
    // free memory
    for(int i=0 ;i<maxpoints; i++)
        delete[] points[i]; 
    delete[] points;
}


