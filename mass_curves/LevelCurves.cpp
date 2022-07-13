/*
 *   DESCRIPTION NEEDED
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

#include "LevelCurves.h"
/*
TODO
- add comments 
- fix notation: x1/x2 used many times for different things
*/


double L2SquareDistance(double x1, double y1, double x2, double y2){
   return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2); 
}
int BisectionForLevelCurves(double x1, double y1, double x2, double y2, 
                            long double (*func)(long double, long double), double K, double tol, 
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
    
    //cout<<"-------------------------------------------"<<endl;
    //printf("x1: %-8.5f y1: %-8.5f M1: %-8.5Lf\n", x1, y1, func(x1,y1));
    //printf("x2: %-8.5f y2: %-8.5f M2: %-8.5Lf\n", x2, y2, func(x2,y2));
    
    if (f1*f2>0){
        return iter;
    }
    
    iter++;
    
    while (1) {
        xm = (x1+x2)/2;
        ym = (y1+y2)/2;

        //cout << "xm: " << xm << " ym: " << ym << endl;

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
    /*
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
    */
    x1 = max(x1,xmin);
    x2 = min(x2,xmax);
    y1 = max(y1,ymin);
    y2 = min(y2,ymax);

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
                   long double (*func)(long double, long double)){
    double k1, k2, k3, k4;
    double intervals[4][2]={{0}};
    int logic_vec[4]={0};

    k1 = func(square[0][0], square[0][1]);
    k2 = func(square[1][0], square[1][1]);
    k3 = func(square[2][0], square[2][1]);
    k4 = func(square[3][0], square[3][1]);
    
    //cout<<"-------------------------------------------"<<endl;
    //printf("x1: %-8.5f y1: %-8.5f M1: %-8.5f\n", square[0][0], square[0][1], k1);
    //printf("x2: %-8.5f y2: %-8.5f M2: %-8.5f\n", square[1][0], square[1][1], k2);
    //printf("x3: %-8.5f y3: %-8.5f M3: %-8.5f\n", square[2][0], square[2][1], k3);
    //printf("x4: %-8.5f y4: %-8.5f M4: %-8.5f\n", square[3][0], square[3][1], k4);
    //cout<<"-------------------------------------------"<<endl;

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

    int side_number=0;
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

    // check single passage (possible at iter=0)
    if (side_number<1e-12){
        for(int i=0; i<4; i++){
            if (logic_vec[i]==1){
                side_number = i+1;
                break;
            }
        }
    }
    
    return side_number;
}

void FindPointsForRootFinder(double xp,  double yp,  double x0,  double y0, double L, 
                             double &x1, double &y1, double &x2, double &y2){
    if (fabs(xp-x0)<1.e-12){
        double dy = y0-yp;
        x1 = x0-L;
        y1 = y0+dy;
        x2 = x0+L;
        y2 = y0+dy;
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

int FindSecondPoint(double xp, double yp, double dx, double dy, int start_direction, 
                     double xmin, double xmax, double ymin, double ymax,
                     long double (*func)(long double, long double), double K, double tol, double &x0, double &y0, double &L){
    printf("Searching level-curve for K=%.5f starting from (x,y)=(%.5f,%.5f)\n", K, xp, yp);
    double a,b,c;
    double dX,dY,molt;
    bool isyaxis;
    int iter=-1; 
    
    molt = 1.;
    while (iter<0){
        dX    = dx*molt;
        dY    = dy*molt;
        if (start_direction==1) {
            a       = max(yp-dY,ymin);
            b       = min(yp+dY,ymax);
            c       = xp-dx;
            isyaxis = 1;
        } else if (start_direction==2){
            a       = max(xp-dX,xmin);
            b       = min(xp+dX,xmax);
            c       = yp+dy;
            isyaxis = 0;
        } else if (start_direction==3){
            a       = max(yp-dY,ymin);
            b       = min(yp+dY,ymax);
            c       = xp+dx;
            isyaxis = 1;
        } else if (start_direction==4){
            a       = max(xp-dX,xmin);
            b       = min(xp+dX,xmax);
            c       = yp-dy;
            isyaxis = 0;
        } else {
            printf("start_direction must be 1,2,3 or 4");
            return 0;
        }
        if (isyaxis) {
            iter = BisectionForLevelCurves(c, a, c, b, func, K, tol, x0, y0);
        } else {
            iter = BisectionForLevelCurves(a, c, b, c, func, K, tol, x0, y0);
        }
        L = fabs(b-a)/2;
        //printf("%f %f %f %f %f %f\n", dx, dy, dX, dY, molt, L);
        molt *= 1.5;
    }
    return iter;
}

int Iterate(double x_init, double y_init, double dx, double dy, double x0, double y0,
                    double xmin, double xmax, double ymin, double ymax, double L, double tol, double **points, int N,
                    long double (*func)(long double, long double), double K){
    double x1_bis, y1_bis, x2_bis, y2_bis;
    double x0_tmp, y0_tmp;

    // variables needed for the square
    double square[4][2];
    double x1,y1,x2,y2;
    int side_number;
    double a,b,c;
    bool isyaxis;

    double xp, yp;
    xp = x_init;
    yp = y_init;
   
    // Find points 
    int iter_bisec;
    int points_counter = 2;
    int last_point_on_boundary = 0;
    
    // variables used for checking that the points for bisection are into the boundary
    int boundary_iters, boundary_iters_max=9, boundary_iters_max_reached;
    double L_tmp, puppet, L_jj;
    
    for(int i=2; i<N; i++){
        if (!last_point_on_boundary){
            
            L_jj=L;
            
            //for ( int jj=1; jj<=2; jj++){
            
                FindPointsForRootFinder(xp,yp,x0,y0,L_jj,x1_bis,y1_bis,x2_bis,y2_bis);
           
                boundary_iters_max_reached = 0;
            
                // check if (x1_bis, y1_bis) is inside the boundary, otherwise change it
                L_tmp = L_jj;
                boundary_iters = 0;
                while (x1_bis<xmin || x1_bis>xmax || y1_bis<ymin || y1_bis>ymax){
                    L_tmp = L_tmp-L_jj/10;
                    FindPointsForRootFinder(xp,yp,x0,y0,L_tmp,x1_bis,y1_bis,puppet,puppet);
                    boundary_iters++;
                    if (boundary_iters>boundary_iters_max){
                        boundary_iters_max_reached = 1;
                        break;
                    }
                }
            
                // check if (x2_bis, y2_bis) is inside the boundary, otherwise change it
                L_tmp = L_jj;
                boundary_iters = 0;
                while (x2_bis<xmin || x2_bis>xmax || y2_bis<ymin || y2_bis>ymax){
                    L_tmp = L_tmp-L_jj/10;
                    FindPointsForRootFinder(xp,yp,x0,y0,L_tmp,puppet,puppet,x2_bis,y2_bis);
                    boundary_iters++;
                    if (boundary_iters>boundary_iters_max){
                        boundary_iters_max_reached = 1;
                        break;
                    }
                }
           
                if (!boundary_iters_max_reached) {
                    // just a check
                    if (x2_bis<xmin || x2_bis>xmax || y2_bis<ymin || y2_bis>ymax ||
                        x1_bis<xmin || x1_bis>xmax || y1_bis<ymin || y1_bis>ymax){
                        printf("+++ Warning! Initial points for bisection outside of boundary at point #%d ++\n", i); 
                    }
                    // bisection
                    iter_bisec = BisectionForLevelCurves(x1_bis, y1_bis, x2_bis, y2_bis, func, K, tol, x0_tmp, y0_tmp);
                } else {
                    iter_bisec = -1;
                }
                
               // L_jj*=2.;
               // if (iter_bisec >= 0) break;
           
           // }
        } else {
            iter_bisec = -1;
        }
        if (iter_bisec>=0){
            xp  = x0;
            yp  = y0;
            x0 = x0_tmp;
            y0 = y0_tmp;
        } else {
            printf("Using squares-method at iteration %d, (xp,yp)=(%9f,%9f)\n", i, xp, yp);
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
        printf("Points found: %d \t ( %f , %f )\n", points_counter,x0,y0);
        points_counter++;
        
        if (x0<xmin || x0>xmax || y0<ymin || y0>ymax || 
            last_point_on_boundary || (fabs(x0-x_init)<=dx && fabs(y0-y_init)<=dy)){
            printf("Stopped at %d iteration\n", i);
            break;
        }
        
        if ((x0-xmin)<dx || (xmax-x0)<dx || (ymax-y0)<dy || (y0-ymin)<dy)
            last_point_on_boundary++;
    }
    return points_counter; 
}

void FindLevelCurve(double x_init, double y_init, double ds, char *start_direction_str, double tol, long double
                   (*func)(long double, long double), double xmin, double xmax, double ymin, double ymax, 
                   char *fname, int maxpoints){
    double x0, y0, K;
    double dx, dy, L;
    int start_direction;
    int iter=0;

    K  = func(x_init,y_init);
    dx = ds*(xmax-xmin);
    dy = ds*(ymax-ymin);
   
    if (!strcmp(start_direction_str,"left")){
        start_direction = 1;
    } else if (!strcmp(start_direction_str, "up")){
        start_direction = 2;
    } else if (!strcmp(start_direction_str, "right")){
        start_direction = 3;
    } else if (!strcmp(start_direction_str, "down")){
        start_direction = 4;
    } else {
        printf("Wrong input! The initial direction can be 'left', 'up', 'right' or 'down'");
        return;
    }
    FindSecondPoint(x_init, y_init, dx, dy, start_direction, xmin, xmax, ymin, ymax, func, K, tol, x0, y0, L);
    
    // allocate points
    double **points; 
    points = new double*[maxpoints]; 
    for (int i=0; i<maxpoints; i++)
        points[i] = new double[2];
    
    points[0][0] = x_init;
    points[0][1] = y_init;
    points[1][0] = x0;
    points[1][1] = y0;

    iter = Iterate(x_init, y_init, dx, dy, x0, y0, xmin, xmax, ymin, ymax, L, tol, points, maxpoints, func, K);
    
    FILE *fp;
    fp=fopen(fname,"w");
    
    // output for fermion-axion stars
    fprintf(fp, "-10 10000\n");
    fprintf(fp, " %d \n\n\n", iter);
    for(int i=0; i<iter; i++)
        fprintf(fp, "%23.15e \t %23.15e \t 100 \n", points[i][1], points[i][0]);
    
    // standard output 
    //for(int i=0; i<iter; i++)
    //    fprintf(fp, "%23.15e \t %23.15e \n", points[i][0], points[i][1]);
    
    fclose(fp);

    cout<<"Ouput printed on '"<<fname<<"'"<<endl;
    
    // free memory
    for(int i=0 ;i<maxpoints; i++)
        delete[] points[i]; 
    delete[] points;
}
