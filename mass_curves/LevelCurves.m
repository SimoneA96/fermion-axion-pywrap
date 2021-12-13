function LevelCurves(x_init, y_init, start_direction, ds, N, onlySquare, saveGIF)
%
%
% Method of perpendicular: 
% see figure that doesn't exist
%
%
% Square construction:
% Name convention for for vertices and sides. The starting direction is 
% perpendicular to the sides and pointing outside the square (e.g.
% start_direction=3 means --->, start_direction=1 means <---).
% 
%                 2
%         2---------------3
%         |               | 
%         |               |
%       1 |               | 3
%         |               |
%         |               |
%         1---------------4
%                 4
%
%
tol         = 1e-10;
DEBUG       = 0;
ShowTestFig = 1;
gifname     = 'curve.gif';
delay_time  = 0.1;

xmin = -1;
xmax =  0.2;
nx   = 500;
ymin = -0.6;
ymax =  0.6;
ny   = 500;

dx = ds*(xmax-xmin);
dy = ds*(ymax-ymin);

x1d = linspace(xmin,xmax,nx);
y1d = linspace(ymin,ymax,ny);

[x,y] = meshgrid(x1d,y1d);
z     = func(x,y);

if DEBUG
    dx = 0.08;
    dy = 0.08;

    xp_vec = [-0.4949, -0.798 , 0.4949,  0.2323];
    yp_vec = [-0.2323,  0.6768, 0.7172, -0.3535];
    x0_vec = [-0.3131, -0.6768, 0.5758,  0.3939];
    y0_vec = [-0.3737,  0.7374, 0.6162, -0.2121];

    for i=1:length(x0_vec)
        x0 = x0_vec(i);
        y0 = y0_vec(i);
        xp = xp_vec(i);
        yp = yp_vec(i);

        z0 = func(x0, y0);

        square = CreateSquare(x0, y0, dx, dy, xmin, xmax, ymin, ymax);
        xs = square(:,1);
        ys = square(:,2);
        x1 = square(1,1);
        x2 = square(4,1); 
        y1 = square(1,2);
        y2 = square(2,2);
        side_number = FindSquareSide(z0, x0, y0, square, xp, yp);
        if side_number==0
            disp(logic_vec)
            error('side not found!')
        end
        if side_number==1 || side_number==3
            a         = y1;
            b         = y2;
            isyaxis   = 1;
        else
            a         = x1;
            b         = x2; 
            isyaxis   = 0;
        end
        figure
        contour(x,y,z, 50)
        hold on 
        scatter(xs,ys, 'x', 'MarkerEdgeColor', 'r')
        scatter(x0,y0, 'o', 'filled', 'MarkerFaceColor', 'g')
        scatter(xp,yp, 'o', 'filled', 'MarkerFaceColor', 'r')
        
        [xc,yc,xd,yd] = FindPointsForRootFinder(xp,yp,x0,y0,0.1);
        xb = 2*x0-xp;
        yb = 2*y0-yp;
        scatter(xc,yc, 'filled', 'MarkerFaceColor', 'b')
        scatter(xd,yd, 'filled', 'MarkerFaceColor', 'b')
        scatter(xb,yb, 'filled', 'MarkerFaceColor', 'c')
        plot([xp,xb], [yp,yb], 'r')
        plot([xc,xd], [yc,yd], 'b')
        pbaspect([1 1 1])
        disp('square:')
        disp(square)
        fprintf('a           : %f\n', a)
        fprintf('b           : %f\n', b)
        fprintf('isyaxis     : %d\n', isyaxis)
        fprintf('side_number : %d\n', side_number)
        disp('---------------------------------------------------')
        disp('Press Enter to go to the next point...')
        pause
    end
    return
end

tic 

xp = x_init;
yp = y_init;
M  = func(xp,yp);
root_function = @(x,y) (func(x,y)-M);

fprintf('Searching level-curve for M=%.3f starting from (x,y)=(%.3f,%.3f)\n', M, xp, yp)

% find second point 
iter = 0;
molt = 1;
while iter<1
    dX = dx*molt;
    dY = dy*molt;
    switch start_direction
        case 1
            a       = max(yp-dY,ymin);
            b       = min(yp+dY,ymax);
            c       = xp-dx;
            isyaxis = 1;
        case 2
            a       = max(xp-dX,xmin);
            b       = min(xp+dX,xmax);
            c       = yp+dy;
            isyaxis = 0;
        case 3
            a       = max(yp-dY,ymin);
            b       = min(yp+dY,ymax);
            c       = xp+dx;
            isyaxis = 1;
        case 4
            a       = max(xp-dX,xmin);
            b       = min(xp+dX,xmax);
            c       = yp-dy;
            isyaxis = 0;
        otherwise
            error('start_direction must be 1,2,3 or 4')
    end
    
    if isyaxis
        [xsc, ysc, iter] = BisectionAlongLine(c, a, c, b, tol, root_function);
    else
        [xsc, ysc, iter] = BisectionAlongLine(a, c, b, c, tol, root_function);
    end
    molt = molt*1.5;
    L    = abs(b-a)/2;
end

points      = zeros(N+2,2);
points(1,:) = [xp, yp ];
points(2,:) = [xsc,ysc];

% squares method
if ShowTestFig
    fig = figure;
    set(gcf,'color','w');
    %set(gcf, 'Renderer', 'painters', 'Position', [78 10 700 500])
    contour(x,y,z, 50)
    hold on 
    scatter(xp,  yp, 80, 'o', 'filled', 'MarkerFaceColor', 'g')
    %xline(ymin, 'HandleVisibility', 'off');
    %xline(ymax, 'HandleVisibility', 'off');
    %yline(xmin, 'HandleVisibility', 'off');
    %yline(xmax, 'HandleVisibility', 'off');
    xlim([xmin xmax])
    ylim([ymin ymax])
    drawnow
    if saveGIF
        SavePhotogram(fig, gifname, 1,  delay_time*3)
    end
    scatter(xsc,ysc, 80, 'x', 'g', 'LineWidth', 2 )
    drawnow
    if saveGIF
        SavePhotogram(fig, gifname, 2,  delay_time)
    end
end

iter = 0;
last_point_on_boundary = 0;
boundary_iters_max = 9;
for i=1:N
    
    if ShowTestFig
        xb = 2*xsc-xp;
        yb = 2*ysc-yp;
        xp_old = xp;
        yp_old = yp;
        xsc_old = xsc;
        ysc_old = ysc;
    end
    
    if ~onlySquare && ~last_point_on_boundary
        [x1_bis,y1_bis,x2_bis,y2_bis] = FindPointsForRootFinder(xp,yp,xsc,ysc,L);
        
        boundary_iters_max_reached = 0;
        
        boundary_iters=0;
        L_tmp = L;
        while (x1_bis<xmin || x1_bis>xmax || y1_bis<ymin || y1_bis>ymax)
            L_tmp = L_tmp-L/10;
            [x1_bis,y1_bis] = FindPointsForRootFinder(xp,yp,xsc,ysc,L_tmp);
            boundary_iters = boundary_iters+1;
            if boundary_iters>boundary_iters_max
                boundary_iters_max_reached = 1;
                break;
            end
        end
        
        boundary_iters=0;
        L_tmp = L;
        while (x2_bis<xmin || x2_bis>xmax || y2_bis<ymin || y2_bis>ymax)
            L_tmp = L_tmp-L/10;
            [~, ~, x2_bis,y2_bis] = FindPointsForRootFinder(xp,yp,xsc,ysc,L_tmp);
            boundary_iters = boundary_iters+1;
            if boundary_iters>boundary_iters_max
                boundary_iters_max_reached = 1;
                break;
            end
        end
        
        if ~boundary_iters_max_reached
            [xsc_tmp,ysc_tmp,iter] = BisectionAlongLine(x1_bis, y1_bis, x2_bis, y2_bis, tol, root_function);
        else
            iter = 0;
        end
        %fprintf('%d iterations in BisectionAlongLine\n', iter)
    end
    
    if last_point_on_boundary
        iter = 0;
    end
    
    if iter>0
        xp  = xsc;
        yp  = ysc;
        xsc = xsc_tmp;
        ysc = ysc_tmp;
    else
        square = CreateSquare(xsc, ysc, dx, dy, xmin, xmax, ymin, ymax);
        x1 = square(1,1);
        x2 = square(4,1); 
        y1 = square(1,2);
        y2 = square(2,2);
        %scatter(x1,y1)
        %scatter(x1,y2)
        %scatter(x2,y2)
        %scatter(x2,y1)
        [side_number, logic_vec] = FindSquareSide(M, xsc, ysc, square, xp, yp);
        if side_number==0
            disp(logic_vec)
            error('side not found!')
        end
        switch side_number
            case 1
                a       = y1;
                b       = y2;
                c       = x1;
                isyaxis = 1;
            case 2
                a       = x1;
                b       = x2;
                c       = y2;
                isyaxis = 0;
            case 3
                a       = y1;
                b       = y2;
                c       = x2;
                isyaxis = 1;
            case 4
                a       = x1;
                b       = x2;
                c       = y1;
                isyaxis = 0;
            otherwise 
                error('Side of the square not found! Try to change the initial direction')
        end

        xp = xsc;
        yp = ysc;

        if isyaxis
            [xsc,ysc] = BisectionAlongLine(c, a, c, b, tol, root_function);
        else
            [xsc,ysc] = BisectionAlongLine(a, c, b, c, tol, root_function);
        end
    end
    
    points(i+2,:) = [xsc, ysc];
    
    if ShowTestFig
        if iter==0
            scatter(xsc_old,ysc_old, 70, 'o', 'MarkerEdgeColor', 'k')
            scatter(square(:,1), square(:,2), 'x', 'MarkerEdgeColor', 'r')
        else
            scatter(x1_bis,y1_bis, 20, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
            scatter(x2_bis,y2_bis, 20, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
            scatter(xb,yb, 'filled', 'MarkerFaceColor', 'c')
            plot([xp_old,xb], [yp_old,yb], 'r')
            plot([x1_bis,x2_bis], [y1_bis,y2_bis], 'b')
        end
        scatter(xsc,ysc, 80, 'x', 'g', 'LineWidth', 2 )
        xlim([xmin xmax])
        ylim([ymin ymax])
        drawnow
        
        if saveGIF
            SavePhotogram(fig, gifname, i+2, delay_time)
        end
    end
    
    if (xsc<xmin) || (xsc>xmax) || (ysc<xmin) || (ysc>ymax) || ...
       last_point_on_boundary || (abs(xsc-x_init)<=dx && abs(ysc-y_init)<=dy)
        fprintf('Stopped at %d iterations\n', i+1)
        break
    end
    
    if (xsc-xmin)<dx || (xmax-xsc)<dx || (ymax-ysc)<dy || (ysc-ymin)<dy
        last_point_on_boundary = last_point_on_boundary + 1;
    end
end

toc 


return 

function z = func(x,y)
z = sin(4*x).*cos(4*y);
%z = x;
%{
xmin = -1;
xmax =  1;
ymin = -1;
ymax =  1;
if x<xmin
    error('x=%f<xmin\n',x);
elseif x>xmax
    error('x=%f>xmax\n',x);
elseif y<ymin
    error('y=%f<ymin\n',y);
elseif y>ymax
    error('y=%f>ymax\n',y);
end
%}
return

function [xroot, yroot, iter] = BisectionAlongLine(x1, y1, x2, y2, tol, myfunc)
verbose = 0;
itermax = 100;
tol2    = tol*tol; % avoid sqrt() in if-statement

f1 = myfunc(x1,y1);
f2 = myfunc(x2,y2);

iter    = 0;
if f1*f2>0
    if verbose
        fprintf('root not found! Returning (xroot,yroot)=(%f,%f)\n',x1,y1)
    end
    xroot = x1;
    yroot = y1;
    return
end

while 1
    
    xm = (x1+x2)/2;
    ym = (y1+y2)/2;
    
    if ((x1-x2)^2+(y1-y2)^2)<tol2 || iter>itermax
        if iter>itermax
            warning('reached itermax=%d!', itermax)
        end
        xroot = xm;
        yroot = ym;
        if verbose
            fprintf('tol=%e\titer=%3d\t(%7.4f,%7.4f)\n',tol,iter,xm,ym)
        end
        return
    end
    
    iter = iter + 1;
    
    fm = myfunc(xm,ym);
    
    if f1*fm<0
        x2 = xm;
        y2 = ym;
    elseif f1*fm>0
        x1 = xm;
        y1 = ym;
        f1 = fm;
    else
        xroot = xm;
        yroot = ym;
        return
    end
end
return

function points = CreateSquare(x0, y0, dx, dy, xmin, xmax, ymin, ymax)
x1  = x0-dx;
x2  = x0+dx;
y1  = y0-dy;
y2  = y0+dy;

x1 = max(xmin,x1);
x2 = min(xmax,x2);
y1 = max(ymin,y1);
y2 = min(ymax,y2);

points      = zeros(4,2);
points(1,:) = [x1, y1];
points(2,:) = [x1, y2];
points(3,:) = [x2, y2];
points(4,:) = [x2, y1];
return

function [side_number, logic_vec] = FindSquareSide(k, x0, y0, square, x_previous, y_previous)
k1 = func(square(1,1), square(1,2));
k2 = func(square(2,1), square(2,2));
k3 = func(square(3,1), square(3,2));
k4 = func(square(4,1), square(4,2));

intervals      = zeros(4,2);
intervals(1,:) = sort([k1,k2]);
intervals(2,:) = sort([k2,k3]);
intervals(3,:) = sort([k3,k4]);
intervals(4,:) = sort([k4,k1]);
logic_vec      = zeros(4,1);
for i=1:4
    if k>intervals(i,1) && k<intervals(i,2)
        logic_vec(i) = 1;
    end
end

side_number = 0;
if logic_vec(1) && logic_vec(3)
    if x0>x_previous
        side_number = 3;
    else
        side_number = 1;
    end
elseif logic_vec(2) && logic_vec(4)
    if y0>y_previous
        side_number = 2;
    else
        side_number = 4;
    end
elseif logic_vec(1) && logic_vec(2)
    if y0>y_previous
        side_number = 2;
    else
        side_number = 1;
    end
elseif logic_vec(2) && logic_vec(3)
    if y0>y_previous
        side_number = 2;
    else
        side_number = 3;
    end
elseif logic_vec(3) && logic_vec(4)
    if y0>y_previous
        side_number = 3;
    else
        side_number = 4;
    end
elseif logic_vec(1) && logic_vec(4)
    if y0>y_previous
        side_number = 1;
    else
        side_number = 4;
    end
end

if sum(logic_vec)==1
    [~, side_number] = max(logic_vec);
end

return

function [x1,y1,x2,y2] = FindPointsForRootFinder(xp,yp,x0,y0,L)
if xp==x0
    % vertical line
    dy = y0-yp;
    x1 = x0-L;
    y1 = y0+dy;
    x2 = x0+L;
    y2 = y0+dy;
else
    % angular coeff of the line passing for (xp,yp) and (x0,y0)
    m = (y0-yp)/(x0-xp);
    %q = y0-m*x0;

    % middle point of (x1,y1) and (x2,y2)
    xm = 2*x0-xp;
    ym = 2*y0-yp;

    % line perpendicular to y=m*x and passing for the middle point
    M = -1/m;
    Q = ym - M*xm;
    
    % find (x1,y1) and (x2,y2)
    a  = (1+M^2);
    b  = 2*M*Q-2*xm-2*ym*M;
    c  = xm^2+ym^2+Q^2-2*ym*Q-L^2;
    x1 = (-b-sqrt(b^2-4*a*c))/(2*a);
    y1 = M*x1+Q;
    x2 = (-b+sqrt(b^2-4*a*c))/(2*a);
    y2 = M*x2+Q;
end
return

function SavePhotogram(fig, filename, iter, delay_time)
drawnow
frame      = getframe(fig); 
im         = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
if iter == 1 
  imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', delay_time); 
else 
  imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', delay_time); 
end 
return


