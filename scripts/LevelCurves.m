%function LevelCurves(DEBUG, k, x0, y0)
function LevelCurves(DEBUG,N)

xmin = -1;
xmax =  1;
nx   = 500;
ymin = -1;
ymax =  1;
ny   = 500;

x1d = linspace(xmin,xmax,nx);
y1d = linspace(ymin,ymax,ny);

[x,y] = meshgrid(x1d,y1d);
z     = func(x,y);

if DEBUG
    dxs = 0.08;
    dys = 0.08;

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
        zp = func(xp,yp);

        square = CreateSquare(x0, y0, dxs, dys, xmin, xmax, ymin, ymax);
        xs = square(:,1);
        ys = square(:,2);
        zs = func(xs,ys);
        x1 = square(1,1);
        x2 = square(4,1); 
        y1 = square(1,2);
        y2 = square(2,2);
        side_number = FindSquareSide(z0, x0, y0, square, xp, yp);
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
        scatter3(xs,ys,zs, 'x', 'MarkerEdgeColor', 'r')
        scatter3(x0,y0,z0, 'o', 'filled', 'MarkerFaceColor', 'g')
        scatter3(xp,yp,zp, 'o', 'filled', 'MarkerFaceColor', 'r')
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
end
points = zeros(N+2,2);
% initial point 
xp  = -1;      % previous point
yp  = 0.08216;
xsc = -0.9599; % center of square
ysc = 0.08617;
dxs = 0.03;
dys = 0.03;
M   = func(xp,yp);
points(1,:) = [xp, yp ];
points(2,:) = [xsc,ysc];

 
figure
contour(x,y,z, 50)
hold on 
scatter(xp,  yp , 'o', 'filled', 'MarkerFaceColor', 'r')
scatter(xsc, ysc, 'o', 'filled', 'MarkerFaceColor', 'g')

for i=1:N
    square = CreateSquare(xsc, ysc, dxs, dys, xmin, xmax, ymin, ymax);
    side_number = FindSquareSide(M, xsc, ysc, square, xp, yp);
    if side_number==1 || side_number==3
        a         = square(1,2);
        b         = square(2,2);
        isyaxis   = 1;
    else
        a         = square(1,1);
        b         = square(4,1); 
        isyaxis   = 0;
    end
    x1 = square(1,1);
    x2 = square(4,1); 
    y1 = square(1,2);
    y2 = square(2,2);
    not_found = 0;
    switch side_number
        case 0
            not_found = 1;
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
    end
        
    if not_found
        error('problem')
    else
        xp = xsc;
        yp = ysc;
        if isyaxis
            f   = @(y)(func(c,y)-M);
            xsc = c;
            ysc = fzero(f, (a+b)/2);
        else
            f   = @(x)(func(x,c)-M);
            xsc = fzero(f, (a+b)/2);
            ysc = c;
        end
    end
    
    points(i+2,:) = [xsc, ysc];
    
    scatter(xsc,ysc, 'x', 'b')
    xlim([-1, 0])
    ylim([0 1])
    pause(0.1)
    
    if tooclose(xsc, xmin) || tooclose(xsc, xmax) || tooclose(ysc, ymax) || tooclose(ysc,ymin)
        break
    end
end

return 

function bool = tooclose(a,b)
eps = 1e-3;
if abs(a-b)<eps
    bool = 1;
else
    bool = 0;
end
return

function z = func(x,y)
z = sin(x.*y);
return

function points = CreateSquare(x0, y0, dx, dy, xmin, xmax, ymin, ymax)
eps = 1e-10;
x1  = x0-dx/2;
x2  = x0+dx/2;
y1  = y0-dy/2;
y2  = y0+dy/2;
if x1<xmin
    x1 = xmin+eps;
end
if x2>xmax
    x2 = xmax-eps;
end
if y1<ymin
    y1 = ymin+eps;
end
if y2>ymax
    y2 = ymax-eps;
end
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

if all(logic_vec==0)
    side_number = 0;
end

if logic_vec(1) && logic_vec(3)
    if x0>x_previous
        side_number = 3;
    else
        side_number = 1;
    end
end

if logic_vec(2) && logic_vec(4)
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
return
