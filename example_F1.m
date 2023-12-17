% Example of injection induces stresses around a fault
% Stresses are computed from inclusion theory
close all
clear variables

%======================
% Reservoir geometry
%======================
width = 20000;  % total width of reservoir (m)
height = 225;  % reservoir height (m)
throw = 75;  % fault throw (m)
dip = 70;  % dip angle (degrees)
theta = dip/180*pi;  % dip angle (rad)
%======================
% Parameters
%======================
alpha = 0.9;  % Biot's coefficient (-)
nu = 0.15;  % Poisson's ratio (-)
mu = 6500;  % Shear modulus (MPa)
lambda = 2*mu*nu/(1-2*nu);  % Lamé's first parameter (MPa)
delta_p = 1;  % uniform pore pressure change in reservoir (MPa)
              % positive delta_p = injection
              % negative delta_p = depletion
D = (1-2*nu)*alpha*delta_p/(2*pi*(1-nu)*mu);  % scaling coefficient for strains and displacements
%======================
% Grid
%======================
% Define x- and y-coordinates of grid
xmin = 9700;
xmax = 10300;
ymin = -250;
ymax = 250;
resolution_x = 200;  % number of cells in x-direction
resolution_y = 200;  % number of cells in y-direction
% Create evenly spaced grid points
[x,y] = meshgrid(linspace(xmin, xmax, resolution_x), linspace(ymin, ymax, resolution_y));

%======================
% Analytical solution
%======================
% Reservoir consists of two rectangles + two triangles

% *** Rectangle 1 ***
p = 0;  % minimum x-coordinate of rectangle 1
q = 0.5*(width - (height + throw)/tan(theta));  % maximum x-coordinate of rectangle 1
r = -0.5*(height + throw);  % minimum y-coordinate of rectangle 1
s = r + height;  % maximum y-coordinate of rectangle 1
[Gxx1, Gyy1, Gxy1] = strain_rectangle(x, y, p, q, r, s);


% *** Triangle 1 *** 
% Shape:
%  ___
% |  /
% | /
% |/
%
x0 = 0.5*(width - (height + throw)/tan(theta));  % minimum x-coordinate of triangle 1
x1 = x0 + height/tan(theta);  % maximum x-coordinate of triangle 1
y0 = -0.5*(height + throw);  % minimum y-coordinate of triangle 1
y1 = y0 + height;  % maximum y-coordinate of triangle 1
% Note that we shift the x- and y-coordinates such that the
% hypotenuse of the triangle goes through the origin of the coordinate
% system. % Additionally, this triangle is upside down compared to the reference
% case. % Hence, we must use f(x,y,p,o,s,r) instead of f(x,y,o,p,r,s) to get the correct
% solution. 
[Gxx2, Gyy2, Gxy2] = strain_triangle(x - x0, y - y0, x1 - x0, 0, y1 - y0, 0, theta);

% *** Triangle 2 *** 
% Shape: 
% 
%   /|
%  / |
% /__|
%
x0 = 0.5*(width + (height + throw)/tan(theta)) - height/tan(theta); % minimum x-coordinate of triangle 2
x1 = 0.5*(width + (height + throw)/tan(theta));  % maximum x-coordinate of triangle 2
y0 = -0.5*(height - throw);  % minimum y-coordinate of triangle 2
y1 = 0.5*(height + throw);  % maximum y-coordinate of triangle 2
% Note that we shift the x- and y-coordinates such that the
% hypotenuse of the triangle goes through the origin of the coordinate
% system
[Gxx3, Gyy3, Gxy3] = strain_triangle(x - x0, y - y0, 0, x1 - x0, 0, y1 - y0, theta);


% *** Rectangle 2 *** 
p = 0.5*(width + (height + throw)/tan(theta));  % minimum x-coordinate of rectangle 2
q = width;  % maximum x-coordinate of rectangle 2
r = -0.5*(height - throw);  % minimum y-coordinate of rectangle 2
s = 0.5*(height + throw);  % maximum y-coordinate of rectangle 2
[Gxx4, Gyy4, Gxy4] = strain_rectangle(x, y, p, q, r, s);

% *** Total strains = sum of all rectangles and triangles *** 
exx = 0.5*D*(Gxx1 + Gxx2 + Gxx3 + Gxx4);  % Horizontal strain
eyy = 0.5*D*(Gyy1 + Gyy2 + Gyy3 + Gyy4);  % Vertical strain
exy = 0.5*D*(Gxy1 + Gxy2 + Gxy3 + Gxy4);  % Shear strain

% Find the grid cells inside the inclusion
% Here, we simply use the trick that the scaled volumetric strain must be larger
% than zero (since volumetric strain is zero outside the inclusion)
krDelta = abs(exx + eyy)/D > 1;

% *** Stresses ***
sigma_xx = lambda*(exx + eyy) + 2*mu*exx - alpha*krDelta*delta_p;
sigma_yy = lambda*(exx + eyy) + 2*mu*eyy - alpha*krDelta*delta_p;
sigma_xy = 2*mu*exy;

% ======================
% Plot stresses
% ======================
subplot(1,3,1), surf(x,y,sigma_xx), view(0,90), shading interp, colormap jet;
xlabel('x (m)'), ylabel('y (m)')
set(gca,'FontSize',14)
cbar = colorbar;
cbar.Label.String = '$\Delta \sigma_{xx}$';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 16;
xlim([xmin, xmax]), ylim([ymin, ymax])
subplot(1,3,2), surf(x,y,sigma_yy), view(0,90), shading interp;
xlabel('x (m)'), ylabel('y (m)')
set(gca,'FontSize',14)
cbar = colorbar;
cbar.Label.String = '$\Delta \sigma_{yy}$';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 16;
xlim([xmin, xmax]), ylim([ymin, ymax])
subplot(1,3,3), surf(x,y,sigma_xy), view(0,90), shading interp;
xlabel('x (m)'), ylabel('y (m)')
set(gca,'FontSize',14)
cbar = colorbar;
cbar.Label.String = '$\Delta \sigma_{xy}$';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 16;
xlim([xmin, xmax]), ylim([ymin, ymax])




function [Gxx, Gyy, Gxy] = strain_rectangle(x, y, p, q, r, s)
% Analytical solution for the strains due to a uniform pore pressure change
% inside a rectangular inclusion

% Scaled horizontal strain
Gxx = atan((y-s)./(x-q)) - atan((y-s)./(x-p)) - atan((y-r)./(x-q)) + atan((y-r)./(x-p));
% Scaled vertical strain
Gyy = atan((x-q)./(y-s)) - atan((x-p)./(y-s)) - atan((x-q)./(y-r)) + atan((x-p)./(y-r));
% Scaled shear strain
Gxy = 0.5*log(((x-q).^2+(y-s).^2).*((x-p).^2+(y-r).^2)./(((x-q).^2+(y-r).^2).*((x-p).^2+(y-s).^2)));
end

function [Gxx,Gyy,Gxy] = strain_triangle(x,y,o,p,r,s,theta)
% Analytical solution for the strains due to a uniform pore pressure change
% inside a triangular inclusion

% Scaled horizontal strain
Gxx = 0.5*log(((y-r).^2+(x-r*cot(theta)).^2)./((y-s).^2+(x-s*cot(theta)).^2))*sin(theta)*cos(theta) ...
    + atan2((r-s).*(x-p),(x-p).^2+(y-r).*(y-s)) ...
    - atan2((s-r).*(y*cot(theta)-x),(x.^2+y.^2+r.*s*csc(theta)^2 - (r+s).*(y+x*cot(theta))))*sin(theta)^2;
% Scaled vertical strain
Gyy = 0.5*sin(theta)*cos(theta)*log(((x-p).^2+(y-p*tan(theta)).^2)./((x-o).^2+(y-o*tan(theta)).^2)) ...
    - atan2((o-p).*(y-r),(y-r).^2+(x-p).*(x-o)) ...
    - atan2((p-o).*(y-x*tan(theta)),(x.^2+y.^2+o.*p*sec(theta)^2-(p+o).*(x+y*tan(theta))))*cos(theta)^2;
% Scaled shear strain
Gxy = 0.5*log(((x-p).^2+(y-s).^2)./((x-p).^2+(y-r).^2)) + 0.5*log(((y-r).^2+(x-r*cot(theta)).^2)./((y-s).^2+(x-s*cot(theta)).^2))*sin(theta)^2 ...
    + atan2((s-r).*(y*cot(theta)-x),(x.^2+y.^2+r.*s*csc(theta)^2-(r+s).*(y+x*cot(theta))))*sin(theta)*cos(theta);
end


