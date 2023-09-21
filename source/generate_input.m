function [alpha,xsize, ysize, Nx, Ny, rho, rhof,rhog, Vp, cs, nu, G, lambda, eta, g, K0, Vw, Vi, mu0, V0, a0 ,b0, L]=generate_input()

%fault parameters
alpha=70;                               % Dip of the fault

%grid parameters
xsize=2000;                             % Horizontal model size, m
ysize=2000;                             % Vertical model size, m
Nx=201;                                 % Horizontal model resolution (number of points)
Ny=201;                                 % Vertical model resolution (number of points)

%material parameters
rho=2400;                               %Density of the rock [kg/m3]
rhof=1150;                              %Density of the fuid [kg/m3]
rhog=200;                               %Density of the gas [kg/m3]
Vp=0;                                   %Far-field loading rate [m/s]
cs=1645;                                %Shear wave speed in rock [m/s]
% cp=2568;                              %Compressional wave speed in rock [km/s]
nu=0.15;                                %Poisson ratio of rock
G=rho*cs*cs;                            %Shear modulus
lambda=2*G*(1 + nu)/3/(1-2*nu)-2/3*G;   %First Lame parameter
eta=G/2/cs;                             %Viscosity
g=9.81;                                 %Gravitation 
K0=0.75;                                %Gradient between sigma_max en simga_min

%slip rate parameters
Vw=1e90;                                %Dynamic weakening velocity
Vi=1e-30;                               %Background slip rate

%initial parameters
mu0=0.6;                                %Friction coefficient 
V0=1e-6;                                %Reference slip rate
a0=0.015;                               %Rate-and-state direct effect
b0=0.01;                               %Rate-and-state evolution effect
L=0.5;                                  %Characteristic slip distance 

save('input.mat','alpha','xsize', 'ysize', 'Nx', 'Ny', 'rho', 'rhof','rhog', 'Vp', 'cs', 'nu', 'G', 'lambda', 'eta', 'g', 'K0', 'Vw', 'Vi', 'mu0', 'V0', 'a0' ,'b0', 'L' );


end

