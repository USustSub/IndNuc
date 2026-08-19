function [rho,lambda,G,eta,K0,a,b,L,mu0,V0]=build_layer_rsf(Nx,Ny,x,z,param)
% Homogeneous BP3-QD material and down-dip friction profile.

G0=param.rho*param.cs^2;
lambda0=2*G0*param.nu/(1-2*param.nu);

rho=zeros(Ny,Nx)+param.rho;
lambda=zeros(Ny,Nx)+lambda0;
G=zeros(Ny,Nx)+G0;

% The central column is the internal fault interface, not bulk material.
fault_ix=(Nx+1)/2;
rho(:,fault_ix)=nan;
lambda(:,fault_ix)=nan;
G(:,fault_ix)=nan;

eta=sqrt(param.rho*G0)/2;

xd=z/sind(param.alpha);
a=zeros(Ny,1)+param.amax;
a(xd<param.H)=param.a0;
transition=xd>=param.H & xd<param.H+param.h;
a(transition)=param.a0+(param.amax-param.a0)...
    .*(xd(transition)-param.H)/param.h;

b=zeros(Ny,1)+param.b0;
L=zeros(Ny,1)+param.L0;
mu0=zeros(Ny,1)+param.f0;
V0=param.V0;
K0=zeros(Ny,1);
end
