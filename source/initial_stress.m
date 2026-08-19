function [sigman0,tau0,Pl0,Pr0]=initial_stress(rho,rhof,rhog,g,alpha,K0,Biot,z,param)
% BP3-QD uniform effective normal stress and steady initial shear stress.

Ny=numel(z);
sigman0=zeros(Ny,1)+param.sigma0;

steady_friction=param.amax*asinh(param.Vinit/(2*param.V0)...
    *exp((param.f0+param.b0*log(param.V0/abs(param.Vinit)))... 
    /param.amax));
tau0=zeros(Ny,1)+param.sigma0*steady_friction+...
    sqrt(param.rho*(param.rho*param.cs^2))/2*param.Vinit;

Pl0=zeros(Ny,1);
Pr0=zeros(Ny,1);
end
