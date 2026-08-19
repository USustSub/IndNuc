function [sigma,tau,U,V,theta]=initial_fault(L,V0,Vi,mu0,eta,Ny,a,b,tau0,sigman0,theta0,param)
% BP3-QD initial slip, slip rate, and state (benchmark equations 24-26).

U=zeros(Ny,1);
V=zeros(Ny,1)+param.Vinit;

friction_at_initial=(tau0-eta*V)./sigman0;
argument=2*V0./V.*sinh(friction_at_initial./a);
if any(argument<=0)
    error('BP3 initial-state inversion produced a non-positive logarithm argument.');
end
theta=L./V0.*exp((a.*log(argument)-mu0)./b);

sigma=sigman0;
tau=tau0-eta*V;
end
