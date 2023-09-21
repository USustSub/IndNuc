function [sigma,tau,U,V,theta]=initial_fault(L,V0,Vi,mu0,eta,Ny,a,b,tau0,sigman0)
    U=zeros(Ny,1);
    V=zeros(Ny,1)+Vi; 
    theta=L/V0*exp(a./b.*log(2*V0/Vi*sinh((tau0-eta*Vi)./a./sigman0))-mu0./b);
    sigma=zeros(Ny,1)+sigman0;
    tau=tau0-eta*V;
end

% fixed Vi
%     U=zeros(Ny,1);
%     V=zeros(Ny,1)+Vi; 
%     theta=L/V0*exp(a./b.*log(2*V0/Vi*sinh((tau0-eta*Vi)./a./sigman0))-mu0./b);
%     sigma=zeros(Ny,1)+sigman0;
%     tau=tau0-eta*V;

% steady state
%     U=zeros(Ny,1);
%     theta=zeros(Ny,1)+L/V0*exp(-0);
%     V=V0*sinh(tau0./a./sigman0).*exp(-(mu0+b.*log(theta*V0/L))./a);
%     sigma=zeros(Ny,1)+sigman0;
%     tau=tau0-eta*V;