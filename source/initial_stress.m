function [sigman0,tau0,Pl0,Pr0]=initial_stress(rho,rhof,rhog,g,alpha,K0,Biot,z,param)
TBt=param.BZb-param.ytop;
SSt=param.TBb-param.ytop;
SSb=param.SSb-param.ytop;
offset=param.offset;

sigmatop=rho*g*2000-4.70e6;
Ptop=rhof*g*2000;
sigmav=rho*g*z+sigmatop;
% sigmav=20.0e6*(z/1000+2).^1.08; % Buijze lower 0.1-0.4 MPa
% Pl0=rhof*g*z+(rhof-rhog)*g*(SSb-offset-z).*(z>=SSt-offset).*(z<=SSb-offset)+Ptop+1.16e6*(z>SSt-offset);
% Pr0=rhof*g*z+(rhof-rhog)*g*(SSb-offset-z).*(z>=SSt).*(z<=SSb)+Ptop+1.16e6*(z>SSt);
Pl0=rhof*g*z+(rhof-rhog)*g*(SSb-offset-z).*(z>=TBt-offset).*(z<=SSb-offset)+Ptop+1.16e6*(z>TBt-offset);
Pr0=rhof*g*z+(rhof-rhog)*g*(SSb-offset-z).*(z>=TBt).*(z<=SSb)+Ptop+1.16e6*(z>TBt);
sigman0=(1+K0)/2.*sigmav+(1-K0)/2*cosd(2*alpha).*sigmav-(z<SSb-offset).*Pl0*Biot-(z>=SSb-offset).*Pr0*Biot;
tau0=(1-K0)/2*sind(2*alpha).*sigmav;

end
