function [sigman0,tau0,Pl0,Pr0]=initial_stress(rho,rhof,rhog,g,alpha,K0,y)

sigmatop=rho*g*2000-4.70e6;
Ptop=rhof*g*2000;
sigmav=rho*g*y+sigmatop;
Pl0=rhof*g*y+(rhof-rhog)*g*(1000-y).*(y>800).*(y<=1000)+Ptop+1.16e6*(y>800);
Pr0=rhof*g*y+(rhof-rhog)*g*(1000-y).*(y>850).*(y<=1050)+Ptop+1.16e6*(y>850);
sigman0=(1+K0)/2*sigmav+(1-K0)/2*cosd(2*alpha)*sigmav-(y<1000).*Pl0-(y>=1000).*Pr0;
tau0=(1-K0)/2*sind(2*alpha)*sigmav;

end
