function []=generate_parameters()
yr=365*24*60*60;
g=9.81;

param.checkpointer=000;

param.output_interval=10;
param.checkpoint_interval=1000;

param.Nt=50000;

param.xsize=2000;
param.ysize=2000;
param.Nx=1001;
param.Ny=1001;

param.alpha=70;
param.dPdtSS=-0.0127;
param.dPdtTB=param.dPdtSS*0.;
param.hdiff=0;
param.Ddiff=0e-5;

param.rhoss=2400;
param.nuss=0.15;
param.K0ss=0.78;
param.rhof=1150;
param.rhog=200;
param.Biot=1.0;

param.ytop = 2000;
param.offset = 50;
%depth ranges of the layers (left side of the fault)
% param.RSt = 2000;             %Rocksalt top
param.RSb = 2730;             %Rocksalt bottom
% param.BZb = 2780;             %Basal zechstein bottom
% param.TBb = 2850;             %Ten Boer bottom
% param.SSb = 3050;             %Slochteren sandstone bottom
param.Cb = 4000;              %Carbonefious bottom

param.ass=0.011;
param.bss=0.01;
param.Lss=0.5e-3;
param.mu0ss=0.6;
param.V0=1e-6;
param.theta0=100e6*yr;
param.Vw=1e-1;
param.muw=0.3;

param.tload=0e0*yr;

param.Vi=1e-25;


save('parameters.mat','param');

end

% VS
% mu0=0.6;
% V0=1e-6;
% a0=0.005;
% b0=0.0045;
% a_max=0.005;
% L=2e-3;

% VW
% mu0=0.3;
% V0=1e-6;
% a0=0.005;
% b0=0.020;
% a_max=0.025;
% L=5e-3;

% VN
% mu0=0.3;
% V0=1e-6;
% a0=0.0075;
% b0=0.0075;
% a_max=0.0075;
% L=1e-3;
