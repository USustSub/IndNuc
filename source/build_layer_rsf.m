function [rho,lambda,G,eta,K0,a,b,L,mu0,V0]=build_layer_rsf(Nx,Ny,x,z,param)
ytop = param.ytop;
offset = param.offset;
%depth ranges of the layers (left side of the fault)
RSt = param.ytop-ytop;             %Rocksalt top
RSb = param.RSb-ytop;             %Rocksalt bottom
BZt = RSb;             %Basal zechstein top
BZb = param.BZb-ytop;             %Basal zechstein bottom
TBt = BZb;             %Ten boer top
TBb = param.TBb-ytop;             %Ten Boer bottom
SSt = TBb;             %Slochteren sandstone top
SSb = param.SSb-ytop;             %Slochteren sandstone bottom
Ct = SSb;              %Carbonefious top
Cb = param.Cb-ytop;              %Carbonefious bottom

%2015 model

    density=[2.1 2.9 2.3 param.rhoss/1000 2.3]*1000;
    Young=[30 70 40 15 40]*1e9;
    Poisson=[0.35 0.25 0.2 param.nuss 0.2];
    lame_2=Young/2./(1+Poisson);
    lame_1=lame_2*2.*Poisson./(1-2*Poisson);

    rocktype=zeros(Ny,Nx);
    rocktype(z<=RSb,x>0)=1;
    rocktype(z>BZt&z<=BZb,x>0)=2;
    rocktype(z>TBt&z<=TBb,x>0)=3;
    rocktype(z>SSt&z<=SSb,x>0)=4;%
    rocktype(z>Ct,x>0)=5;
    rocktype(z<=RSb-offset,x<=0)=1;
    rocktype(z>BZt-offset&z<=BZb-offset,x<=0)=2;
    rocktype(z>TBt-offset&z<=TBb-offset,x<=0)=3;
    rocktype(z>SSt-offset&z<=SSb-offset,x<=0)=4;
    rocktype(z>Ct-offset,x<=0)=5;
%     rocktype(:,:)=4;

    rho=density(rocktype);
    lambda=lame_1(rocktype);
    G=lame_2(rocktype);

    rho(:,(Nx+1)/2)=nan;
    lambda(:,(Nx+1)/2)=nan;
    G(:,(Nx+1)/2)=nan;

    eta=sqrt(density(4)*lame_2(4))/2;

%hetrogenious fricitonal properties code
AA=[0.002 0.018 0.002 param.ass 0.002];
BB=[0.0001 0.02 0.0001 param.bss 0.0001];
% AA=[0.002 0.008 0.002 0.006 0.002];
% BB=[0.0001 0.01 0.0001 0.004 0.0001];
LL=[0.5e-3 0.5e-3 0.5e-3 param.Lss 0.5e-3];
mu0mu0=[0.5 0.6 0.4 param.mu0ss 0.5];

% rocktype(:,:)=4;
%     rocktype=zeros(Ny,Nx);
%     rocktype(z<=RSb,x>0)=1;
%     rocktype(z>BZt&z<=BZb,x>0)=2;
%     rocktype(z>TBt&z<=TBb,x>0)=3;
%     rocktype(z>SSt&z<=SSb,x>0)=4;
%     rocktype(z>Ct,x>0)=5;
%     rocktype(z<=RSb-offset,x<=0)=1;
%     rocktype(z>BZt-offset&z<=BZb-offset,x<=0)=2;
%     rocktype(z>TBt-offset&z<=TBb-offset,x<=0)=3;
%     rocktype(z>SSt-offset&z<=SSb-offset,x<=0)=4;
%     rocktype(z>Ct-offset,x<=0)=5;

rockl=rocktype(:,(Nx+1)/2-1);
rockr=rocktype(:,(Nx+1)/2+1);

% a = zeros(Ny,1)+0.01;
al = AA(rockl);
ar = AA(rockr);
a = min(ar,al)';
% b = zeros(Ny,1)+0.008;
bl = BB(rockl);
br = BB(rockr);
b = min(br,bl)';
L=LL(rockr)'; % 25.3.11
mu0l = mu0mu0(rockl);
mu0r = mu0mu0(rockr);
mu0 = max(mu0l,mu0r)';
% mu0 = (mu0l'+mu0r')/2;
% mu0 = zeros(Ny,1)+0.6;

V0=param.V0;

K0K0=[1 0.9 param.K0ss param.K0ss param.K0ss];
K0l = K0K0(rockl);
K0r = K0K0(rockr);
K0 = (K0l'+K0r')/2;

% geometry = struct('TBt',TBt,'SSt',SSt,'SSb',SSb,'offset',offset);

end

% Use Groningen Velocity Model 2017 - Groningen full elastic velocity model
% Autor(s) Remco Romijn
% rock type,rho,Vp,Vs,nu,lambda,G
% rocksalt 4400 2486 2.09
% zechstein 5900 3238 2.81
% ten boer
% sandstone 3900 2286 2.46
% carboniferous Vp=0.541*Z + 2572.3 Vs = 0.927 * Vp - 1547.313 2.65
% at depth 3200 m, 4300 2442 2.65
% Poisson 0.2655 0.2845 0.2383 0.2620

%     density=[2.09 2.81 2.46 2.46 2.65]*1000;
%     speed_P=[4400 5900 3900 3900 4300];
%     speed_S=[2486 3238 2286 2286 2442];
%     lame_1=density.*(speed_P.*speed_P-2*speed_S.*speed_S);
% %     lame_2=density.*speed_S.*speed_S;
% %     Poisson_ratio=lame_1./(lame_1+lame_2)/2;
%     Poisson_ratio=[0.15 0.15 0.15 0.15 0.15];
%     lame_2=lame_1.*(1/2./Poisson_ratio-1);
%
%     rocktype=zeros(Ny,Nx);
%     rocktype(z<=RSb,x>0)=1;
%     rocktype(z>BZt&z<=BZb,x>0)=2;
%     rocktype(z>TBt&z<=TBb,x>0)=3;
%     rocktype(z>SSt&z<=SSb,x>0)=4;
%     rocktype(z>Ct,x>0)=5;
%     rocktype(z<=RSb-offset,x<=0)=1;
%     rocktype(z>BZt-offset&z<=BZb-offset,x<=0)=2;
%     rocktype(z>TBt-offset&z<=TBb-offset,x<=0)=3;
%     rocktype(z>SSt-offset&z<=SSb-offset,x<=0)=4;
%     rocktype(z>Ct-offset,x<=0)=5;
%     rocktype(:,:)=4;
%
%     rho=density(rocktype);
%     Vp=speed_P(rocktype);
%     Vs=speed_S(rocktype);
%     lambda=lame_1(rocktype);
%     G=lame_2(rocktype);
%
%     eta=density(4)*speed_S(4)/2;


%first the a-parameter over depth
% a = zeros(Ny,1)+0.015;
% a(z <= RSb) = 0.00447;                  %Zechstein rocksalt (halite)
% a(z > BZt & z <= BZb) = 0.06895;        %Basal zechstein
% a(z > TBt & z <= TBb) = 0.00305;        %Ten Boer
% a(z > SSt & z <= SSb) = 0.04065;        %Slochteren Sandstone
% a(z > Ct) = 0.02538;                    %Carboniferous member

%and for b
% b = zeros(Ny,1); bl = zeros(Ny,1); br = zeros(Ny,1);
% br(z <= RSb) = 0.001;                  %Zechstein rocksalt (halite)
% br(z > BZt & z <= BZb) = 0.02;        %Basal zechstein
% br(z > TBt & z <= TBb) = 0.001;        %Ten Boer
% br(z > SSt & z <= SSb) = 0.02;        %Slochteren Sandstone
% br(z > Ct) = 0.001;                    %Carboniferous member
% bl(z <= RSb-offset) = 0.001;                  %Zechstein rocksalt (halite)
% bl(z > BZt-offset & z <= BZb-offset) = 0.02;        %Basal zechstein
% bl(z > TBt-offset & z <= TBb-offset) = 0.001;        %Ten Boer
% bl(z > SSt-offset & z <= SSb-offset) = 0.02;        %Slochteren Sandstone
% bl(z > Ct-offset) = 0.001;                    %Carboniferous member
% b=min(br,bl);

%and for L
% L=zeros(Ny,1)+1e-3;
% L(z <= RSb) = 0.001;                  %Zechstein rocksalt (halite)
% L(z > BZt & z <= BZb) = 0.006;        %Basal zechstein
% L(z > TBt & z <= TBb) = 0.001;        %Ten Boer
% L(z > SSt & z <= SSb) = 0.006;        %Slochteren Sandstone
% L(z > Ct) = 0.001;                    %Carboniferous member

%and for mu0
% mu0=zeros(Ny,1); mu0l=zeros(Ny,1); mu0r=zeros(Ny,1);
% mu0r(z <= RSb) = 0.5;                  %Zechstein rocksalt (halite)
% mu0r(z > BZt & z <= BZb) = 0.6;        %Basal zechstein
% mu0r(z > TBt & z <= TBb) = 0.5;        %Ten Boer
% mu0r(z > SSt & z <= SSb) = 0.6;        %Slochteren Sandstone
% mu0r(z > Ct) = 0.5;                    %Carboniferous member
% mu0l(z <= RSb-offset) = 0.5;                  %Zechstein rocksalt (halite)
% mu0l(z > BZt-offset & z <= BZb-offset) = 0.6;        %Basal zechstein
% mu0l(z > TBt-offset & z <= TBb-offset) = 0.5;        %Ten Boer
% mu0l(z > SSt-offset & z <= SSb-offset) = 0.6;        %Slochteren Sandstone
% mu0l(z > Ct-offset) = 0.5;                    %Carboniferous member
% mu0=max(mu0l,mu0r);
