tic;
% clear all;
% modified from rsd2d_qd_v_pressure, add fluid pressure
% to reservoir setup + gravity
% correct bc with total stress continuity
% add flash heating
% adjust dPdtSS
% updated initial condition; modulized param; add Biot 25.3.6
% change TB pressure 25.3.25
% fixed GP, lambdaP, fixed top BC and pressure at transition 25.5.5
% add dPdtSS diffusion 1D 25.4.3, introduce Ddiff 25.5.6
% parametric version

colordot='r.-';
% generate_parameters();
load('parameters.mat');
% load('randomvars.mat');
param=convert_thickness(param);
save('parameters','param','-append');
checkpointer=param.checkpointer;
yr=365*24*60*60;
g=9.81;

output_interval=param.output_interval;
checkpoint_interval=param.checkpoint_interval;

Nt=param.Nt;

% Defining numerical model
alpha=param.alpha;
sina=sind(alpha);
cosa=cosd(alpha);
xsize=param.xsize; % Horizontal model size, m
ysize=param.ysize; % Vertical model size, m
Nx=param.Nx; % Horizontal model resolution (number of points)
Ny=param.Ny; % Vertical model resolution
N=(Nx+1)*(Ny+1)*2;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m
dz=dy*sina;

if (~checkpointer)
  x=-xsize/2:dx:xsize/2; % Vector of coordinates of basic grid points, m
  y=(0:dy:ysize)'; % Vertical coordinates of basic grid points, m
  xp=(-xsize/2-dx/2):dx:(xsize/2+dx/2);
  yp=((-dy/2):dy:(ysize+dy/2))';
  Xuy=y*cosa+xp;
  Yuy=y*sina+xp*0;
  Xux=yp*cosa+x;
  Yux=yp*sina+x*0;
  Xtau=y*cosa+x;
  Ytau=y*sina+x*0;
  Xsigma=yp(2:Ny,1)*cosa+xp(1,2:Nx);
  Ysigma=yp(2:Ny,1)*sina+xp(1,2:Nx)*0;
  z=y*sina;
  save('coord.mat','x','y','xp','yp','Xtau','Xux','Xuy','Ytau','Yux','Yuy','Xsigma','Ysigma','z');
else
  load('coord.mat');
  % load('parameters.mat');
end

[rho,lambda,G,eta,K0,a,b,L,mu0,V0]=build_layer_rsf(Nx,Ny,x,z,param);
lambdaP=movmean(movmean(lambda,2,2,"omitnan",'Endpoints','discard'),2,1,'Endpoints','discard');
GP=movmean(movmean(G,2,2,"omitnan",'Endpoints','discard'),2,1,'Endpoints','discard');
save('parameters','rho','lambda','G','lambdaP','GP','eta','K0','a','b','L','mu0','V0','-append');

tload=param.tload;
Vw=param.Vw;
muw=param.muw;
Vi=param.Vi;
TBt=param.BZb-param.ytop;
SSt=param.TBb-param.ytop;
SSb=param.SSb-param.ytop;
offset=param.offset;
Biot=param.Biot;

[sigman0,tau0,Pl0,Pr0]=initial_stress(param.rhoss,param.rhof,param.rhog,g,param.alpha,K0,Biot,z,param);

Pl=Pl0;
Pr=Pr0;
P=(z<SSb-offset).*Pl+(z>=SSb-offset).*Pr;

% dPdtSS=-0e-2;
% dPdtTB=0;
% dPdt=build_dPdt(Ny,z,param);
% dPdt.L=zeros(Ny,1);
% dPdt.R=zeros(Ny,1);

ux=zeros(Ny+1,Nx);
uy=zeros(Ny,Nx+1);
vx=zeros(Ny+1,Nx);
vy=zeros(Ny,Nx+1);
% vpx=zeros(Ny+1,Nx+1);
% vpy=zeros(Ny+1,Nx+1);
tauqs=zeros(Ny,Nx);
sigmaqs=zeros(Ny-1,Nx-1);
% U=zeros(Ny,1);
% V=zeros(Ny,1);
t=0;
t2=-1;
dt=1e+0;
dt_max=1e15;

dPdt=build_dPdt(Ny,z,t2,param);

if (~checkpointer)
    %initial condition
    [sigma,tau,U,V,theta]=initial_fault(L,V0,Vi,mu0,eta,Ny,a,b,tau0,sigman0,param.theta0);
    %numerical solver
    LH=build_LH(lambda,G,sina,cosa,Nx,Ny,N,dx,dy);
    RH=build_RH(lambda,G,sina,cosa,dPdt,Biot,Nx,Ny,N,dx,dy,y,V,z,dz,param);
    ksi=build_ksi(G(:,(Nx+1)/2),L,dy,a,b,sigman0);
    save('initiation.mat','LH','RH','ksi','-v7.3');
    %output
    fid=fopen('output.txt','w+');
    flag=0;
else
    %initial condition
    load(['data_', int2str(checkpointer), '.mat']);
    %numerical solver
    load('initiation.mat');
    if(flag==2)%
%         dPdt=build_dPdt(Ny,z,t2,param);
        % RH=build_RH(lambda,G,sina,cosa,dPdt,Biot,Nx,Ny,N,dx,dy,y,V,z,dz,param);
        dt_max=1e6;
    end
    %output
    fid=fopen('output.txt','a+');
end



global Um Vm taum sigmam thetam Pm dtm tm tm2 taumall sigmamall uymall vymall uxmall vxmall
Um=zeros(Ny,Nt/output_interval);
Vm=zeros(Ny,Nt/output_interval);
taum=zeros(Ny,Nt/output_interval);
sigmam=zeros(Ny,Nt/output_interval);
Pm=zeros(Ny,Nt/output_interval);
thetam=zeros(Ny,Nt/output_interval);
dtm=zeros(1,Nt/output_interval);
tm=zeros(1,Nt/output_interval);
tm2=zeros(1,Nt/output_interval);
% taumall=zeros(Ny,Nx,Nt/output_interval);
% sigmamall=zeros(Ny-1,Nx-1,Nt/output_interval);
% uymall=zeros(Ny,Nx+1,Nt/output_interval);
% vymall=zeros(Ny,Nx+1,Nt/output_interval);
% uxmall=zeros(Ny+1,Nx,Nt/output_interval);
% vxmall=zeros(Ny+1,Nx,Nt/output_interval);


dLH=decomposition(LH);
toc;

options = optimset('TolFun', 5,'TolX', 0);

for it=1:Nt
    if (flag==1)
%         dPdt=build_dPdt(Ny,z,t2,param);
%         RH=build_RH(lambda,G,sina,cosa,dPdt,Biot,Nx,Ny,N,dx,dy,y,V,z,dz,param);
%         save('initiation.mat','RH','-append');
        dt=1e0;
        dt_max=1e6;
        t2=0;
        flag=2;
    end

    dPdt=build_dPdt(Ny,z,t2,param);
    RH=build_RH(lambda,G,sina,cosa,dPdt,Biot,Nx,Ny,N,dx,dy,y,V,z,dz,param);
%     plot(z,dPdt.L)
%     hold on
%     plot(z,dPdt.R)
% %     hold off
%     drawnow limitrate;

    if (it<1)
        for iy=1:Ny
%             V(iy)=fzero(@(VV) sigma(iy)*a(iy)*asinh(VV/(2*V0)*exp((mu0+b(iy)*log(V0*theta(iy)/L))/a(iy))) + eta*VV - (tauqs(iy,(Nx+1)/2)+tau0(iy)), V(iy));
            V(iy)=fzero(@(VV) sigma(iy)*a(iy)*asinh(VV/(2*V0)*exp((mu0(iy)+b(iy)*log(V0*theta(iy)/L(iy)))/a(iy))) / (1+L(iy)/Vw/theta(iy)) + eta*VV - (tauqs(iy,(Nx+1)/2)+tau0(iy)), V(iy));
        end
    else
%         [V,fval,exitflag]=bisection(@(VV) sigma.*a.*asinh(VV/(2*V0).*exp((mu0+b.*log(V0*theta/L))./a)) + eta*VV - (tauqs(:,(Nx+1)/2)+tau0),0,V*2,zeros(Ny,1),options);
        [V,fval,exitflag]=bisection(@(VV) sigma.*((a.*asinh(VV/(2*V0).*exp((mu0+b.*log(V0*theta./L))./a))-muw)./(1+L/Vw./theta)+muw) + eta*VV - (tauqs(:,(Nx+1)/2)+tau0),0,max(V)*2+zeros(Ny,1),zeros(Ny,1),options);
    end

%     semilogy(V);drawnow;

    V=max(V,1e-40);
%     V(1)=V(2); V(end)=V(end-1);
    %dt=max(0.1*L/V,1e-150);
%     dt=min(max(min(0.13*L./V),1e-150),1.2*dt);
    dt=min(max(min(ksi(2:Ny-1).*L(2:Ny-1)./V(2:Ny-1)),1e-150),1.2*dt);
    dt=min(dt_max,dt);
    if (flag==0&&t+dt>=tload)
        dt=tload-t;
        flag=1;
    end

%     theta=theta+dt*(1-V.*theta./L);
    expo=V*dt./L>1e-6;
    theta=expo.*(L./V.*(1-exp(-V*dt./L))+theta.*exp(-V*dt./L))+(~expo).*(theta+dt*(1-V.*theta./L));
    tau=tauqs(:,(Nx+1)/2)+tau0-eta*V;
    U=U+dt*V;


    %     RH((1:Ny)*2)=-Vp;
    %     RH((Nx*(Ny+1)+(1:Ny))*2)=Vp;
    RH(((Nx-1)/2*(Ny+1)+(2:Ny-1))*2)=V(2:Ny-1);

    S=dLH\RH;

    vpx=reshape(S(1:2:end),Ny+1,Nx+1);
    vpy=reshape(S(2:2:end),Ny+1,Nx+1);

    vy=vpy(1:Ny,:);
    vx=vpx(:,1:Nx);
    uy=uy+vy*dt;
    ux=ux+vx*dt;
    [~,duydy]=gradient(uy,dx,dy);
    [duxdx,~]=gradient(ux,dx,dy);

    tauqs=G/sina.*(diff(uy,1,2)/dx+(1-2*cosa*cosa)*diff(ux,1,1)/dy+cosa*(movmean(duxdx,2,1,'Endpoints','discard')-movmean(duydy,2,2,'Endpoints','discard')));
    tauqs(:,(Nx+1)/2)=(tauqs(:,(Nx+1)/2-1)+tauqs(:,(Nx+1)/2+1))/2;
    %     tauqs=G/sina*(diff(uy,1,2)/dx+cosa*movmean(duxdx,2,1,'Endpoints','discard'));
    %     tauqs(:,(Nx+1)/2)=(tauqs(:,(Nx+1)/2-1)+tauqs(:,(Nx+1)/2+1))/2;
    %     tauqs=tauqs+G/sina*((1-2*cosa*cosa)*diff(ux,1,1)/dy+cosa*(-movmean(duydy,2,2,'Endpoints','discard')));
    sigmaqs=(lambdaP+2*GP).*diff(ux(2:Ny,:),1,2)/dx+lambdaP.*diff(uy(:,2:Nx),1,1)/dy-2*GP*cosa.*movmean(movmean(diff(ux,1,1)/dy,2,2,'Endpoints','discard'),2,1,'Endpoints','discard');
%     sigmatemp=(sigmaqs(:,(Nx-1)/2)+sigmaqs(:,(Nx+1)/2))/2;
%     sigmatemp=max(sigmaqs(:,(Nx-1)/2)+Pl0,sigmaqs(:,(Nx+1)/2)+Pr0);
%     sigma=sigman0+max(Pl0,Pr0)-[sigmatemp(1);movmean(sigmatemp,2,'Endpoints','discard');sigmatemp(end)];
    sigmal=[sigmaqs(1,(Nx-1)/2);movmean(sigmaqs(:,(Nx-1)/2),2,'Endpoints','discard');sigmaqs(end,(Nx-1)/2)];
    sigmar=[sigmaqs(1,(Nx+1)/2);movmean(sigmaqs(:,(Nx+1)/2),2,'Endpoints','discard');sigmaqs(end,(Nx+1)/2)];
    sigmatemp=min(sigmal,sigmar);
%     sigmatemp=(z<SSb-offset).*sigmal+(z>=SSb-offset).*sigmar; %avoid updip error
    sigma=sigman0-sigmatemp;

    Pl=dPdt.L*dt+Pl;
    Pr=dPdt.R*dt+Pr;
    P=(z<SSb-offset).*Pl+(z>=SSb-offset).*Pr;


    fprintf(fid,'it=%d, t=%f, t2=%f, dt=%e, maxV=%e, minV=%e, maxU=%f\n',checkpointer+it,t/yr,t2/yr,dt,max(V),min(V),max(U));

    if(mod(it,output_interval)==0)
        write_memory(it,output_interval,U,V,tau,sigma,P,theta,dt,t,t2,tauqs,sigmaqs,uy,vy,ux,vx);
    end

    t=t+dt;
    if (flag==2)
        t2=t2+dt;
    end

%     if(mod(it,output_interval)==0)
%         plot(z,tau/1e6)
%         hold on
%         plot(z,sigma/1e6)
%         plot(z,P/1e6)
%         plot(z,log10(V))
%         plot(z,tau./sigma*10)
%         hold off
%         drawnow limitrate;
%     end


    if(mod(it,checkpoint_interval)==0)
        save(['data_', int2str(checkpointer+it), '.mat'],'U','V','tau','sigma','P','theta','dt','t','t2','tauqs','sigmaqs','uy','vy','ux','vx','flag');
        save('dataall.mat','Um','Vm','taum','sigmam','Pm','thetam','dtm','tm','tm2','-v7.3');
        disp(it);
        toc;
    end
end

save(['data_',datestr(now,'yyyy-mm-dd-HH-MM-ss'),'.mat'],'Um','Vm','taum','sigmam','Pm','thetam','dtm','tm','tm2','-v7.3');
fclose(fid);
% toc;

% music=load('handel.mat');
% sound(music.y,music.Fs);
