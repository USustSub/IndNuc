tic;   %start the counting of running time of the model.
%% Seismo-thermomechanical numerical model with;
% fluid pressure
% to reservoir setup + gravity
% correct bc with total stress continuity
% flash heating
% adjustable dPdt

clear all;

addpath('.\source\')
addpath('.\visualization\')

colordot='r.';              %Red lines with . = only dots or .- = line+dot
generate_input();           %Call the input file
load('input.mat');          %Load in the parameters in this matlab code

checkpointer=0;             %If you don't want to start at 0, use checkpointer
yr=365*24*60*60;            %Not standard in matlab.

output_interval=10;         %The interval saved in the output.txt file
checkpoint_interval=1000;   %The checkpoint interval, for all the details saved

Nt=1000;                    %timestep.

%% Defining the numerical model

sina=sind(alpha);
cosa=cosd(alpha);
N=(Nx+1)*(Ny+1)*2;
dx=xsize/(Nx-1);            % Horizontal grid step, m
dy=ysize/(Ny-1);            % Vertical grid step, m
x=-xsize/2:dx:xsize/2;      % Vector of coordinates of basic grid points, m
y=(0:dy:ysize)';            % Vertical coordinates of basic grid points, m
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
save('coord.mat','x','y','xp','yp','Xtau','Xux','Xuy','Ytau','Yux','Yuy','Xsigma','Ysigma');

%arrays to store the values of displacement in x-direction, displacement in y-direction,
%vertical velocity, horizontal velocity ,shear stress, normal stress
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

%% Time and time interval

t=0;
t2=0;
dt=1e+0;

%% Frictional parameters a and b
% a and b, the parameters that discribe the relative influence of direct and
% evolutiary effects.

%homogenious fricitonal properties
%a = zeros(Ny,1)+ a0 ;
b=zeros(Ny,1)+b0;

%depth ranges of the layers (left side of the fault)
RSt = 2000-ysize;             %Rocksalt top
RSb = 2730-ysize;             %Rocksalt bottom
BZt = 2730-ysize;             %Basal zechstein top
BZb = 2780-ysize;             %Basal zechstein bottom
TBt = 2780-ysize;             %Ten boer top
TBb = 2850-ysize;             %Ten Boer bottom
SSt = 2850-ysize;             %Slochteren sandstone top
SSb = 3050-ysize;             %Slochteren sandstone bottom
Ct = 3050-ysize;              %Carbonefious top
Cb = 4000-ysize;              %Carbonefious bottom

%hetrogenious fricitonal properties code
%first the a-parameter over depth
a = zeros(Ny,1);
a(y <= RSb) = 0.00447;                  %Zechstein rocksalt (halite)
a(y > BZt & y <= BZb) = 0.06895;        %Basal zechstein
a(y > TBt & y <= TBb) = 0.00305;        %Ten Boer
a(y > SSt & y <= SSb) = 0.04065;        %Slochteren Sandstone
a(y > Ct) = 0.02538;                    %Carboniferous member

%and for b
b = zeros(Ny,1);

b(y <= RSb) = -0.0059;                  %Zechstein rocksalt (halite)
b(y > BZt & y <= BZb) = 0.07209;        %Basal zechstein
b(y > TBt & y <= TBb) = -0.00093;       %Ten Boer
b(y > SSt & y <= SSb) = 0.03796;        %Slochteren Sandstone
b(y > Ct) = 0.02347;                    %Carboniferous member


%% Stress input

[sigman0,tau0,Pl0,Pr0]=initial_stress(rho,rhof,rhog,g,alpha,K0,y);

Pl=Pl0;  %initial pressure for l
Pr=Pr0;  %initial pressure for r
P=(y<1000).*Pl+(y>=1000).*Pr; %pressure calculation

dPdt = -0e-2;    %pressure rate before depletion

%% Start building the RHS and LHS to solve the PDE's

if (~checkpointer)
    %initial condition
    [sigma,tau,U,V,theta]=initial_fault(L,V0,Vi,mu0,eta,Ny,a,b,tau0,sigman0);
    %numerical solver
    LH=build_LH(lambda,G,sina,cosa,Nx,Ny,N,dx,dy);
    RH=build_RH(lambda,G,sina,cosa,dPdt,Nx,Ny,N,dx,dy,y,V);
    ksi=build_ksi(G,L,dy,a,b,sigman0);
    save('initiation.mat','LH','RH','ksi','-v7.3');
    %output
    fid=fopen('output.txt','w+');
else
    %initial condition
    load(['data_', int2str(checkpointer), '.mat']);
    %numerical solver
    load('initiation.mat');
    %output
    fid=fopen('output.txt','a+');
end

global Um Vm taum sigmam thetam Pm dtm tm taumall sigmamall uymall vymall uxmall vxmall
Um=zeros(Ny,Nt/output_interval);
Vm=zeros(Ny,Nt/output_interval);
taum=zeros(Ny,Nt/output_interval);
sigmam=zeros(Ny,Nt/output_interval);
Pm=zeros(Ny,Nt/output_interval);
thetam=zeros(Ny,Nt/output_interval);
dtm=zeros(1,Nt/output_interval);
tm=zeros(1,Nt/output_interval);
taumall=zeros(Ny,Nx,Nt/output_interval);
sigmamall=zeros(Ny-1,Nx-1,Nt/output_interval);
uymall=zeros(Ny,Nx+1,Nt/output_interval);
vymall=zeros(Ny,Nx+1,Nt/output_interval);
uxmall=zeros(Ny+1,Nx,Nt/output_interval);
vxmall=zeros(Ny+1,Nx,Nt/output_interval);


dLH=decomposition(LH);
toc;
dt_max=1e6;
time=0;
tload=10*yr;    %time you want to have heating (in s)

%% Calculations

options = optimset('TolFun', 5,'TolX', 0);
for it=1:Nt
    if (time==1)
        dPdt=-0.0127;
        RH=build_RH(lambda,G,sina,cosa,dPdt,Nx,Ny,N,dx,dy,y,V);
        dt=1e0;
        dt_max=1e6;
        t2=0;
        time=2;
    end
    % rate-and-state solver
        for iy=1:Ny
%             V(iy)=fzero(@(VV) sigma(iy)*a(iy)*asinh(VV/(2*V0)*exp((mu0+b(iy)*log(V0*theta(iy)/L))/a(iy))) + eta*VV - (tauqs(iy,(Nx+1)/2)+tau0(iy)), V(iy));
            V(iy)=fzero(@(VV) sigma(iy)*a(iy)*asinh(VV/(2*V0)*exp((mu0+b(iy)*log(V0*theta(iy)/L))/a(iy))) / (1+L/Vw/theta(iy)) + eta*VV - (tauqs(iy,(Nx+1)/2)+tau0(iy)), V(iy));
        end

    V=max(V,1e-40);
    V(1)=V(2); V(end)=V(end-1);
    %dt=max(0.1*L/V,1e-150);

    % adaptive time stepping
%     dt=min(max(min(0.13*L./V),1e-150),1.2*dt);
    dt=min(max(min(ksi(2:Ny-1)*L./V(2:Ny-1)),1e-150),1.2*dt);
    dt=min(dt_max,dt);
    %t=t+dt;

    if (time==0&&t+dt>=tload)
        dt=tload-t;
        time=1;
    end

    % aging law
    theta=theta+dt*(1-V.*theta/L);
    tau=tauqs(:,(Nx+1)/2)+tau0-eta*V;
    U=U+dt*V;


    %     RH((1:Ny)*2)=-Vp;
    %     RH((Nx*(Ny+1)+(1:Ny))*2)=Vp;
    RH(((Nx-1)/2*(Ny+1)+(2:Ny-1))*2)=V(2:Ny-1);

%     RH(((1:((Nx-1)/2-1))*(Ny+1)+40)*2)=dPdt/dy*dx*dx/G;
%     RH((((Nx+1)/2:(Nx-1))*(Ny+1)+40)*2)=dPdt/dy*dx*dx/G;
%     RH(((1:((Nx-1)/2-1))*(Ny+1)+50)*2)=-dPdt/dy*dx*dx/G;
%     RH((((Nx+1)/2:(Nx-1))*(Ny+1)+50)*2)=-dPdt/dy*dx*dx/G;

    %solver
    S=dLH\RH;

    vpx=reshape(S(1:2:end),Ny+1,Nx+1);
    vpy=reshape(S(2:2:end),Ny+1,Nx+1);

    vy=vpy(1:Ny,:);
    vx=vpx(:,1:Nx);

    %calculate slip
    uy=uy+vy*dt;
    ux=ux+vx*dt;
    [~,duydy]=gradient(uy,dx,dy);
    [duxdx,~]=gradient(ux,dx,dy);

    %calculation of stress (definition of strain + Hooke's law)
    tauqs=G/sina*(diff(uy,1,2)/dx+(1-2*cosa*cosa)*diff(ux,1,1)/dy+cosa*(movmean(duxdx,2,1,'Endpoints','discard')-movmean(duydy,2,2,'Endpoints','discard')));
    tauqs(:,(Nx+1)/2)=(tauqs(:,(Nx+1)/2-1)+tauqs(:,(Nx+1)/2+1))/2;
    %     tauqs=G/sina*(diff(uy,1,2)/dx+cosa*movmean(duxdx,2,1,'Endpoints','discard'));
    %     tauqs(:,(Nx+1)/2)=(tauqs(:,(Nx+1)/2-1)+tauqs(:,(Nx+1)/2+1))/2;
    %     tauqs=tauqs+G/sina*((1-2*cosa*cosa)*diff(ux,1,1)/dy+cosa*(-movmean(duydy,2,2,'Endpoints','discard')));
    sigmaqs=(lambda+2*G)*diff(ux(2:Ny,:),1,2)/dx+lambda*diff(uy(:,2:Nx),1,1)/dy-2*G*cosa*movmean(movmean(diff(ux,1,1)/dy,2,2,'Endpoints','discard'),2,1,'Endpoints','discard');
%     sigmatemp=(sigmaqs(:,(Nx-1)/2)+sigmaqs(:,(Nx+1)/2))/2;
%     sigmatemp=max(sigmaqs(:,(Nx-1)/2)+Pl0,sigmaqs(:,(Nx+1)/2)+Pr0);
%     sigma=sigman0+max(Pl0,Pr0)-[sigmatemp(1);movmean(sigmatemp,2,'Endpoints','discard');sigmatemp(end)];
    sigmal=[sigmaqs(1,(Nx-1)/2);movmean(sigmaqs(:,(Nx-1)/2),2,'Endpoints','discard');sigmaqs(end,(Nx-1)/2)];
    sigmar=[sigmaqs(1,(Nx+1)/2);movmean(sigmaqs(:,(Nx+1)/2),2,'Endpoints','discard');sigmaqs(end,(Nx+1)/2)];
%     sigmatemp=(y<1000).*sigmal+(y>=1000).*sigmar;
    sigmatemp=min(sigmal,sigmar);
    sigma=sigman0-sigmatemp;

    %calculate pressure
    Pl=dPdt*dt*(y>800).*(y<=1000)+Pl;
    Pr=dPdt*dt*(y>850).*(y<=1050)+Pr;
    P=(y<1000).*Pl+(y>=1000).*Pr;

    %output.txt
    fprintf(fid,'it=%d, t=%f, dt=%e, maxV=%e, minV=%e, maxU=%f\n',checkpointer+it,t2/yr,dt,max(V),min(V),max(U));

    if(mod(it,output_interval)==0)
        write_memory(it,output_interval,U,V,tau,sigma,P,theta,dt,t,tauqs,sigmaqs,uy,vy,ux,vx);
    end

    t = t+dt;
    if (time==2)
        t2 = t2+dt;
    end




    %to plot the proces, but the numbers are weird, to be able to plot this
    %into one figure! So the numbers don't add up!!
%     if(mod(it,output_interval)==0)
%         plot(tau/1e6)
%         hold on
%         plot(sigma/1e6)
%         plot(log10(V))
%         hold off
%         drawnow limitrate;
%     end


    if(mod(it,checkpoint_interval)==0)
        save(['data_', int2str(checkpointer+it), '.mat'],'U','V','tau','sigma','P','theta','dt','t','tauqs','sigmaqs','uy','vy','ux','vx');
        save('dataall.mat','Um','Vm','taum','sigmam','Pm','thetam','dtm','tm','-v7.3');
        disp(it);
        toc;
        plot_memory(yr,colordot);  %to plot the state param etc.
        drawnow                    %plot it now!
    end
end

% save(['data_',datestr(now,'yyyy-mm-dd-HH-MM-ss'),'.mat'],'Um','Vm','taum','sigmam','Pm','thetam','dtm','tm','-v7.3');
fclose(fid);
% toc;

%% Make or save a plot
saveas(gcf, 'Time_plot_figure_V8_M2.png')

plot_memory(yr,colordot);
saveas(gcf, 'Time_plot_V8_M2.fig')

plot_slip(Vm, Nt, output_interval, Ny, y);
saveas(gcf, 'Slip_plot_V8_M2.fig')

plot_tausimga(sigmam, taum, Nt, output_interval, Ny, y);
saveas(gcf, 'Ratio_tau_sigma_plot_V8_M2.fig')
