tic;

project_root=fileparts(mfilename('fullpath'));
addpath(project_root,fullfile(project_root,'source'),...
    fullfile(project_root,'visualization'));

generate_parameters();
load('parameters.mat','param');

yr=365*24*60*60;
checkpointer=param.checkpointer;
output_interval=param.output_interval;
checkpoint_interval=param.checkpoint_interval;
Nt=param.Nt;

alpha=param.alpha;
sina=sind(alpha);
cosa=cosd(alpha);
xsize=param.xsize;
ysize=param.ysize;
% The grid module owns Nx/Ny now: with stretching they follow from the core
% cell counts plus the arms, not from xsize/element_size. With no stretching
% requested it returns exactly the old uniform grid, Nx and Ny included.
g=build_stretched_grid(param);
Nx=g.Nx;
Ny=g.Ny;
N=(Nx+1)*(Ny+1)*2;
dx=g.dx_core;      % core spacing, for reporting and surface-station x
dy=g.dy_core;
dz=dy*sina;

if mod(Nx,2)==0
    error('BP3 requires an odd Nx so the fault lies on the central column.');
end

if ~checkpointer
    x=g.x;
    y=g.y;
    xp=g.xp;
    yp=g.yp;
    Xuy=y*cosa+xp;
    Yuy=y*sina+xp*0;
    Xux=yp*cosa+x;
    Yux=yp*sina+x*0;
    Xtau=y*cosa+x;
    Ytau=y*sina+x*0;
    Xsigma=yp(2:Ny,1)*cosa+xp(1,2:Nx);
    Ysigma=yp(2:Ny,1)*sina+xp(1,2:Nx)*0;
    z=y*sina;
    xd=y;
    save('coord.mat','x','y','xp','yp','Xtau','Xux','Xuy','Ytau',...
        'Yux','Yuy','Xsigma','Ysigma','z','xd','g');
else
    load('coord.mat');
end
% Spacings between adjacent nodes of each staggered field, for the stress
% recovery below. Uniform grid: every entry equals dx or dy.
sx_uy=diff(xp(:))';        % between uy columns -> at x positions,  1 x Nx
sy_ux=diff(yp(:));         % between ux rows    -> at y positions,  Ny x 1
sx_ux=diff(x(:))';         % between ux columns -> 1 x (Nx-1)
sy_uy=diff(y(:));          % between uy rows    -> (Ny-1) x 1

% Exact operators for the parts of the recovery that were hand-merged: the two
% node derivatives that were MATLAB gradient() (a two-cell secant inside, a
% two-point one-sided difference at the ends), the two midpoint-to-node
% interpolations that were 50/50 movmean, and the cell-centre-to-fault-node map
% whose END NODES were a copy of the nearest cell centre -- so sigma at the free
% surface came from half a cell down. Built once, because this feeds the
% timestep loop. recovery_operators.m records which direction of averaging was
% already exact and is deliberately left alone.
R=recovery_operators(x,y,xp,yp);

[rho,lambda,G,eta,K0,a,b,L,mu0,V0]=...
    build_layer_rsf(Nx,Ny,x,z,param);
lambdaP=movmean(movmean(lambda,2,2,'omitnan','Endpoints','discard'),...
    2,1,'Endpoints','discard');
GP=movmean(movmean(G,2,2,'omitnan','Endpoints','discard'),...
    2,1,'Endpoints','discard');

[sigman0,tau0,Pl0,Pr0]=initial_stress(0,0,0,0,alpha,K0,0,z,param);
P=zeros(Ny,1);

ux=zeros(Ny+1,Nx);
uy=zeros(Ny,Nx+1);
vx=zeros(Ny+1,Nx);
vy=zeros(Ny,Nx+1);
tauqs=zeros(Ny,Nx);
sigmaqs=zeros(Ny-1,Nx-1);
t=0;
dt=param.dt0;

fault_ix=(Nx+1)/2;
Gfault=(G(:,fault_ix-1)+G(:,fault_ix+1))/2;
% Local down-dip spacing at each fault node, not the core value: with a
% stretched y the fault can extend into the coarsened arm.
mfault=grid_metrics(g);
ksi=build_ksi(Gfault,L,mfault.uy.hyc(:),a,b,sigman0);
% At a domain ending at Wf, use the bottom fault node as the imposed
% deep-creep driver; all shallower nodes remain rate-and-state.
creep_mask=xd>=param.Wf;
rsf_mask=~creep_mask;

if ~checkpointer
    [sigma,tau,U,V,theta]=initial_fault(L,V0,param.Vinit,mu0,eta,...
        Ny,a,b,tau0,sigman0,0,param);
    V(creep_mask)=param.VL;
    LH=build_LH(lambda,G,sina,cosa,Nx,Ny,N,g,param);
    save('initiation.mat','LH','ksi','-v7.3');
    fid=fopen('output.txt','w+');
else
    load(['data_',int2str(checkpointer),'.mat']);
    load('initiation.mat','LH','ksi');
    fid=fopen('output.txt','a+');
end

global Um Vm taum sigmam Pm thetam dtm tm tm2
nout=ceil(Nt/output_interval);
Um=zeros(Ny,nout);
Vm=zeros(Ny,nout);
taum=zeros(Ny,nout);
sigmam=zeros(Ny,nout);
Pm=zeros(Ny,nout);
thetam=zeros(Ny,nout);
dtm=zeros(1,nout);
tm=zeros(1,nout);
tm2=zeros(1,nout);

surface_x=[-32e3,-16e3,-8e3,dx/2,-dx/2,8e3,16e3,32e3];
surface_side=[-1,-1,-1,1,-1,1,1,1];
surface_disp1=zeros(numel(surface_x),nout);
surface_disp2=zeros(numel(surface_x),nout);
surface_vel1=zeros(numel(surface_x),nout);
surface_vel2=zeros(numel(surface_x),nout);

dLH=decomposition(LH);
options=optimset('TolFun',param.friction_tolerance,'TolX',0);
toc;

for it=1:Nt
    drive=tauqs(rsf_mask,fault_ix)+tau0(rsf_mask);
    % Bracket the internal slip rate, which follows sign(Vp) = -motion_sign,
    % not motion_sign itself. Keying it off Vp keeps the bracket on the same
    % side as the drive even if a run overrides Vp directly.
    %
    % The width is twice the previous step's fastest rate-and-state point, not
    % a fixed ceiling. That contains the new root unless max|V| doubles in a
    % single step, whose measured worst case is 1.18 over ~700k steps of the
    % 100/50/25 m sets -- so a failure below is now a warning that the
    % velocity jumped, not an arbitrary cap being hit. The old fixed 10 m/s
    % stopped 50m_dip60_rev and yc40_dip60_rev mid-rupture at 9.97 m/s, where
    % Uphoff's dip-60 reverse itself reaches 7.8 m/s at the trace.
    % abs() first: the magnitude must stay positive, the branch carries sign.
    vb=2*max(abs(V(rsf_mask)));
    if param.Vp>0
        lower=zeros(nnz(rsf_mask),1);
        upper=zeros(nnz(rsf_mask),1)+vb;
    else
        lower=zeros(nnz(rsf_mask),1)-vb;
        upper=zeros(nnz(rsf_mask),1);
    end

    [Vrsf,~,exitflag]=bisection(@(VV) ...
        sigma(rsf_mask).*a(rsf_mask).*asinh(VV/(2*V0)...
        .*exp((mu0(rsf_mask)+b(rsf_mask)...
        .*log(V0*theta(rsf_mask)./L(rsf_mask)))./a(rsf_mask)))...
        +eta*VV-drive,lower,upper,zeros(nnz(rsf_mask),1),options);
    if any(exitflag<0) || any(~isfinite(Vrsf))
        error(['BP3 friction solve failed at iteration %d: the root left the '...
            'bracket of %g m/s, i.e. max|V| more than doubled in one step.'],...
            it,vb);
    end
    V(rsf_mask)=Vrsf;
    V(creep_mask)=param.VL;

    active=rsf_mask;
    speed=max(abs(V(active)),realmin);
    dt_stability=min(ksi(active).*L(active)./speed);
    dt=min([param.dtmax,param.dt_growth*dt,dt_stability,param.tfinal-t]);
    if dt<=0
        break;
    end

    q=abs(V)*dt./L;
    expo=q>1e-6;
    theta=expo.*(L./abs(V).*(1-exp(-q))+theta.*exp(-q))+...
        (~expo).*(theta+dt.*(1-abs(V).*theta./L));
    tau=tauqs(:,fault_ix)+tau0-eta*V;
    U=U+dt*V;

    RH=build_RH(lambda,G,sina,cosa,[],0,Nx,Ny,N,dx,dy,y,V,z,dz,param);
    S=dLH\RH;
    vpx=reshape(S(1:2:end),Ny+1,Nx+1);
    vpy=reshape(S(2:2:end),Ny+1,Nx+1);
    vy=vpy(1:Ny,:);
    vx=vpx(:,1:Nx);
    uy=uy+vy*dt;
    ux=ux+vx*dt;

    % Stress recovery. Every difference is divided by the spacing it actually
    % spans; the node derivatives and the midpoint-to-node interpolations come
    % from R (built above) instead of gradient() and movmean, which were exact
    % only on a uniform mesh -- and, at the boundary rows, not even there.
    %
    % Still plain differences, and correct as they stand: diff(uy,1,2)./sx_uy
    % and diff(ux,1,1)./sy_ux span a single interval between half-nodes. They
    % carry a (hp-hm)/4 position offset relative to the node they are assigned
    % to (about 1.3 m on the benchmark fault), which needs a wider stencil to
    % remove -- see recovery_operators.m.
    duxdx=R.Pyp2y*(ux*R.Dx.');       % d(ux)/dx at (y,x)
    duydy=(R.Dy*uy)*R.Pxp2x.';       % d(uy)/dy at (y,x)
    tauqs=G/sina.*(diff(uy,1,2)./sx_uy+...
        (1-2*cosa*cosa)*diff(ux,1,1)./sy_ux+...
        cosa*(duxdx-duydy));
    tauqs(:,fault_ix)=(tauqs(:,fault_ix-1)+tauqs(:,fault_ix+1))/2;

    % The two movmeans here go NODES -> MIDPOINTS, and xp/yp are the numerical
    % midpoints by construction, so 50/50 is exact. Left alone on purpose.
    sigmaqs=(lambdaP+2*GP).*diff(ux(2:Ny,:),1,2)./sx_ux+...
        lambdaP.*diff(uy(:,2:Nx),1,1)./sy_uy-...
        2*GP*cosa.*movmean(movmean(diff(ux,1,1)./sy_ux,2,2,...
        'Endpoints','discard'),2,1,'Endpoints','discard');
    sigmal=R.Psig*sigmaqs(:,(Nx-1)/2);
    sigmar=R.Psig*sigmaqs(:,(Nx+1)/2);
    sigma=sigman0-(sigmal+sigmar)/2;

    t=t+dt;
    fprintf(fid,['it=%d, t=%f yr, dt=%e s, max|V|=%e m/s, ',...
        'min sigma=%e Pa\n'],checkpointer+it,t/yr,dt,max(abs(V)),min(sigma));

    if mod(it,output_interval)==0
        write_memory(it,output_interval,U,V,tau,sigma,P,theta,dt,t,0,...
            tauqs,sigmaqs,uy,vy,ux,vx);
        io=it/output_interval;
        un=mean(ux(1:2,:),1);
        vn=mean(vx(1:2,:),1);
        us=interp1(x,un,surface_x,'linear');
        vs=interp1(x,vn,surface_x,'linear');
        uts=interp1(xp,uy(1,:),surface_x,'linear');
        vts=interp1(xp,vy(1,:),surface_x,'linear');
        if param.load_side_boundaries
            rigid_rate=zeros(size(surface_side));
        else
            rigid_rate=-surface_side*param.Vp/2;
        end
        surface_disp1(:,io)=(us+uts*cosa+rigid_rate*t*cosa)';
        surface_disp2(:,io)=(uts*sina+rigid_rate*t*sina)';
        surface_vel1(:,io)=(vs+vts*cosa+rigid_rate*cosa)';
        surface_vel2(:,io)=(vts*sina+rigid_rate*sina)';

        if param.live_plot && (io==1 || mod(io,param.live_plot_interval)==0)
            plot_bp3_live(param,xd,tm(1:io),Vm(:,1:io),Um(:,1:io),...
                taum(:,1:io),sigmam(:,1:io),dtm(1:io));
        end
    end

    if mod(it,checkpoint_interval)==0
        checkpoint_id=checkpointer+it;
        save(['data_',int2str(checkpoint_id),'.mat'],'U','V','tau',...
            'sigma','theta','dt','t','tauqs','sigmaqs','uy','vy','ux','vx');
        save('dataall.mat','Um','Vm','taum','sigmam','Pm','thetam',...
            'dtm','tm','tm2','-v7.3');
        disp(it);
        toc;
    end

    if t>=param.tfinal
        break;
    end
end

nwritten=floor(it/output_interval);
Um=Um(:,1:nwritten);
Vm=Vm(:,1:nwritten);
taum=taum(:,1:nwritten);
sigmam=sigmam(:,1:nwritten);
Pm=Pm(:,1:nwritten);
thetam=thetam(:,1:nwritten);
dtm=dtm(1:nwritten);
tm=tm(1:nwritten);
tm2=tm2(1:nwritten);
surface_disp1=surface_disp1(:,1:nwritten);
surface_disp2=surface_disp2(:,1:nwritten);
surface_vel1=surface_vel1(:,1:nwritten);
surface_vel2=surface_vel2(:,1:nwritten);
save(['data_BP3_QD_',datestr(now,'yyyy-mm-dd-HH-MM-ss'),'.mat'],...
    'Um','Vm','taum','sigmam','thetam','dtm','tm','xd','param',...
    'surface_x','surface_disp1','surface_disp2','surface_vel1',...
    'surface_vel2','-v7.3');
fclose(fid);
if nwritten>0
    write_bp3_outputs(param,xd,tm,Um,Vm,taum,sigmam,thetam,...
        surface_x,surface_disp1,surface_disp2,surface_vel1,surface_vel2);
end
