function g=build_stretched_grid(param)
%BUILD_STRETCHED_GRID Uniform-core / power-law-stretched grid with metrics.
%
%   g = BUILD_STRETCHED_GRID(param) returns the node and half-node coordinates
%   for the sheared BP3 grid together with the mapping metric ds/dzeta at both,
%   so build_LH can be written on a UNIFORM computational grid and carry the
%   geometry in the metric rather than in variable spacings.
%
%   The map (after Pranger's garnet, experiments/SEAS/BP4_QD) is uniform over
%   the first n_core intervals and then a power-law arm of exponent r that
%   lands exactly on the outer boundary:
%
%       s(z) = b*z                             z <= n_core/n_tot
%       s(z) = b*z + (H_tot - b)*chi^r         otherwise
%       b    = n_tot*H_core/n_core
%       chi  = n_tot*(z - n_core/n_tot)/(n_tot - n_core)
%
%   s is continuous at the seam (s = H_core there) and exact at s(1) = H_tot.
%   The metric is analytic, and because chi^(r-1) vanishes at the seam for
%   r >= 2 the CELL SIZE is continuous across it too; only its slope jumps.
%
%   param fields (all optional, defaults reproduce the uniform grid exactly):
%     element_size      spacing inside the uniform core
%     xsize, ysize      total extents: x spans +/-xsize/2, y spans [0 ysize]
%     x_core            half-width of the uniform core in x   (default xsize/2)
%     y_core            extent of the uniform core in y       (default ysize)
%     nx_stretch        cells in EACH x arm                   (default 0)
%     ny_stretch        cells in the y arm                    (default 0)
%     stretch_r         exponent                              (default 2)
%
%   With nx_stretch = ny_stretch = 0 this returns bit-for-bit the grid
%   Main_code.m builds today, which is the first regression test.

es=param.element_size;
xsize=param.xsize;
ysize=param.ysize;
r=getfield_default(param,'stretch_r',2);
x_core=getfield_default(param,'x_core',xsize/2);
y_core=getfield_default(param,'y_core',ysize);
nxs=getfield_default(param,'nx_stretch',0);
nys=getfield_default(param,'ny_stretch',0);

% The y seam and the fault: history, not a live constraint.
%
% Coarsening ON the fault was once a wrong answer rather than a resolution
% trade-off, because build_LH's hand-merged centred first and mixed derivatives
% were only first-order once hm ~= hp, which pinned truncation error at the
% seam and let it accumulate with the load. Measured at dip 60, dx 100 m,
% Wf = 40 km, seam at y_core = 20 km, ON THE PRE-2e OPERATOR:
%
%   nys   max cell in Wf   sigma error on the fault   first event
%   250   100 m   (uniform)          --               133.2 yr
%   133   251 m                  +0.8 MPa             147.5 yr
%    90   410 m                  +1.5 MPa             159.8 yr
%    54   739 m                  +2.6 MPa             179.0 yr
%    36  1140 m                  +3.6 MPa             195.3 yr
%    27  1538 m                  +4.4 MPa             none in 200 yr
%
% Increment 2e split those coefficients and took the nys=36 nucleation error
% from 62 yr to 0.05 yr. The five-case benchmark now runs at y_core = 20 km and
% matches the references to within a few years, and on the fixed operator that
% configuration beats the release code's y_core = 40 km on every measure.
%
% The old error() guard, and the two warnings that replaced it, are therefore
% GONE. They quoted the pre-2e magnitudes above as if they were current, and
% they fired on every production run. The same applies to the "y_core = Wf
% exactly is not enough, leave a 5 km buffer" rule: that came from a +0.11 ->
% +1.76 MPa measurement on the unfixed operator.
%
% What remains true is the ordinary statement, and it needs no trap: cells
% below y_core are coarser than the core, so y_core sets deep-fault
% resolution. Choose it deliberately.

nxc=round(x_core/es);        % core cells per x arm
nyc=round(y_core/es);        % core cells in y
nxt=nxc+nxs;                 % cells per x arm, total
nyt=nyc+nys;

Nx=2*nxt+1;
Ny=nyt+1;
if mod(Nx,2)==0
    error('build_stretched_grid:evenNx','Nx must be odd; got %d.',Nx);
end

% Computational coordinates. x is symmetric on [-1,1], y one-sided on [0,1].
% x and y are THE nodes; xp and yp are their numerical midpoints.
%
% They used to be the map evaluated at computational half-nodes, which is the
% same thing on a uniform grid but NOT on a stretched one: for a quadratic map
%     s(z+h/2) - (s(z)+s(z+h))/2 = -(h^2/8)*s''
% so the half-node sits off the midpoint by an amount proportional to the
% curvature of the grading -- 0 in the uniform core, a constant -4.13 m through
% the arm at y_core = 20 km / nys = 36, stepping at the seam. Because the
% offset is a STEP, diff(yp) departs from the divisor hyc that build_LH uses by
% -3.8 % in the single cell at the seam, exactly where the normal-stress bump
% begins. Defining the half-nodes AS the midpoints makes diff(yp) == hyc
% identically, so the operator and the stress recovery in Main_code (which
% divides by diff(yp)) agree by construction instead of approximately.
%
% On a uniform grid this is a no-op TO ROUNDOFF, not bit-for-bit: the map is
% linear there so midpoint and half-node coincide analytically, but the two
% arithmetic routes differ by ~1 ulp (7e-12 m on a 45 km domain, 2e-16
% relative). That is far inside test_LH_equivalence's 1e-12 tolerance, and the
% m.uniform guard in build_LH still forces the graded centre weights to exactly
% zero so nnz cannot drift.
dzx=2/(Nx-1);
dzy=1/(Ny-1);
zx=linspace(-1,1,Nx);
zy=linspace(0,1,Ny)';

[x,bx]=map_sym(zx,x_core,xsize/2,nxc,nxs,r);
[y,by]=map_one(zy,y_core,ysize,nyc,nys,r);
xp=midpoints(x);
yp=midpoints(y);

g=struct();
g.x=x; g.xp=xp; g.y=y; g.yp=yp;
% bx,by are ds/dzeta at the NODES. The half-node metric is gone with the
% half-node coordinates: nothing consumed g.bxp/g.byp.
g.bx=bx; g.by=by;
g.dzx=dzx; g.dzy=dzy;                     % uniform computational spacings
g.Nx=Nx; g.Ny=Ny;
g.nxc=nxc; g.nyc=nyc; g.nxs=nxs; g.nys=nys; g.r=r;
g.stretched=(nxs>0)||(nys>0);
% Physical spacings, for reporting and for build_ksi (which needs the local
% on-fault value, not a single global dy).
g.dx=diff(x); g.dy=diff(y);
g.dy_core=es; g.dx_core=es;
end


function p=midpoints(v)
%MIDPOINTS Staggered nodes: the numerical midpoints of v, with one ghost
%   half-cell beyond each end. Orientation is preserved, so a row of x nodes
%   gives a row of xp and a column of y nodes gives a column of yp.
row=isrow(v);
v=v(:).';
p=[v(1)-0.5*(v(2)-v(1)), 0.5*(v(1:end-1)+v(2:end)), v(end)+0.5*(v(end)-v(end-1))];
if ~row
    p=p.';
end
end

function v=getfield_default(s,name,d)
if isfield(s,name) && ~isempty(s.(name))
    v=s.(name);
else
    v=d;
end
end


function [s,beta]=map_one(z,Hc,Ht,nc,ns,r)
%MAP_ONE One-sided map on z in [0,1], extendable a little outside for ghosts.
if ns==0 || Ht==Hc
    s=Ht*z;
    beta=Ht*ones(size(z));
    return;
end
nt=nc+ns;
b=nt*Hc/nc;
zs=nc/nt;
chi=nt*(z-zs)/(nt-nc);
inner=z<=zs;
s=zeros(size(z));
beta=zeros(size(z));
s(inner)=b*z(inner);
beta(inner)=b;
% Chain rule directly: d/dz of (Ht-b)*chi^r is (Ht-b)*r*chi^(r-1)*dchi/dz with
% dchi/dz = nt/(nt-nc). garnet writes this as b*(1+(Ht/alpha-1)*r*chi^(r-1))
% with alpha = Ht*Hc*(nt-nc)/(nc*(Ht-Hc)); the two are algebraically identical.
s(~inner)=b*z(~inner)+(Ht-b)*chi(~inner).^r;
beta(~inner)=b+(Ht-b)*r*chi(~inner).^(r-1)*nt/(nt-nc);
end


function [s,beta]=map_sym(z,Hc,Ht,nc,ns,r)
%MAP_SYM The same map folded about the origin, for z in [-1,1].
%   s(z) = sign(z)*f(|z|) so s is odd; beta = f'(|z|) is even, since the two
%   sign factors from the chain rule cancel.
[a,ba]=map_one(abs(z),Hc,Ht,nc,ns,r);
s=sign(z).*a;
beta=ba;
end
