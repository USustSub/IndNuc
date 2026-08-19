function m=grid_metrics(g)
%GRID_METRICS Local spacings the stencils need, one entry per node.
%
%   The operator is assembled on a UNIFORM computational grid; all geometry
%   enters here. Two maps compose to make the physical grid:
%
%     shear      X = x + y*cos(psi),  Z = y*sin(psi)
%     stretch    x = X(zeta),         y = Y(eta)
%
%   and their Jacobians multiply,
%
%     d(X,Z)/d(zeta,eta) = [[1 cos(psi)];[0 sin(psi)]] * diag(bx,by),
%
%   so the shear algebra already hard-coded in build_LH is untouched and only
%   d/dx -> (1/bx) d/dzeta, d/dy -> (1/by) d/deta change. Because the stretch
%   is separable, the spacings depend on ix or iy alone, never on both -- that
%   is what keeps this a vector per direction rather than a field per node.
%
%   Returned, for each of the two staggered variables:
%     hxm,hxp   distance to the x-neighbour below/above  (indexed by ix)
%     hym,hyp   distance to the y-neighbour below/above  (indexed by iy)
%     hxc,hyc   the centred spacing (hm+hp)/2, for flux-form second derivatives
%
%   uy lives at (xp, y); ux lives at (x, yp) -- hence separate sets.
%
%   A GENERAL curvilinear mesh would break separability and need the full
%   metric tensor with cross terms; the intent is that only this function
%   changes when that day comes, and build_LH keeps consuming m.<var>.<h>.

m=struct();
m.uy=spacings(g.xp,g.y,g.Nx+1,g.Ny);
m.ux=spacings(g.x,g.yp,g.Nx,g.Ny+1);
m.uniform=~g.stretched;
end


function s=spacings(xv,yv,nx,ny)
s.hxm=pad_diff(xv(:),nx);
s.hxp=[s.hxm(2:end);s.hxm(end)];
s.hym=pad_diff(yv(:),ny);
s.hyp=[s.hym(2:end);s.hym(end)];
s.hxc=(s.hxm+s.hxp)/2;
s.hyc=(s.hym+s.hyp)/2;
end


function h=pad_diff(v,n)
%PAD_DIFF Backward spacing at each of n nodes; the first is copied from the
%   second so boundary rows never see a NaN. Boundary rows are algebraic
%   (Dirichlet or one-sided) and do not use a two-sided spacing anyway.
d=diff(v(:));
if numel(d)<n-1
    d(end+1:n-1)=d(end);
end
h=[d(1);d(1:n-1)];
end
