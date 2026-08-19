function R=recovery_operators(x,y,xp,yp)
%RECOVERY_OPERATORS Exact staggered-grid operators for the stress recovery.
%
%   R = RECOVERY_OPERATORS(x,y,xp,yp) returns sparse matrices that replace the
%   hand-merged `gradient` calls and 50/50 `movmean` averages in Main_code's
%   stress recovery. Each one is exact where the thing it replaces was only
%   first order (or, at the free surface, zeroth order) on a graded mesh.
%
%   Applied as, with ux (Ny+1 x Nx) at (yp,x) and uy (Ny x Nx+1) at (y,xp):
%
%       d(ux)/dx at (y ,x)  =  R.Pyp2y*(ux*R.Dx.')
%       d(uy)/dy at (y ,x)  =  (R.Dy*uy)*R.Pxp2x.'
%       sigma    at  y      =  R.Psig*sigmaqs(:,column)
%
%   WHY EACH ONE EXISTS
%
%   Dy, Dx -- first derivative at the NODES.  MATLAB's gradient(f,coord) uses
%   the two-cell secant (f(i+1)-f(i-1))/(coord(i+1)-coord(i-1)) (see
%   toolbox/matlab/datafun/gradient.m), which drops the f(i) term entirely and
%   is only first order once hm ~= hp; its error is ((hp-hm)/2)*f''. That is
%   the same defect increment 2e removed from build_LH, left behind here. The
%   exact three-point weights are what FDWEIGHTS returns, and they are
%   identical to the distance-weighted average of the two adjacent staggered
%   slopes -- so this is consistent with the staggered operator, not a
%   competing discretisation.
%
%   At the two end nodes gradient uses a TWO-point one-sided difference, first
%   order on any mesh including a uniform one. These use the three-point
%   one-sided stencil instead, so the change is NOT a no-op on a uniform grid:
%   it is second order at the boundary where the old code was first order.
%
%   Pyp2y, Pxp2x -- staggered midpoints back onto the nodes. Direction matters
%   and only one direction was wrong:
%
%     nodes -> midpoints   xp = midpoints(x) BY CONSTRUCTION, so a 50/50
%                          average is EXACT. Main_code's lambdaP/GP and the
%                          cos(alpha) term of sigmaqs go this way and are
%                          left alone.
%     midpoints -> nodes   y(j) is NOT the midpoint of yp(j), yp(j+1) on a
%                          graded mesh -- it sits (hp-hm)/4 away. A 50/50
%                          average therefore samples the field at the wrong
%                          place, first order in the grading.
%
%   The end rows come out at exactly 1/2 on their own, because midpoints()
%   builds the ghost half-cells symmetrically about the end nodes.
%
%   Psig -- cell-centred sigma onto the fault nodes. sigmaqs row k lives at
%   yp(k+1), so the interior nodes need the same distance weighting as above.
%   The ends needed more: the old code COPIED the nearest cell centre, so
%   sigma at the free-surface fault node was taken from half a cell down --
%   50 m at 100 m resolution, a zeroth-order error and by far the largest of
%   the three, present on a uniform mesh too. Here both ends are one-sided
%   quadratic extrapolations from the three nearest centres.
%
%   NOT CHANGED, deliberately: the single-interval staggered differences
%   diff(uy,1,2)./sx_uy and diff(ux,1,1)./sy_ux. These land at the midpoint of
%   two half-nodes, which is (hp-hm)/4 from the node they are assigned to, so
%   they carry the same kind of first-order position offset -- about 1.3 m on
%   the benchmark fault, worth ~0.1 %. Removing it needs a wider stencil than
%   two points, which changes the recovery's footprint; that is a separate
%   decision, not a bug fix.

x=x(:).'; xp=xp(:).'; y=y(:); yp=yp(:);
Nx=numel(x); Ny=numel(y);
if numel(xp)~=Nx+1 || numel(yp)~=Ny+1
    error('recovery_operators:size', ...
        'expected numel(xp)=Nx+1 and numel(yp)=Ny+1; got %d,%d for Nx=%d,Ny=%d.', ...
        numel(xp),numel(yp),Nx,Ny);
end

R=struct();
R.Dy=node_derivative(y);
R.Dx=node_derivative(x);
R.Pyp2y=mid_to_node(y);
R.Pxp2x=mid_to_node(x);
R.Psig=centres_to_nodes(y,yp);
end


function D=node_derivative(v)
%NODE_DERIVATIVE First derivative at each node from its own three-point
%   stencil: centred inside, one-sided at the two ends.
v=v(:); n=numel(v);
if n<3
    error('recovery_operators:tooFewNodes','need at least 3 nodes; got %d.',n);
end
ii=zeros(3*n,1); jj=ii; vv=ii; c=0;
for k=1:n
    if k==1
        idx=[1 2 3];
    elseif k==n
        idx=[n-2 n-1 n];
    else
        idx=[k-1 k k+1];
    end
    w=fdweights(v(k),v(idx),1);
    ii(c+1:c+3)=k; jj(c+1:c+3)=idx; vv(c+1:c+3)=w; c=c+3;
end
D=sparse(ii,jj,vv,n,n);
end


function P=mid_to_node(v)
%MID_TO_NODE Linear interpolation from the n+1 midpoints of v onto its n
%   nodes, weighted by the distance to each neighbour rather than 50/50.
v=v(:); n=numel(v);
h=diff(v);
hh=[h(1);h;h(end)];                 % n+1 spacings, ghost half-cells mirrored
w=hh(2:n+1)./(hh(1:n)+hh(2:n+1));   % weight on midpoint j; 1/2 at both ends
P=sparse([(1:n).';(1:n).'],[(1:n).';(2:n+1).'],[w;1-w],n,n+1);
end


function P=centres_to_nodes(v,vp)
%CENTRES_TO_NODES sigmaqs rows -- at vp(2:n), i.e. the n-1 cell centres --
%   onto the n nodes v. Distance weighted inside, one-sided quadratic
%   extrapolation at the two ends, where no centre brackets the node.
v=v(:); vp=vp(:); n=numel(v); m=n-1;
if m<3
    error('recovery_operators:tooFewCentres', ...
        'need at least 3 cell centres for the end extrapolations; got %d.',m);
end
h=diff(v);
ii=[]; jj=[]; vv=[];
we=fdweights(v(1),vp(2:4),0);       % free surface, from the three nearest
ii=[ii;1;1;1]; jj=[jj;1;2;3]; vv=[vv;we(:)];
for k=2:n-1
    w=h(k)/(h(k-1)+h(k));           % weight on the centre ABOVE the node
    ii=[ii;k;k]; jj=[jj;k-1;k]; vv=[vv;w;1-w];
end
% Columns 1..m of sigmaqs sit at vp(2..n), so the three centres nearest the
% deep end are vp(n-2:n) -- NOT vp(n+1), which is the ghost half-node past the
% last row and carries no sigma at all.
wl=fdweights(v(n),vp(n-2:n),0);     % deep end, from the three nearest
ii=[ii;n;n;n]; jj=[jj;m-2;m-1;m]; vv=[vv;wl(:)];
P=sparse(ii,jj,vv,n,m);
end
