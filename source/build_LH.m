function LH=build_LH(lambda,G,sina,cosa,Nx,Ny,N,g,param)
%BUILD_LH Elastostatic stiffness on the sheared, optionally stretched grid.
%
%   Geometry enters only through grid_metrics. The shear Jacobian and the
%   stretch Jacobian multiply, so every cos(alpha)/sin(alpha) term below is
%   unchanged and only the spacings become position dependent: dx and dy are
%   reassigned per node at the top of each variable's block. On a uniform grid
%   they take the same value everywhere and this reduces exactly to the
%   original constant-spacing operator, which is the regression test.
%
%   ACCURACY ON A STRETCHED GRID -- READ THIS BEFORE TRUSTING A STRETCHED RUN.
%   The second-derivative stencils use the asymmetric FLUX form rather than the
%   symmetric -2/+1/+1 weights (increment 2b), and those are correct and
%   conservative on any spacing.
%
%   The centred FIRST and MIXED derivatives carried an extra error term. Over
%   unequal intervals,
%
%       (u+ - u-)/(hm+hp) = u' + (hp-hm)/2 * u'' + ...
%
%   For a smooth map hp-hm = dzeta^2*s'', so that form is still SECOND-ORDER
%   CONVERGENT -- but its leading error coefficient is proportional to s'', the
%   CURVATURE OF THE MESH GRADING. Measured consequences on the dip-60
%   y-coarsening series, with the seam inside Wf (y_core = 20 km):
%     - the error is pinned at the seam, because stretch_r = 2 makes s'' jump
%       there, and grows with the grading ratio;
%     - sigma on the rate-and-state fault rose by up to +5 MPa and nucleation
%       was delayed by 62 yr, or suppressed entirely at nys <= 27;
%     - stretch_r = 3 makes s'' continuous at the seam and removed only ~1/3
%       of it, confirming the map exponent is not the remedy;
%     - every one of these errors is EXACTLY zero when hm == hp, which is why
%       test_LH_equivalence passes at 3e-15 and cannot see any of it, and they
%       all carry a cosa factor, so the dip-90 symmetry regression is blind too.
%
%   INCREMENT 2d. The exact three-point non-uniform weights are exact for
%   quadratics, so the u'' term is identically absent. The mixed derivatives
%   become the full 3x3 tensor product instead of four corners. Measured: this
%   changed the seam bump by ~1 %. It is correct, but it was NOT the problem.
%
%   INCREMENT 2e -- THE BIG ONE. The cross-variable coupling is two operators,
%       u2 equation:  (lambda+G) ( d1 d2 - cos(alpha) d1^2 ) u1
%       u1 equation:  (lambda+G) ( d1 d2 - cos(alpha) d2^2 ) u2
%   (Shang, solver.tex). They used to be SUMMED into one coefficient per node,
%   which is free on a uniform grid but hides that each needs its own
%   spacing-dependent weights. The second-derivative half is evaluated at a
%   HALF NODE from four columns/rows, and its weights were hard-coded as
%   (1,-1,-1,1)/(2h^2) -- the uniform set. On the nys = 36 seam spacing the
%   correct weights are not even antisymmetric and the hard-coded set is 25 %
%   wrong on a pure quadratic. That is an order of magnitude larger than
%   anything increment 2d touched, and it is where the y-coarsening bump lives.
%   Now split, with source/fdweights.m supplying the exact weights. Verified to
%   reproduce the old merged coefficients to 3e-16 on a uniform grid.
%
%   All stencils now come from ONE helper, source/fdweights.m -- weights on
%   arbitrarily spaced nodes by polynomial exactness. onesided3 and centred3
%   are thin wrappers over it.
%
%   INCREMENT 2c: the last two row families that assumed LOCALLY uniform
%   spacing are now derived for arbitrary spacing.
%
%     - Fault traction continuity (ix==fault_ix normal, ix==fault_ix+1 shear)
%       differenced across the fault WITHOUT dividing by a spacing, which is
%       only valid when the cells either side are equal. Each side is now
%       divided by its own interval and the row rescaled by the centred
%       spacing, so the row keeps its previous magnitude.
%     - The free-surface rows used the uniform one-sided weights
%       (-3,4,-1)/(2dy). They now use the general three-point one-sided
%       weights for spacings (h1,h2), which collapse to those on a uniform
%       grid.
%
%   Both reduce EXACTLY to the previous expressions when h1==h2, so
%   test_LH_equivalence still pins them against the frozen reference.
    m=grid_metrics(g);
    % One-sided three-point d/dy weights at the free surface, for the two
    % staggered y grids. uy sits on g.y (first node ON the surface); ux sits on
    % g.yp (nodes straddling it).
    [wy1,wy2,wy3]=onesided3(m.uy.hyp(1),m.uy.hyp(2));
    hux1=m.ux.hyp(1);   % ux node spacing across the surface
    % The core-width guard is now belt-and-braces rather than load-bearing --
    % the stencils above are correct on a non-uniform mesh -- but a core
    % narrower than the stencils it contains is pathological, so keep failing
    % loudly.
    if g.nxs>0 && g.nxc<4
        error('build_LH:narrowCoreX', ...
            ['x core is %d cells per side; the fault stencils reach 2 cells, ' ...
             'so a stretched x grid needs at least 4.'],g.nxc);
    end
    if g.nys>0 && g.nyc<4
        error('build_LH:narrowCoreY', ...
            ['y core is %d cells; the one-sided free-surface stencil reaches ' ...
             '2 cells, so a stretched y grid needs at least 4.'],g.nyc);
    end
    % TRIPLETS PER ROW, not stencil width -- the stencil is 17 NODES either way
    % (9 own-variable + 8 cross-variable), as in the stencil figure.
    %   before 2d/2e : 5 second-derivative + 4 mixed corners + 8 coupling = 17
    %   after        : 5 + 9 (mixed 3x3) + 12 (coupling split) = 26
    % The extra triplets land on columns that are already in the stencil, and
    % sparse() sums duplicates, so nnz per row stays 17. 28 gives a little
    % slack over the counted 26 because the overflow check below is fatal.
    ntrip=28;
    I=zeros(ntrip*N,1);
    J=zeros(ntrip*N,1);
    LL=zeros(ntrip*N,1);
    ik=1;
    for ix=1:Nx+1
        for iy=1:Ny+1
            kux=((ix-1)*(Ny+1)+iy-1)*2+1;
            kuy=kux+1;
            if (iy<Ny+1)
                % uy lives at (xp, y). dx,dy are the centred spacings, which is
                % what first derivatives need; hxm/hxp/hym/hyp are the one-sided
                % neighbour distances, which second derivatives need.
                dx=m.uy.hxc(ix); dy=m.uy.hyc(iy);
                hxm=m.uy.hxm(ix); hxp=m.uy.hxp(ix);
                hym=m.uy.hym(iy); hyp=m.uy.hyp(iy);
                if (ix==1) % far-left boundary uy=0
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
                elseif (ix==Nx+1) % far-right boundary uy=0
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=1;ik=ik+1;
                elseif (iy==1 && ix~=(Nx+1)/2 && ix~=(Nx+1)/2+1)
                    % free surface: sigma_zz=0
                    % Both fault columns are excluded here so each falls
                    % through to its own branch below -- the slip-rate jump at
                    % ix==(Nx+1)/2 and shear-traction continuity at
                    % ix==(Nx+1)/2+1 -- the same pairing the interior rows use.
                    % Giving the surface row two sigma_zz=0 rows instead left
                    % the trace unconstrained: its elastic slip then crept away
                    % from the friction-solved value without bound, reaching
                    % 50x by 48 yr at dx=100 m. Applying the jump to the minus
                    % side alone fixes that but breaks the x -> -x reflection
                    % symmetry; both sides together fix it and keep the dip-90
                    % antisymmetry exact. build_RH must loop iy=1:Ny.
                    % (dx/G)*sigma_yy=0 at x^2=0:
                    % d2(u2)+lambda/(lambda+2G)d1(u1)
                    % -2G*cos(alpha)/(lambda+2G)d1(u2)=0.
                    % The normal derivative is one-sided at the surface;
                    % the other terms are evaluated at x^2=0.
                    lambda_surface=mean(lambda(1,[ix-1,ix]),'omitnan');
                    G_surface=mean(G(1,[ix-1,ix]),'omitnan');
                    normal_scale=dx/G_surface;
                    nsc=normal_scale*(lambda_surface+2*G_surface);
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=nsc*wy1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=nsc*wy2;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+4;LL(ik)=nsc*wy3;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=normal_scale*G_surface*cosa/dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=-normal_scale*G_surface*cosa/dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kux;LL(ik)=normal_scale*lambda_surface/(2*dx);ik=ik+1;
                    I(ik)=kuy;J(ik)=kux+2;LL(ik)=normal_scale*lambda_surface/(2*dx);ik=ik+1;
                    I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=-normal_scale*lambda_surface/(2*dx);ik=ik+1;
                    I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-normal_scale*lambda_surface/(2*dx);ik=ik+1;
                elseif (iy==Ny && ix==(Nx+1)/2) % creeping bottom fault jump
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
                elseif (iy==Ny) % far-bottom boundary uy=0 away from fault
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2) % fault diff(vy)=Vy
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=-2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2*(Ny+1)*2;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2+1) % fault continous tau?
                    GA=G(iy,ix-2);GB=G(iy,ix);GC=(GA+GB)/2;
                    % Each side's d(uy)/dx is a one-sided difference over ITS
                    % OWN interval -- the minus side spans columns ix-2..ix-1,
                    % so its spacing is hxm(ix-1), not hxm(ix). Rescaled by the
                    % centred spacing dx so the row keeps the magnitude it had
                    % when both sides were equal (where sm=sp=1).
                    sm=dx/m.uy.hxm(ix-1);
                    sp=dx/m.uy.hxp(ix);
                    I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2;LL(ik)=sm*GA/GC;ik=ik+1; % uy3
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=-sm*GA/GC;ik=ik+1; % uy4
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-sp*GB/GC;ik=ik+1; % uy5
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=sp*GB/GC;ik=ik+1; % uy6
                    % The cos(alpha)*d(uy)/dy terms, one per side. cA and cB
                    % already carry the 1/(2dy), so the bracketed weights are
                    % just the derivative stencil: centred (-1,0,+1) in the
                    % interior, one-sided (-3,+4,-1) at the free surface where
                    % iy-1 does not exist. Duplicate (I,J) pairs are summed by
                    % sparse(), so overlapping with uy4/uy5 above is fine.
                    cA=cosa/2/dy*dx*GA/GC;    % minus side, column ix-1
                    cB=-cosa/2/dy*dx*GB/GC;   % plus side, column ix
                    % Without the 1/(2dy): the one-sided weights carry their own
                    % spacing, so cA0*wy1 == -3*cA when h1==h2.
                    cA0=cosa*dx*GA/GC;
                    cB0=-cosa*dx*GB/GC;
                    if (iy==1)
                        I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=cA0*wy1;ik=ik+1;
                        I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cA0*wy2;ik=ik+1;
                        I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+4;LL(ik)=cA0*wy3;ik=ik+1;
                        I(ik)=kuy;J(ik)=kuy;LL(ik)=cB0*wy1;ik=ik+1;
                        I(ik)=kuy;J(ik)=kuy+2;LL(ik)=cB0*wy2;ik=ik+1;
                        I(ik)=kuy;J(ik)=kuy+4;LL(ik)=cB0*wy3;ik=ik+1;
                    else
                        I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cA;ik=ik+1; % uy1
                        I(ik)=kuy;J(ik)=kuy-2;LL(ik)=-cB;ik=ik+1; % uy2
                        I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cA;ik=ik+1; % uy7
                        I(ik)=kuy;J(ik)=kuy+2;LL(ik)=cB;ik=ik+1; % uy8
                    end
                    I(ik)=kuy;J(ik)=kux;LL(ik)=cosa/2*GB/GC;ik=ik+1; % ux3
                    I(ik)=kuy;J(ik)=kux+2;LL(ik)=cosa/2*GB/GC;ik=ik+1; % ux6
                    I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=-cosa/2*(GA+GB)/GC+(1-2*cosa*cosa)/dy*dx*(GA-GB)/GC;ik=ik+1; % ux2
                    I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-cosa/2*(GA+GB)/GC+(1-2*cosa*cosa)/dy*dx*(GB-GA)/GC;ik=ik+1; % ux5
                    I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2;LL(ik)=cosa/2*GA/GC;ik=ik+1; % ux1
                    I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2+2;LL(ik)=cosa/2*GA/GC;ik=ik+1; % ux4
%                                     I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2;LL(ik)=(1-2*cosa*cosa)/dy*dx+cosa/2;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2+2;LL(ik)=-(1-2*cosa*cosa)/dy*dx+cosa/2;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux;LL(ik)=-(1-2*cosa*cosa)/dy*dx+cosa/2;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux+2;LL(ik)=(1-2*cosa*cosa)/dy*dx+cosa/2;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=-cosa;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-cosa;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy-2;LL(ik)=cosa/2/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-cosa/2/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cosa/2/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cosa/2/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2;LL(ik)=1;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=-1;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2;LL(ik)=(1-2*cosa*cosa)/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2+2;LL(ik)=-(1-2*cosa*cosa)/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy;LL(ik)=-1;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux;LL(ik)=-(1-2*cosa*cosa)/dy*dx;ik=ik+1;
%                                     I(ik)=kuy;J(ik)=kux+2;LL(ik)=(1-2*cosa*cosa)/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kux+(Ny+1)*2;LL(ik)=cosa/4;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kux+(Ny+1)*2+2;LL(ik)=cosa/4;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=-cosa/2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-cosa/2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kux-3*(Ny+1)*2;LL(ik)=cosa/4;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kux-3*(Ny+1)*2+2;LL(ik)=cosa/4;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2-2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2+2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
%                 elseif (iy==1)
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-1;ik=ik+1;
%                 elseif (iy==Ny)
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-2;LL(ik)=-1;ik=ik+1;
                else % uy-Navier
                    GA=G(iy,ix);GB=G(iy+1,ix);GC=(GA+GB)/2;
                    lambdaA=lambda(iy,ix);lambdaB=lambda(iy+1,ix);lambdaC=(lambdaA+lambdaB)/2;
                    % Second derivatives in asymmetric form (increment 2b). The
                    % row scale stays dx^2 = hxc^2, so the cross terms below are
                    % unchanged. Uniform grid: hxm=hxp=dx and hym=hyp=dy, giving
                    % back -2 / +1 / +1 and dx^2/dy^2 exactly.
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-dx*(1/hxm+1/hxp) ...
                        -dx*dx/(GC*dy)*((lambdaA+2*GA)/hym+(lambdaB+2*GB)/hyp);ik=ik+1; % uy5
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=dx/hxm;ik=ik+1; % uy4
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=dx/hxp;ik=ik+1; % uy6
                    I(ik)=kuy;J(ik)=kuy-2;LL(ik)=dx*dx*(lambdaA+2*GA)/(GC*dy*hym);ik=ik+1; % uy2
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=dx*dx*(lambdaB+2*GB)/(GC*dy*hyp);ik=ik+1; % uy8
                    % INCREMENT 2d: mixed derivative as an EXACT non-uniform
                    % tensor product. The four-corner form is the plain centred
                    % difference in each direction, which carries an
                    % O(hp-hm)*f'' error whose coefficient is the CURVATURE of
                    % the mesh grading -- discontinuous at the seam, which is
                    % what pinned the normal-stress bump on the fault.
                    % centred3 removes that term. The stencil becomes 3x3; the
                    % five new weights are identically zero on a uniform grid,
                    % so the frozen-reference regression is untouched.
                    [axm,ax0,axp]=centred3(hxm,hxp);
                    [aym,ay0,ayp]=centred3(hym,hyp);
                    if m.uniform, ax0=0; ay0=0; end
                    sx=hxm+hxp;   % restores the previous +/-1 x-difference scale
                    Mv=[cosa*dx*(lambdaC+GC+2*GA)/GC, ...
                        cosa*dx*(lambdaC+GC+2*GC)/GC, ...
                        cosa*dx*(lambdaC+GC+2*GB)/GC];
                    wy=[aym ay0 ayp]; wx=[axm ax0 axp];
                    for ey=1:3
                        for ex=1:3
                            cxy=-Mv(ey)/2*wy(ey)*sx*wx(ex);
                            if cxy~=0
                                I(ik)=kuy;
                                J(ik)=kuy+(ex-2)*(Ny+1)*2+(ey-2)*2;
                                LL(ik)=cxy;ik=ik+1;   % uy1/3/7/9, +2/4/5/6/8 if graded
                            end
                        end
                    end
                    if (ix==2 || ix==Nx)
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaA+GC)/GC;ik=ik+1; % ux2
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-1/dy*dx*(lambdaB+GC)/GC;ik=ik+1; % ux6
                        I(ik)=kuy;J(ik)=kux;LL(ik)=-1/dy*dx*(lambdaA+GC)/GC;ik=ik+1; % ux3
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=1/dy*dx*(lambdaB+GC)/GC;ik=ik+1; % ux7
                    else
                        % INCREMENT 2e. The u1 coupling in the u2 equation is
                        %     (lambda+G) * ( d1 d2 - cos(alpha) d1^2 ) u1
                        % i.e. TWO operators. They used to be summed into one
                        % coefficient at ux2/ux3/ux6/ux7, which is free on a
                        % uniform grid but hides the fact that each needs its
                        % own spacing-dependent weights. Split here.
                        %
                        % (i) the mixed derivative d1 d2 u1: staggered, one
                        % interval each way, so it is already right on any mesh.
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaA+GC)/GC;ik=ik+1; % ux2
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-1/dy*dx*(lambdaB+GC)/GC;ik=ik+1; % ux6
                        I(ik)=kuy;J(ik)=kux;LL(ik)=-1/dy*dx*(lambdaA+GC)/GC;ik=ik+1; % ux3
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=1/dy*dx*(lambdaB+GC)/GC;ik=ik+1; % ux7
                        % (ii) -cos(alpha)*d1^2 u1: a SECOND derivative at the
                        % half node xp(ix), from the four ux columns
                        % x(ix-2..ix+1), averaged over the two y levels. The old
                        % code hard-coded the uniform weights (1,-1,-1,1)/(2dx^2)
                        % as +/-cosa*(lam+G)/GC/4. Those are 25 % wrong on the
                        % nys=36 seam spacing and are not even antisymmetric
                        % there. fdweights gives the exact set for any spacing
                        % and reproduces (1,-1,-1,1)/(2dx^2) when uniform.
                        wd2x=fdweights(g.xp(ix),[g.x(ix-2) g.x(ix-1) g.x(ix) g.x(ix+1)],2);
                        d2cA=-cosa*(lambdaA+GA)/GC*dx*dx/2;   % level iy
                        d2cB=-cosa*(lambdaB+GB)/GC*dx*dx/2;   % level iy+1
                        I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2;LL(ik)=d2cA*wd2x(1);ik=ik+1; % ux1
                        I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2+2;LL(ik)=d2cB*wd2x(1);ik=ik+1; % ux5
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=d2cA*wd2x(2);ik=ik+1; % ux2
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=d2cB*wd2x(2);ik=ik+1; % ux6
                        I(ik)=kuy;J(ik)=kux;LL(ik)=d2cA*wd2x(3);ik=ik+1; % ux3
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=d2cB*wd2x(3);ik=ik+1; % ux7
                        I(ik)=kuy;J(ik)=kux+(Ny+1)*2;LL(ik)=d2cA*wd2x(4);ik=ik+1; % ux4
                        I(ik)=kuy;J(ik)=kux+(Ny+1)*2+2;LL(ik)=d2cB*wd2x(4);ik=ik+1; % ux8
                    end
                end
            else
                I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
            end
            if (ix<Nx+1)
                % ux lives at (x, yp); same split as above.
                dx=m.ux.hxc(ix); dy=m.ux.hyc(iy);
                hxm=m.ux.hxm(ix); hxp=m.ux.hxp(ix);
                hym=m.ux.hym(iy); hyp=m.ux.hyp(iy);
                if (iy==1 && ix==(Nx+1)/2)
                    % Fault-line ghost at (x=0, y=-dy/2): on the fault plane,
                    % half a cell above the free surface. sigma_xz there IS the
                    % fault shear traction, so imposing sigma_xz=0 destroyed the
                    % trace node once a rupture reached the surface. Zero
                    % curvature along the fault instead -- constrains no
                    % derivative, so no traction at the singular corner. See the
                    % same branch in $HOME/BP3/source/build_LH.m.
                    %
                    % Non-uniform safe: a pure difference of three collinear
                    % values, no spacing enters, so it is identical on a
                    % stretched mesh. (The y core always starts at the surface,
                    % so these three nodes are uniformly spaced anyway.)
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=-2;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+4;LL(ik)=1;ik=ik+1;
                elseif (iy==1) % free surface: sigma_xz=0
                    % (dx/G)*sigma_xy=0 at x^2=0:
                    % d2(u1)+(1-2*cos(alpha)^2)d1(u2)
                    % +cos(alpha)[d2(u2)-d1(u1)]=0.
                    % d2(u1) and d1(u1) use the ux values straddling the
                    % surface; d2(u2) is one-sided at the surface.
                    G_surface=mean(G(1,max(1,min(Nx,[ix-1,ix]))),'omitnan');
                    shear_scale=dx/sina;
                    % d2(u1) straddles the surface, so its interval is the ux
                    % node spacing there, not the centred hyc. The cosa*d2(u2)
                    % term is one-sided on the uy grid and split over the two
                    % columns, hence cosa/2 times the one-sided weights --
                    % which is -3cosa/(4dy), cosa/dy, -cosa/(4dy) when uniform.
                    cw1=cosa/2*wy1; cw2=cosa/2*wy2; cw3=cosa/2*wy3;
                    I(ik)=kux;J(ik)=kux;LL(ik)=-shear_scale/hux1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=shear_scale/hux1;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy;LL(ik)=shear_scale*(cw1-(1-2*cosa*cosa)/dx);ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+2;LL(ik)=shear_scale*cw2;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+4;LL(ik)=shear_scale*cw3;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=shear_scale*(cw1+(1-2*cosa*cosa)/dx);ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=shear_scale*cw2;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2+4;LL(ik)=shear_scale*cw3;ik=ik+1;
                    if ix==1
                        I(ik)=kux;J(ik)=kux;LL(ik)=shear_scale*cosa/2/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux+2;LL(ik)=shear_scale*cosa/2/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=-shear_scale*cosa/2/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-shear_scale*cosa/2/dx;ik=ik+1;
                    elseif ix==Nx
                        I(ik)=kux;J(ik)=kux;LL(ik)=-shear_scale*cosa/2/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux+2;LL(ik)=-shear_scale*cosa/2/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=shear_scale*cosa/2/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux-(Ny+1)*2+2;LL(ik)=shear_scale*cosa/2/dx;ik=ik+1;
                    else
                        I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=-shear_scale*cosa/4/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-shear_scale*cosa/4/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=shear_scale*cosa/4/dx;ik=ik+1;
                        I(ik)=kux;J(ik)=kux-(Ny+1)*2+2;LL(ik)=shear_scale*cosa/4/dx;ik=ik+1;
                    end
                elseif (iy==Ny+1) % far-bottom boundary ux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=1;ik=ik+1;
                elseif (ix==1) % left boundary ux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                elseif (ix==Nx) % right boundary ux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2) % fault continous sigma?
                    GA=G(iy,ix-1);GB=G(iy,ix+1);GC=(GA+GB)/2;
                    lambdaA=lambda(iy,ix-1);lambdaB=lambda(iy,ix+1);lambdaC=(lambdaA+lambdaB)/2;
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=-(lambda+2*G)/G;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
                    % Normal-traction continuity: each side's d(ux)/dx over its
                    % own interval, the row rescaled by dx=hxc so that on a
                    % uniform grid dx/hxm=dx/hxp=1 and these are the previous
                    % weights exactly.
                    I(ik)=kux;J(ik)=kux;LL(ik)=-dx*((lambdaA+2*GA)/hxm ...
                        +(lambdaB+2*GB)/hxp)/GC;ik=ik+1; % ux3
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=dx*(lambdaB+2*GB)/(GC*hxp);ik=ik+1; % ux4
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=dx*(lambdaA+2*GA)/(GC*hxm);ik=ik+1; % ux2
                    % 2-interval centred d/dy, so it gets the exact weights like
                    % every other. It vanishes when GA==GB, which is every
                    % homogeneous run including BP3 -- that is not a reason to
                    % leave it wrong, since a layered model would wake it up and
                    % the error would then be silent and mesh dependent.
                    % CA is set so CA*byp reproduces the previous +K exactly.
                    CA=2*(GA-GB)*cosa/GC*dx;
                    [bym,by0,byp]=centred3(hym,hyp);
                    if m.uniform, by0=0; end
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=CA*byp;ik=ik+1; % ux5 added
                    I(ik)=kux;J(ik)=kux;  LL(ik)=CA*by0;ik=ik+1; % ux3 added (graded only)
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=CA*bym;ik=ik+1; % ux1 added
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2-2;LL(ik)=2*cosa/dy*dx/4;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-2*cosa/dy*dx/4;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2-2;LL(ik)=-2*cosa/dy*dx/4;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2+2;LL(ik)=2*cosa/dy*dx/4;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy;LL(ik)=-lambdaA/GC/dy*dx;ik=ik+1; % uy3
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=lambdaB/GC/dy*dx;ik=ik+1; % uy4
                    I(ik)=kux;J(ik)=kuy-2;LL(ik)=lambdaA/GC/dy*dx;ik=ik+1; % uy1
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-lambdaB/GC/dy*dx;ik=ik+1; % uy2
                else % ux-Navier
                    GA=G(iy-1,ix);GB=G(iy,ix);GC=(GA+GB)/2;
                    lambdaA=lambda(iy-1,ix);lambdaB=lambda(iy,ix);lambdaC=(lambdaA+lambdaB)/2;
                    % Second derivatives in asymmetric form (increment 2b);
                    % reduces to -2 / +1 / +1 and dx^2/dy^2 on a uniform grid.
                    I(ik)=kux;J(ik)=kux;LL(ik)=-(lambdaC+2*GC)/GC*dx*(1/hxm+1/hxp) ...
                        -dx*dx/(GC*dy)*(GA/hym+GB/hyp);ik=ik+1; % ux5
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambdaC+2*GC)/GC*dx/hxm;ik=ik+1; % ux4
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=(lambdaC+2*GC)/GC*dx/hxp;ik=ik+1; % ux6
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=dx*dx*GA/(GC*dy*hym);ik=ik+1; % ux2
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=dx*dx*GB/(GC*dy*hyp);ik=ik+1; % ux8
                    % INCREMENT 2d, ux block -- same exact tensor product as the
                    % uy mixed derivative above. See centred3.m.
                    [axm,ax0,axp]=centred3(hxm,hxp);
                    [aym,ay0,ayp]=centred3(hym,hyp);
                    if m.uniform, ax0=0; ay0=0; end
                    sx=hxm+hxp;
                    Mv=[cosa*dx*(lambdaA+GA+2*GC)/GC, ...
                        cosa*dx*(lambdaC+GC+2*GC)/GC, ...
                        cosa*dx*(lambdaB+GB+2*GC)/GC];
                    wy=[aym ay0 ayp]; wx=[axm ax0 axp];
                    for ey=1:3
                        for ex=1:3
                            cxy=-Mv(ey)/2*wy(ey)*sx*wx(ex);
                            if cxy~=0
                                I(ik)=kux;
                                J(ik)=kux+(ex-2)*(Ny+1)*2+(ey-2)*2;
                                LL(ik)=cxy;ik=ik+1;   % ux1/3/7/9, +2/4/5/6/8 if graded
                            end
                        end
                    end
                    if (iy==2 || iy==Ny)
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaC+GB)/GC;ik=ik+1; % uy6
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-1/dy*dx*(lambdaC+GA)/GC;ik=ik+1; % uy4
                        I(ik)=kux;J(ik)=kuy;LL(ik)=-1/dy*dx*(lambdaC+GB)/GC;ik=ik+1; % uy5
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=1/dy*dx*(lambdaC+GA)/GC;ik=ik+1; % uy3
                    else
                        % INCREMENT 2e, mirror of the uy block. The u2 coupling
                        % in the u1 equation is
                        %     (lambda+G) * ( d1 d2 - cos(alpha) d2^2 ) u2
                        % -- again two operators previously summed into one
                        % coefficient. Split, and the second derivative gets
                        % exact non-uniform weights.
                        %
                        % (i) mixed derivative d1 d2 u2: staggered, one interval
                        % each way, correct on any mesh.
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaC+GB)/GC;ik=ik+1; % uy6
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-1/dy*dx*(lambdaC+GA)/GC;ik=ik+1; % uy4
                        I(ik)=kux;J(ik)=kuy;LL(ik)=-1/dy*dx*(lambdaC+GB)/GC;ik=ik+1; % uy5
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=1/dy*dx*(lambdaC+GA)/GC;ik=ik+1; % uy3
                        % (ii) -cos(alpha)*d2^2 u2 at the half node yp(iy), from
                        % the four uy rows y(iy-2..iy+1), averaged over the two
                        % columns ix, ix+1. This is THE term the y-coarsening
                        % bump lives on: the hard-coded (1,-1,-1,1)/(2dy^2) is
                        % 25 % wrong at the nys=36 seam. Material follows the
                        % existing per-row A/B assignment.
                        wd2y=fdweights(g.yp(iy),[g.y(iy-2) g.y(iy-1) g.y(iy) g.y(iy+1)],2);
                        d2eA=-cosa*(lambdaA+GA)/GC*dx*dx/2;   % rows iy-2 and iy
                        d2eB=-cosa*(lambdaB+GB)/GC*dx*dx/2;   % rows iy-1 and iy+1
                        I(ik)=kux;J(ik)=kuy-2*2;LL(ik)=d2eA*wd2y(1);ik=ik+1; % uy1
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2*2;LL(ik)=d2eA*wd2y(1);ik=ik+1; % uy2
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=d2eB*wd2y(2);ik=ik+1; % uy3
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=d2eB*wd2y(2);ik=ik+1; % uy4
                        I(ik)=kux;J(ik)=kuy;LL(ik)=d2eA*wd2y(3);ik=ik+1; % uy5
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=d2eA*wd2y(3);ik=ik+1; % uy6
                        I(ik)=kux;J(ik)=kuy+2;LL(ik)=d2eB*wd2y(4);ik=ik+1; % uy7
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=d2eB*wd2y(4);ik=ik+1; % uy8
                    end
                end
            else
                I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
            end
        end
    end
    % A silent overflow of the preallocation would corrupt the operator, so
    % fail loudly instead. Bump ntrip above if a new term trips this.
    if ik-1>numel(I)
        error('build_LH:tripletOverflow', ...
            ['emitted %d triplets but preallocated %d (%d per unknown); ' ...
             'raise ntrip.'],ik-1,numel(I),ntrip);
    end

    LH=sparse(I(1:ik-1),J(1:ik-1),LL(1:ik-1));
end
