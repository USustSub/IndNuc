function RH=build_RH(lambda,G,sina,cosa,dPdt,Biot,Nx,Ny,N,dx,dy,y,V,z,dz,param)
% BP3-QD right-hand side: prescribed tangential velocity jump on the fault.
%
% GRID STRETCHING: nothing here depends on the mesh, and dx/dy are unused. Every
% row this function writes to is algebraic -- a prescribed velocity on a row
% whose LH coefficients are +/-1 (side and bottom boundaries, the fault jump) --
% so no spacing enters. That remains true only while those rows keep unit
% scaling; if a future change scales them by a local dx, the matching factor has
% to be applied here too.

RH=zeros(N,1);
fault_ix=(Nx+1)/2;

% Symmetric far-field plate loading. In the oblique basis, uy is the
% fault-parallel velocity component. The LH boundary row averages the ghost
% and interior values, so its right-hand side is twice the face velocity.
if isfield(param,'load_side_boundaries') && param.load_side_boundaries
    left_velocity=-param.Vp/2;
    right_velocity=param.Vp/2;
    for iy=1:Ny
        kuy_left=(iy-1)*2+2;
        kuy_right=(Nx*(Ny+1)+iy-1)*2+2;
        RH(kuy_left)=2*left_velocity;
        RH(kuy_right)=2*right_velocity;
    end
end

% For a compact total-velocity domain, continue the rigid velocities along
% the two halves of the bottom boundary. Otherwise the deep-creep jump is
% concentrated at the single fault/bottom node while its neighbors are
% pinned to zero, creating a large artificial stress concentration.
if isfield(param,'load_bottom_boundaries') && param.load_bottom_boundaries
    bottom_left=-param.VL/2;
    bottom_right=param.VL/2;
    for ix=2:Nx
        if ix==fault_ix
            continue;
        end
        kuy=((ix-1)*(Ny+1)+Ny-1)*2+2;
        if ix<fault_ix
            RH(kuy)=bottom_left;
        else
            RH(kuy)=bottom_right;
        end
    end
end

% iy starts at 1: the surface fault node carries the jump row too, matching the
% build_LH branch that excludes both fault columns from the sigma_zz=0 case.
% The bottom fault intersection (iy==Ny) is the deep-creep driver.
for iy=1:Ny
    kuy=((fault_ix-1)*(Ny+1)+iy-1)*2+2;
    RH(kuy)=V(iy);
end
end
