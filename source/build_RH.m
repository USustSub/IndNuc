function RH=build_RH(lambda,G,sina,cosa,dPdt,Biot,Nx,Ny,N,dx,dy,y,V,z,dz,param)
    TBt=param.BZb-param.ytop;
    SSt=param.TBb-param.ytop;
    SSb=param.SSb-param.ytop;
    offset=param.offset;
    sina=1; %%

    dPdt.L=dPdt.L*Biot;
    dPdt.R=dPdt.R*Biot;
    % d2Pdt.L=[0;diff(movmean(dPdt.L,2,"Endpoints","discard"));0];
    % d2Pdt.R=[0;diff(movmean(dPdt.R,2,"Endpoints","discard"));0];
    d2Pdt.L=[0;diff(dPdt.L)];
    d2Pdt.R=[0;diff(dPdt.R)];

    RH=zeros(N,1);
    for ix=1:Nx+1
        for iy=1:Ny+1
            kux=((ix-1)*(Ny+1)+iy-1)*2+1;
            kuy=kux+1;
            if (iy<Ny+1)
                if (ix==1)
                elseif (ix==Nx+1)
                elseif (iy==1)
                elseif (iy==Ny)
                elseif (ix==(Nx+1)/2)
                    RH(kuy)=V(iy);
                elseif (ix==(Nx+1)/2+1)
%                 elseif (iy==1)
%                     RH(kuy)=0;
%                 elseif (iy==Ny)
%                     RH(kuy)=0;
%                     if (ix<=(Nx+1)/2)
%                         RH(kuy)=-Vp;
%                     else
%                         RH(kuy)=Vp;
%                     end
                else
                    GA=G(iy,ix);GB=G(iy+1,ix);GC=(GA+GB)/2;
%                     RH(kuy)=-rho*g*sina*dx*dx/G*0;
%                     if (z(iy)+dz/2>SSt&&z(iy)-dz/2<=SSt&&ix>=(Nx+1)/2+1)%
%                         RH(kuy)=dPdt/dy*dx*dx/GC*sina;
%                     end
%                     if (z(iy)+dz/2>SSb&&z(iy)-dz/2<=SSb&&ix>=(Nx+1)/2+1)
%                         RH(kuy)=-dPdt/dy*dx*dx/GC*sina;
%                     end
%                     if (z(iy)+dz/2>SSt-offset&&z(iy)-dz/2<=SSt-offset&&ix<=(Nx+1)/2)
%                         RH(kuy)=dPdt/dy*dx*dx/GC*sina;
%                     end
%                     if (z(iy)+dz/2>SSb-offset&&z(iy)-dz/2<=SSb-offset&&ix<=(Nx+1)/2)
%                         RH(kuy)=-dPdt/dy*dx*dx/GC*sina;
%                     end
                    % if (z(iy)>TBt&&z(iy-1)<=TBt&&ix>=(Nx+1)/2+1)%
                    %     RH(kuy)=dPdtTB/dy*dx*dx/GC*sina;
                    % end
                    % if (z(iy)>SSt&&z(iy-1)<=SSt&&ix>=(Nx+1)/2+1)%
                    %     RH(kuy)=(dPdt-dPdtTB)/dy*dx*dx/GC*sina;
                    % end
                    % if (z(iy+1)>SSb&&z(iy)<=SSb&&ix>=(Nx+1)/2+1)
                    %     RH(kuy)=-dPdt/dy*dx*dx/GC*sina;
                    % end

                    if (ix>=(Nx+1)/2+1)
                        RH(kuy)=(d2Pdt.R(iy)+d2Pdt.R(iy+1))/2/dy*dx*dx/GC;
                    end
                    if (ix<=(Nx+1)/2)
                        RH(kuy)=(d2Pdt.L(iy)+d2Pdt.L(iy+1))/2/dy*dx*dx/GC;
                    end

                    % if (z(iy)>TBt-offset&&z(iy-1)<=TBt-offset&&ix<=(Nx+1)/2)
                    %     RH(kuy)=dPdtTB/dy*dx*dx/GC*sina;
                    % end
                    % if (z(iy)>SSt-offset&&z(iy-1)<=SSt-offset&&ix<=(Nx+1)/2)
                    %     RH(kuy)=(dPdt-dPdtTB)/dy*dx*dx/GC*sina;
                    % end
                    % if (z(iy+1)>SSb-offset&&z(iy)<=SSb-offset&&ix<=(Nx+1)/2)
                    %     RH(kuy)=-dPdt/dy*dx*dx/GC*sina;
                    % end
                end
            end
            if (ix<Nx+1)
                if (iy==1)
                elseif (iy==Ny+1)
                elseif (ix==1)
                elseif (ix==Nx)
                elseif (ix==(Nx+1)/2)
                    GA=G(iy,ix-1);GB=G(iy,ix+1);GC=(GA+GB)/2;
%                     if (z(iy-1)>SSt-offset&&z(iy-1)<=SSt)%
%                         RH(kux)=-dPdt*dx/GC;
%                     end
%                     if (z(iy)>SSb-offset&&z(iy)<=SSb)
%                         RH(kux)=dPdt*dx/GC;
%                     end
                    % RH(kux)=0;
                    % if (z(iy-1)>TBt-offset&&z(iy-1)<=SSt-offset)%
                    %     RH(kux)=RH(kux)-dPdtTB*dx/GC;
                    % end
                    % if (z(iy-1)>SSt-offset&&z(iy-1)<=SSb-offset)%
                    %     RH(kux)=RH(kux)-dPdt*dx/GC;
                    % end
                    % if (z(iy)>TBt&&z(iy)<=SSt)
                    %     RH(kux)=RH(kux)+dPdtTB*dx/GC;
                    % end
                    % if (z(iy)>SSt&&z(iy)<=SSb)
                    %     RH(kux)=RH(kux)+dPdt*dx/GC;
                    % end

                    RH(kux)=(dPdt.R(iy)-dPdt.L(iy)+dPdt.R(iy-1)-dPdt.L(iy-1))/2*dx/GC;
                else
                    GA=G(iy-1,ix);GB=G(iy,ix);GC=(GA+GB)/2;
%                     RH(kux)=rho*g*sina*cosa*dx*dx/G*0;
%                     if (z(iy)+dz/2>SSb&&z(iy)-dz/2<=SSb&&ix>(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/GC*sina*cosa;
%                     end
%                     if (z(iy)+dz/2>SSb-offset&&z(iy)-dz/2<=SSb-offset&&ix<(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/GC*sina*cosa;
%                     end
%                     if (z(iy)+dz/2>SSt&&z(iy)-dz/2<=SSt&&ix>(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/GC*sina*cosa;
%                     end
%                     if (z(iy)+dz/2>SSt-offset&&z(iy)-dz/2<=SSt-offset&&ix<(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/GC*sina*cosa;
%                     end
                    % if (z(iy)>TBt&&z(iy-1)<=TBt&&ix>(Nx+1)/2+1)%okay
                    %     RH(kux)=-dPdtTB/dy*dx*dx/GC*sina*cosa;
                    % end
                    % if (z(iy)>SSt&&z(iy-1)<=SSt&&ix>(Nx+1)/2+1)%okay
                    %     RH(kux)=(-dPdt+dPdtTB)/dy*dx*dx/GC*sina*cosa;
                    % end
                    % if (z(iy)>SSb&&z(iy-1)<=SSb&&ix>(Nx+1)/2+1)
                    %     RH(kux)=dPdt/dy*dx*dx/GC*sina*cosa;
                    % end

                    if (ix>(Nx+1)/2+1)
                        RH(kux)=-d2Pdt.R(iy)/dy*dx*dx/GC*cosa;
                    end
                    if (ix<(Nx+1)/2+1)
                        RH(kux)=-d2Pdt.L(iy)/dy*dx*dx/GC*cosa;
                    end

                    % if (z(iy)>TBt-offset&&z(iy-1)<=TBt-offset&&ix<(Nx+1)/2+1)%okay
                    %     RH(kux)=-dPdtTB/dy*dx*dx/GC*sina*cosa;
                    % end
                    % if (z(iy)>SSt-offset&&z(iy-1)<=SSt-offset&&ix<(Nx+1)/2+1)%okay
                    %     RH(kux)=(-dPdt+dPdtTB)/dy*dx*dx/GC*sina*cosa;
                    % end
                    % if (z(iy)>SSb-offset&&z(iy-1)<=SSb-offset&&ix<(Nx+1)/2+1)
                    %     RH(kux)=dPdt/dy*dx*dx/GC*sina*cosa;
                    % end
                end
            end
        end
    end

end

% case 1
%                     if (y(iy)==1050&&ix>(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa;
%                     end
%                     if (y(iy)==1000&&ix<(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa;
%                     end
%                     if (y(iy-1)==850&&ix>(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa;
%                     end
%                     if (y(iy-1)==800&&ix<(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa;
%                     end
% case 2
%                     if (y(iy)==850&&ix>(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy)==1050&&ix>(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy)==800&&ix<(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy)==1000&&ix<(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy-1)==850&&ix>(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy-1)==1050&&ix>(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy-1)==800&&ix<(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
%                     if (y(iy-1)==1000&&ix<(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa/2;
%                     end
% case 3
%                     if (y(iy)==1050&&ix>(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa;
%                     end
%                     if (y(iy)==1000&&ix<(Nx+1)/2+1)
%                         RH(kux)=dPdt/dy*dx*dx/G*sina*cosa;
%                     end
%                     if (y(iy)==850&&ix>(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa;
%                     end
%                     if (y(iy)==800&&ix<(Nx+1)/2+1)
%                         RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa;
%                     end
