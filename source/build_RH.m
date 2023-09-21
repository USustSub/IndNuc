function RH=build_RH(lambda,G,sina,cosa,dPdt,Nx,Ny,N,dx,dy,y,V)
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
                else % make it more accurate
%                     RH(kuy)=-rho*g*sina*dx*dx/G*0;
                    if (y(iy)==850&&ix>=(Nx+1)/2+1)
                        RH(kuy)=dPdt/dy*dx*dx/G*sina;
                    end
                    if (y(iy)==1050&&ix>=(Nx+1)/2+1)
                        RH(kuy)=-dPdt/dy*dx*dx/G*sina;
                    end                  
                    if (y(iy)==800&&ix<=(Nx+1)/2)
                        RH(kuy)=dPdt/dy*dx*dx/G*sina;
                    end
                    if (y(iy)==1000&&ix<=(Nx+1)/2)
                        RH(kuy)=-dPdt/dy*dx*dx/G*sina;
                    end
                end
            end
            if (ix<Nx+1)
                if (iy==1)
                elseif (iy==Ny+1)
                elseif (ix==1)
                elseif (ix==Nx)
                elseif (ix==(Nx+1)/2)
                    if (y(iy)>800&&y(iy)<=850)
                        RH(kux)=-dPdt*dx/G;
                    end
                    if (y(iy)>1000&&y(iy)<=1050)
                        RH(kux)=dPdt*dx/G;
                    end
                else
%                     RH(kux)=rho*g*sina*cosa*dx*dx/G*0;
                    if (y(iy)==1050&&ix>(Nx+1)/2+1)
                        RH(kux)=dPdt/dy*dx*dx/G*sina*cosa;
                    end
                    if (y(iy)==1000&&ix<(Nx+1)/2+1)
                        RH(kux)=dPdt/dy*dx*dx/G*sina*cosa;
                    end
                    if (y(iy)==850&&ix>(Nx+1)/2+1)
                        RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa;
                    end
                    if (y(iy)==800&&ix<(Nx+1)/2+1)
                        RH(kux)=-dPdt/dy*dx*dx/G*sina*cosa;
                    end
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