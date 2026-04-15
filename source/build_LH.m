function LH=build_LH(lambda,G,sina,cosa,Nx,Ny,N,dx,dy)
% G=zeros(Ny+1,Nx+1)+G0;
% lambda=zeros(Ny+1,Nx+1)+lambda0;
    I=zeros(17*N,1);
    J=zeros(17*N,1);
    LL=zeros(17*N,1);
    ik=1;
    for ix=1:Nx+1
        for iy=1:Ny+1
            kux=((ix-1)*(Ny+1)+iy-1)*2+1;
            kuy=kux+1;
            if (iy<Ny+1)
                if (ix==1) % left boundary duy=0
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=-1;ik=ik+1;%
                elseif (ix==Nx+1) % right boundary duy=0
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=-1;ik=ik+1;%
                elseif (iy==1) % top boundary duy=0
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-1;ik=ik+1;%
                elseif (iy==Ny) % bottom boundary uy=0
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-2;LL(ik)=1;ik=ik+1;%
                elseif (ix==(Nx+1)/2) % fault diff(vy)=Vy
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=-2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2*(Ny+1)*2;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2+1) % fault continous tau?
                    GA=G(iy,ix-2);GB=G(iy,ix);GC=(GA+GB)/2;
                    I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2;LL(ik)=1*GA/GC;ik=ik+1; % uy3
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=-1*GA/GC;ik=ik+1; % uy4
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-1*GB/GC;ik=ik+1; % uy5
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1*GB/GC;ik=ik+1; % uy6
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cosa/2/dy*dx*GA/GC;ik=ik+1; % uy1
                    I(ik)=kuy;J(ik)=kuy-2;LL(ik)=cosa/2/dy*dx*GB/GC;ik=ik+1; % uy2
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cosa/2/dy*dx*GA/GC;ik=ik+1; % uy7
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-cosa/2/dy*dx*GB/GC;ik=ik+1; % uy8
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
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-2-1/dy/dy*dx*dx*(lambdaA+2*GA+lambdaB+2*GB)/GC;ik=ik+1; % uy5
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=1;ik=ik+1; % uy4
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1; % uy6
                    I(ik)=kuy;J(ik)=kuy-2;LL(ik)=1/dy/dy*dx*dx*(lambdaA+2*GA)/GC;ik=ik+1; % uy2
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=1/dy/dy*dx*dx*(lambdaB+2*GB)/GC;ik=ik+1; % uy8
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=cosa/dy*dx*(lambdaC+GC+2*GA)/GC/4;ik=ik+1; % uy3
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=-cosa/dy*dx*(lambdaC+GC+2*GB)/GC/4;ik=ik+1; % uy9
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cosa/dy*dx*(lambdaC+GC+2*GA)/GC/4;ik=ik+1; % uy1
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cosa/dy*dx*(lambdaC+GC+2*GB)/GC/4;ik=ik+1; % uy7
                    if (ix==2 || ix==Nx)
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaA+GC)/GC;ik=ik+1; % ux2
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-1/dy*dx*(lambdaB+GC)/GC;ik=ik+1; % ux6
                        I(ik)=kuy;J(ik)=kux;LL(ik)=-1/dy*dx*(lambdaA+GC)/GC;ik=ik+1; % ux3
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=1/dy*dx*(lambdaB+GC)/GC;ik=ik+1; % ux7
                    else
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaA+GC)/GC+cosa*(lambdaA+GA)/GC/4;ik=ik+1; % ux2
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-1/dy*dx*(lambdaB+GC)/GC+cosa*(lambdaB+GB)/GC/4;ik=ik+1; % ux6
                        I(ik)=kuy;J(ik)=kux;LL(ik)=-1/dy*dx*(lambdaA+GC)/GC+cosa*(lambdaA+GA)/GC/4;ik=ik+1; % ux3
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=1/dy*dx*(lambdaB+GC)/GC+cosa*(lambdaB+GB)/GC/4;ik=ik+1; % ux7
                        I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2;LL(ik)=-cosa*(lambdaA+GA)/GC/4;ik=ik+1; % ux1
                        I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2+2;LL(ik)=-cosa*(lambdaB+GB)/GC/4;ik=ik+1; % ux5
                        I(ik)=kuy;J(ik)=kux+(Ny+1)*2;LL(ik)=-cosa*(lambdaA+GA)/GC/4;ik=ik+1; % ux4
                        I(ik)=kuy;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-cosa*(lambdaB+GB)/GC/4;ik=ik+1; % ux8
                    end
                end
            else
                I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
            end
            if (ix<Nx+1)
                if (iy==1) % top boundary dux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=-1;ik=ik+1;%
                elseif (iy==Ny+1) % bottom boundary dux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=-1;ik=ik+1;%
                elseif (ix==1) % left boundary ux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                elseif (ix==Nx) % right boundary ux=0
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2) % fault continous sigma?
                    GA=G(iy,ix-1);GB=G(iy,ix+1);GC=(GA+GB)/2;
                    lambdaA=lambda(iy,ix-1);lambdaB=lambda(iy,ix+1);lambdaC=(lambdaA+lambdaB)/2;
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=-(lambda+2*G)/G;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kux;J(ik)=kux;LL(ik)=-(lambdaA+2*GA+lambdaB+2*GB)/GC;ik=ik+1; % ux3
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=(lambdaB+2*GB)/GC;ik=ik+1; % ux4
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambdaA+2*GA)/GC;ik=ik+1; % ux2
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=2*(GA-GB)*cosa/GC/2/dy*dx;ik=ik+1; % ux5 added
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=-2*(GA-GB)*cosa/GC/2/dy*dx;ik=ik+1; % ux1 added
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
                    I(ik)=kux;J(ik)=kux;LL(ik)=-2*(lambdaC+2*GC)/GC-1/dy/dy*dx*dx*(GA+GB)/GC;ik=ik+1; % ux5
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambdaC+2*GC)/GC;ik=ik+1; % ux4
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=(lambdaC+2*GC)/GC;ik=ik+1; % ux6
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=1/dy/dy*dx*dx*GA/GC;ik=ik+1; % ux2
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=1/dy/dy*dx*dx*GB/GC;ik=ik+1; % ux8
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2-2;LL(ik)=cosa/dy*dx*(lambdaA+GA+2*GC)/GC/4;ik=ik+1; % ux3
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-cosa/dy*dx*(lambdaB+GB+2*GC)/GC/4;ik=ik+1; % ux9
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2-2;LL(ik)=-cosa/dy*dx*(lambdaA+GA+2*GC)/GC/4;ik=ik+1; % ux1
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2+2;LL(ik)=cosa/dy*dx*(lambdaB+GB+2*GC)/GC/4;ik=ik+1; % ux7
                    if (iy==2 || iy==Ny)
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaC+GB)/GC;ik=ik+1; % uy6
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-1/dy*dx*(lambdaC+GA)/GC;ik=ik+1; % uy4
                        I(ik)=kux;J(ik)=kuy;LL(ik)=-1/dy*dx*(lambdaC+GB)/GC;ik=ik+1; % uy5
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=1/dy*dx*(lambdaC+GA)/GC;ik=ik+1; % uy3
                    else
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=1/dy*dx*(lambdaC+GB)/GC+cosa/dy/dy*dx*dx*(lambdaA+GA)/GC/4;ik=ik+1; % uy6
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-1/dy*dx*(lambdaC+GA)/GC+cosa/dy/dy*dx*dx*(lambdaB+GB)/GC/4;ik=ik+1; % uy4
                        I(ik)=kux;J(ik)=kuy;LL(ik)=-1/dy*dx*(lambdaC+GB)/GC+cosa/dy/dy*dx*dx*(lambdaA+GA)/GC/4;ik=ik+1; % uy5
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=1/dy*dx*(lambdaC+GA)/GC+cosa/dy/dy*dx*dx*(lambdaB+GB)/GC/4;ik=ik+1; % uy3
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=-cosa/dy/dy*dx*dx*(lambdaB+GB)/GC/4;ik=ik+1; % uy8
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2*2;LL(ik)=-cosa/dy/dy*dx*dx*(lambdaA+GA)/GC/4;ik=ik+1; % uy2
                        I(ik)=kux;J(ik)=kuy+2;LL(ik)=-cosa/dy/dy*dx*dx*(lambdaB+GB)/GC/4;ik=ik+1; % uy7
                        I(ik)=kux;J(ik)=kuy-2*2;LL(ik)=-cosa/dy/dy*dx*dx*(lambdaA+GA)/GC/4;ik=ik+1; % uy1
                    end
                end
            else
                I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
            end
        end
    end
    LH=sparse(I(1:ik-1),J(1:ik-1),LL(1:ik-1));
end
