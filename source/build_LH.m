function LH=build_LH(lambda,G,sina,cosa,Nx,Ny,N,dx,dy)
    I=zeros(17*N,1);
    J=zeros(17*N,1);
    LL=zeros(17*N,1);
    ik=1;
    for ix=1:Nx+1
        for iy=1:Ny+1
            kux=((ix-1)*(Ny+1)+iy-1)*2+1;
            kuy=kux+1;
            if (iy<Ny+1)
                if (ix==1)
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=-1;ik=ik+1;
                elseif (ix==Nx+1)
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=-1;ik=ik+1;
                elseif (iy==1)
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2;LL(ik)=1;ik=ik+1;
                elseif (iy==Ny)
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-2;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2)
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=-2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=2;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2*(Ny+1)*2;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2+1)
                    I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
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
                    I(ik)=kuy;J(ik)=kux+(Ny+1)*2;LL(ik)=cosa/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kux+(Ny+1)*2+2;LL(ik)=cosa/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=-cosa/2;ik=ik+1;
                    I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-cosa/2;ik=ik+1;
                    I(ik)=kuy;J(ik)=kux-3*(Ny+1)*2;LL(ik)=cosa/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kux-3*(Ny+1)*2+2;LL(ik)=cosa/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2-2;LL(ik)=-cosa/4/dy*dx;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-2*(Ny+1)*2+2;LL(ik)=cosa/4/dy*dx;ik=ik+1;
%                 elseif (iy==1)
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy+2;LL(ik)=-1;ik=ik+1;
%                 elseif (iy==Ny)
%                     I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
%                     I(ik)=kuy;J(ik)=kuy-2;LL(ik)=-1;ik=ik+1;
                else
                    I(ik)=kuy;J(ik)=kuy;LL(ik)=-2-2/dy/dy*dx*dx*(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2;LL(ik)=1;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-2;LL(ik)=1/dy/dy*dx*dx*(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+2;LL(ik)=1/dy/dy*dx*dx*(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=-cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2-2;LL(ik)=-cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    I(ik)=kuy;J(ik)=kuy-(Ny+1)*2+2;LL(ik)=cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    if (ix==2 || ix==Nx)
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=1/dy*dx*(lambda+G)/G;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-1/dy*dx*(lambda+G)/G;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux;LL(ik)=-1/dy*dx*(lambda+G)/G;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=1/dy*dx*(lambda+G)/G;ik=ik+1;
                    else
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2;LL(ik)=1/dy*dx*(lambda+G)/G+cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux-(Ny+1)*2+2;LL(ik)=-1/dy*dx*(lambda+G)/G+cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux;LL(ik)=-1/dy*dx*(lambda+G)/G+cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux+2;LL(ik)=1/dy*dx*(lambda+G)/G+cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2;LL(ik)=-cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux-2*(Ny+1)*2+2;LL(ik)=-cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux+(Ny+1)*2;LL(ik)=-cosa*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kuy;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-cosa*(lambda+G)/G/4;ik=ik+1;
                    end
                end
            else
                I(ik)=kuy;J(ik)=kuy;LL(ik)=1;ik=ik+1;
            end
            if (ix<Nx+1)
                if (iy==1)
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=-1;ik=ik+1;
                elseif (iy==Ny+1)
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=-1;ik=ik+1;
                elseif (ix==1)
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                elseif (ix==Nx)
                    I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
                elseif (ix==(Nx+1)/2)
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=-(lambda+2*G)/G;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kux;J(ik)=kux;LL(ik)=-2*(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2-2;LL(ik)=2*cosa/dy*dx/4;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-2*cosa/dy*dx/4;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2-2;LL(ik)=-2*cosa/dy*dx/4;ik=ik+1;
%                     I(ik)=kux;J(ik)=kux-(Ny+1)*2+2;LL(ik)=2*cosa/dy*dx/4;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy;LL(ik)=-lambda/G/dy*dx;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=lambda/G/dy*dx;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy-2;LL(ik)=lambda/G/dy*dx;ik=ik+1;
                    I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-lambda/G/dy*dx;ik=ik+1;
                else
                    I(ik)=kux;J(ik)=kux;LL(ik)=-2*(lambda+2*G)/G-2/dy/dy*dx*dx;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2;LL(ik)=(lambda+2*G)/G;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-2;LL(ik)=1/dy/dy*dx*dx;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+2;LL(ik)=1/dy/dy*dx*dx;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2-2;LL(ik)=cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    I(ik)=kux;J(ik)=kux+(Ny+1)*2+2;LL(ik)=-cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2-2;LL(ik)=-cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    I(ik)=kux;J(ik)=kux-(Ny+1)*2+2;LL(ik)=cosa/dy*dx*(lambda+3*G)/G/4;ik=ik+1;
                    if (iy==2 || iy==Ny)
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=1/dy*dx*(lambda+G)/G;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-1/dy*dx*(lambda+G)/G;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy;LL(ik)=-1/dy*dx*(lambda+G)/G;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=1/dy*dx*(lambda+G)/G;ik=ik+1;
                    else
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2;LL(ik)=1/dy*dx*(lambda+G)/G+cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2;LL(ik)=-1/dy*dx*(lambda+G)/G+cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy;LL(ik)=-1/dy*dx*(lambda+G)/G+cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy-2;LL(ik)=1/dy*dx*(lambda+G)/G+cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2+2;LL(ik)=-cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy+(Ny+1)*2-2*2;LL(ik)=-cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy+2;LL(ik)=-cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                        I(ik)=kux;J(ik)=kuy-2*2;LL(ik)=-cosa/dy/dy*dx*dx*(lambda+G)/G/4;ik=ik+1;
                    end
                end
            else
                I(ik)=kux;J(ik)=kux;LL(ik)=1;ik=ik+1;
            end
        end
    end
    LH=sparse(I(1:ik-1),J(1:ik-1),LL(1:ik-1));
end