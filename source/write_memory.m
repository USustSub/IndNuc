function write_memory(it,output_interval,U,V,tau,sigma,P,theta,dt,t,t2,tauqs,sigmaqs,uy,vy,ux,vx)
global Um Vm taum sigmam Pm thetam dtm tm tm2 taumall sigmamall uymall vymall uxmall vxmall
Um(:,it/output_interval)=U;
Vm(:,it/output_interval)=V;
taum(:,it/output_interval)=tau;
sigmam(:,it/output_interval)=sigma;
Pm(:,it/output_interval)=P;
thetam(:,it/output_interval)=theta;
dtm(it/output_interval)=dt;
tm(it/output_interval)=t;
tm2(it/output_interval)=t2;
% taumall(:,:,it/output_interval)=tauqs;
% sigmamall(:,:,it/output_interval)=sigmaqs;
% uymall(:,:,it/output_interval)=uy;
% vymall(:,:,it/output_interval)=vy;
% uxmall(:,:,it/output_interval)=ux;
% vxmall(:,:,it/output_interval)=vx;
end
