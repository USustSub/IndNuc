function plot_screenshot(vx,tauqs,vy,sigmaqs)
load('coord.mat','Xtau','Xux','Xuy','Ytau','Yux','Yuy','Xsigma','Ysigma');
figure(11)
subplot(2,2,1)
pcolor(Xux,Yux,vx);axis ij;colorbar
subplot(2,2,2)
pcolor(Xtau,Ytau,tauqs);axis ij;colorbar
subplot(2,2,3)
pcolor(Xuy,Yuy,vy);axis ij;colorbar
subplot(2,2,4)
pcolor(Xsigma,Ysigma,sigmaqs);axis ij;colorbar
end