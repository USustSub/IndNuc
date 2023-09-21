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
hold on

    subplot(2,2,1)
    xlabel('Xux, Yux, vx');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
    subplot(2,2,2)
    xlabel('Xtau,Ytau,tauqs');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
    subplot(2,2,3)
    xlabel('Xuy,Yuy,vy');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
    subplot(2,2,4)
    xlabel('Xsigma,Ysigma,sigmaqs');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
end

