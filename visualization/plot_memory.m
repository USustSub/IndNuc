function plot_memory(yr,colordot)
global tm Vm taum thetam sigmam Pm
figure(13)
subplot(2,2,1)
semilogy(tm/yr,Vm,colordot);
hold on
subplot(2,2,2)
plot(tm/yr,taum./sigmam,colordot);
hold on
subplot(2,2,3)
semilogy(tm/yr,thetam,colordot);
hold on
subplot(2,2,4)
plot(tm/yr,Pm/1e6,colordot);
hold on
    subplot(2,2,1)
    xlabel('time (yr)');
    ylabel('slip rate V (m/s)');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
    subplot(2,2,2)
    xlabel('time (yr)');
    ylabel('shear stress \tau (MPa)');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
    subplot(2,2,3)
    xlabel('time (yr)');
    ylabel('state \theta (s)');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
    subplot(2,2,4)
    xlabel('time (yr)');
    ylabel('P (MPa)');
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',18);
%     xlim([0 4000])
end
