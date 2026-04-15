function plot_movie(Nt,y)
global Vm taum
figure(4)
for it=1:1:Nt
    semilogy(y,Vm(:,it));
    axis([-inf inf 1e-18 1e2])
    pause(0.0001)
end

figure(5)
for it=1:10:Nt
    plot(y,taum(:,it));
    hold on
%     axis([-inf inf 1e-18 1e2])
    pause(0.0001)
end
end