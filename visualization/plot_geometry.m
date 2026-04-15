function plot_geometry(rho,G,lambda)
% load('coord.mat','Xtau','Ytau');
figure(11)
subplot(2,2,1)
pcolor(Xtau,Ytau,rho);axis ij;colorbar;shading interp
title('\rho [kg/m^3]')
subplot(2,2,2)
pcolor(Xtau,Ytau,G/1e9);axis ij;colorbar;shading interp
title('G [GPa]')
subplot(2,2,4)
pcolor(Xtau,Ytau,lambda/1e9);axis ij;colorbar;shading interp
title('\lambda [GPa]')
subplot(2,2,3)
pcolor(Xtau,Ytau,lambda/(lambda+G)/2);axis ij;colorbar;shading interp
title('\nv [-]')
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperPosition', [-1.5 0 29.7+3 11+0]);
set(gcf, 'PaperSize', [29.7 11]);
% saveas(gcf, 'nov8geo.pdf');
end