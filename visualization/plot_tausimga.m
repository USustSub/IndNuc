function plot_tausimga(sigmam, taum, Nt, output_interval, Ny, y) 
% PLOT the shear and the effective normal stress together
%   This shows the temporal evolution of the ratio shear stress and
%   effective normal stress. 

%to plot the colomns over a convenient interval.

st = 100;  %stepsize of the plot, can adjust self!

colomncount  = (Nt/output_interval);
colomnIndices = 1:st:colomncount;  

%make an array of the ratio with the .
ratio = taum./sigmam;

figure; %start with the figure! 
plot(y + 2000, ratio(: , colomnIndices));  

%plot statements
xlabel ('Depth [m]');
ylabel ('\tau/\sigma_n');
title('Ratio shear and normal eff stress over the depth')
set(gca,'FontSize',18);






%colomncount  = (Nt/output_interval);
%colomnIndices = 1:step:colomncount;

%plot the depth over the Ny value
% depth_start = 2000;
% depth_eind = 4000;
% depth_array = linspace(depth_start,depth_eind,Ny);
% depth_start:(depth_eind - depth_start)/(Ny-1):depth_eind;
 
%make an array of the ratio
%ratio = taum./sigmam;

%figure; %start with the figure! 

%for loop over the lines that will be plotted
% for i = colomnIndices
%      hold on;                   %to plot the diffferent lines in one figure
%     plot(y+2000, ratio(:,colomnIndices));    
% end

%why does it plot 1500 t/m -1500 when the ratio does not exceed 1?

%plot statements
% xlabel ('Depth [m]');
% ylabel ('\tau/\sigma');
% title('Ratio shear and normal eff stress over the depth')
% set(gca,'FontSize',18);


end

