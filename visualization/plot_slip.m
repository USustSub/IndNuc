function plot_slip(Vm, Nt, output_interval, Ny, y) 
%plot slip rate and the tau/sigma 
%   This function can be used to plot the slip rates over depths for different 
%running times of the model. Thus it shows the temporal evolution of the
%slip rate. 

%Slip: Vm is an (Nx:Nt/output_interval) array. 

%define the indiced of the coloms to plot from Vm 
%colomnIndices = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];

st = 100;   %step interval to plot
colomncount  = (Nt/output_interval);
colomnIndices = 1:st:colomncount;

%plot the selected colomns as seperate lines
figure;

%plot
plot(y+2000, Vm(:,colomnIndices));   

% finalize the plot
xlabel ('Depth (m)');
ylabel ('Slip rate (V)');
title('Slip rate over the depth ')
set(gca,'FontSize',18);

end

