%% PlotFinalRun.m
% By Jack Hannum, University of South Carolina
% Plot results of a VFD driving the Emrax 228
close all; clc;
name = 'CapDegradation';
data = readtable('cap.csv');
time_seconds = data.Column1;
temperature_celsius = data.Column2;
lifetime_hours = data.Column3;


%% Plot Capacitor Degradation
CapDegradation = figure('name','Capacitor Degradation','units','inches');
CapDegradation.Position = [1 1 3 3];
plot(time_seconds,lifetime_hours,'LineWidth',2,'DisplayName','VSI Capacitor');
grid on;
xlabel('Time (s)');
ylabel('Lifetime (hours)');
legend('location','northeast','NumColumns',3);
xlim([0 60]);
savefig(CapDegradation,'CapDegradation.fig');
exportgraphics(CapDegradation,strcat(name,'.pdf'),'ContentType','vector');