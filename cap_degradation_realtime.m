
% Given parameters % analytical values
clc;
clear all;
L0 = 10000;  % Expected lifetime for full rated voltage and temperature
V0 = 400;    % Updated rated voltage
T0 = 27;     % Rated temperature
Vx = 400;    % V_X is 400V

% Create an array of temperatures from 0 to 30°C
Tx = 0:60;

% Calculate L_x for each temperature with n = 3
n = 3;
Lx = L0 * (Vx/V0).^(-n) .* 2.^((T0 - Tx)/10);

% Plot L_x vs. T_x
plot(Tx, Lx, 'b-', 'LineWidth', 2);
xlabel('Temperature (°C)');
ylabel('Lifetime(hours)');
%title('Lifetime vs. Temperature');
grid on;
%% plots for overlapping real time and analytical
% Read data from the CSV file
data = readtable('cap.csv');

% Extract columns from the table
time_seconds = data.Column1;
temperature_celsius = data.Column2;
lifetime_hours = data.Column3;

% Plot lifetime vs. time
figure;
subplot(2, 1, 1);
plot(time_seconds, lifetime_hours, 'b-', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Lifetime (hours)');
%title('Lifetime vs. Time');
grid on;

% Plot lifetime vs. temperature
subplot(2, 1, 2);
plot(temperature_celsius, lifetime_hours, 'r-', 'LineWidth', 2);

xlabel('Temperature (°C)');
ylabel('Lifetime (hours)');
%title('Lifetime vs. Temperature');
grid on;

% Adjust subplot spacing
subplot(2, 1, 1);
subplot(2, 1, 2);

% Show the plots
sgtitle('Plots for Real time data');
%% plot for realtime
% Given parameters
L0 = 10000;  % Expected lifetime for full rated voltage and temperature
V0 = 400;    % Updated rated voltage
T0 = 27;     % Rated temperature
Vx = 400;    % V_X is 400V

% Create an array of temperatures from 0 to 30°C
Tx = 0:60;

% Calculate L_x for each temperature with n = 3
n = 3;
Lx = L0 * (Vx/V0).^(-n) .* 2.^((T0 - Tx)/10);

% Read data from the CSV file
data = readtable('cap.csv');

% Extract columns from the table
time_seconds = data.Column1;
temperature_celsius = data.Column2;
lifetime_hours = data.Column3;

% Plot both sets of data on the same graph
figure;
plot(Tx, Lx, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical Data');
hold on;
plot(temperature_celsius, lifetime_hours, 'r-', 'LineWidth', 2, 'DisplayName', 'Real time data');
xlabel('Temperature (°C)');
ylabel('Lifetime (hours)');
title('Lifetime vs Temperature');
grid on;
legend('Location', 'NorthEast');
hold off;
