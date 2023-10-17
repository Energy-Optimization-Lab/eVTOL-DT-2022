% %% Imperix Cap lifetime model for fixed values //L2=L1×2^(T1−(T2+ΔT))/10 
clc;
% Given data
T1 = 100; % Upper category temperature + temperature increase caused by rated ripple current (in °C)
T2 = linspace(0, 60, 100); % Array of ambient temperatures (in °C), you can adjust the range and number of points
delT = 65; % Temperature increase caused by ripple current (in °C)
L1 = 10000; % Guaranteed life at temperature T1 (in hours)

% Calculate L2 using the provided equations
L2 = L1 * exp((T1 - (T2 + delT)) / 10);

% Plotting
plot(T2, L2, 'b', 'LineWidth', 2);
xlabel('Ambient Temperature (°C)');
ylabel('Expected Life (hours)');
title('Expected Life vs Ambient Temperature');
grid on;

%% T1 ABD delT value change
% Given data
clear;
% Given data
T1_values = [80, 85, 90]; % Upper category temperature + temperature increase caused by rated ripple current (in °C)
delT_values = [40, 55, 70]; % Temperature increase caused by ripple current (in °C)// internal temperature
L1=10000;
% Array of ambient temperatures (T2) in Celsius
T2 = linspace(0, 60, 100);

figure; % Create a new figure

hold on; % Enable hold for multiple plots in the same figure

% Loop through different values of T1 and delT
for i = 1:length(T1_values)
    T1 = T1_values(i);
    delT = delT_values(i); % Use the corresponding delT for each T1
    
    % Calculate L2 using the provided equations
    L2 = L1 * exp((T1 - (T2 + delT)) / 10);
    
    % Plot L2 against T2
    plot(T2, L2, 'LineWidth', 2);
end

xlabel('Ambient Temperature (°C)');
ylabel('Expected Lifetime (hours)');
title('Expected Life of Imperix Capacitor Bank');
grid on;

% Create legends for the three combinations of T1 and delT
legends = cell(length(T1_values), 1);
for i = 1:length(T1_values)
    legends{i} = sprintf('T1 = %d°C, ΔT = %d°C', T1_values(i), delT_values(i));
end

% Add legends to the plot
legend(legends, 'Location', 'best');

hold off; % Disable hold to finish plotting



%% General lifetime model  𝐿_𝑥=𝐿_0×(𝑉_𝑥/𝑉_0 )^(−𝑛)×exp⁡[(𝐸_𝑎/𝐾_𝐵 )(1/𝑇_𝑥 −1/𝑇_0 )]
% % Define the parameters

clc;
clear;
L0 = 10000; % expected lifetime for full rated voltage and temperature
V0 = 400;  % updated rated voltage
T0 = 60; % rated temperature
Ea = 0.94; % activation energy
Kb = 8.62e-5; % Boltzmann's constant

% Create a vector of actual applied voltages 
Vx_values = [350, 400, 450]; % Values of Vx
Tx = linspace(0,60, 10); % Create a vector of maximum ambient temperatures 

% Initialize a matrix to store the calculated lifetimes for each Vx value
Lx_values = zeros(length(Vx_values), length(Tx));

% Calculate the expected operating lifetime for each combination of Vx and Tx
n = 3; % exponent for voltage dependence
for i = 1:length(Vx_values)
    Vx = Vx_values(i);
    Lx = zeros(1, length(Tx)); % initialize the vector of Lx for current Vx
    for j = 1:length(Tx)
%         Lx(j) = L0 * (Vx/V0)^(-n) * exp((Ea/Kb) * (1/Tx(j) - 1/T0));
 Lx(j) = L0 * (Vx/V0)^(-n) * (2.^((T0-Tx(j))/10)); % // L_x=L_0×(V_x/V_0 )^(-n)×2^((T_0-T_x)/10)
    end
    Lx_values(i, :) = Lx; % Store the calculated lifetimes for the current Vx
end

% Plot Lx vs Tx for different Vx values
figure;
for i = 1:length(Vx_values)
    plot(Tx , Lx_values(i, :), 'DisplayName', sprintf('Vx = %.1f V', Vx_values(i)));
    hold on;
end
hold off;

xlabel('Maximum ambient temperature (°C)');
ylabel('Expected operating lifetime (h)');
title('Ambient temperature changes and Vo=400V ');
legend('Location', 'best');
grid on;

