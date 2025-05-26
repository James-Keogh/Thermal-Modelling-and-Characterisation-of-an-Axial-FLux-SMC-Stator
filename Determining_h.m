clear all
close all
clc

%% === Importing Experimental Data ===


Experimental_data_coil_1_insulated_core = readmatrix("Formatted_Insulated_Core_10A_DC_Run_3.xlsx");

L = length(Experimental_data_coil_1_insulated_core(:,1));

A10_Sensor_1to4 = Experimental_data_coil_1_insulated_core(6:L, 1:4)';
A10_T_amb = Experimental_data_coil_1_insulated_core(6:L, 5)';

A10_Current = Experimental_data_coil_1_insulated_core(6:L, 7);
A10_Voltage = Experimental_data_coil_1_insulated_core(6:L, 9);
A10_Power = Experimental_data_coil_1_insulated_core(6:L, 8);

A10_Functional_sensors = [A10_Sensor_1to4(1,:); A10_Sensor_1to4(3,:); A10_Sensor_1to4(4,:)];
A10_Mean_Sensor = [];

% Creating a Mean value of all three sensors

for n = 1:length(A10_Functional_sensors)
    A10_Mean_Sensor(end+1)= mean(A10_Functional_sensors(:,n));
end

T_amb = mean(A10_T_amb);                        
t = linspace(1,length(A10_T_amb),length(A10_T_amb));

% Finding the peak temperature / point of switching off power

[max_temp, t_peak] = max(A10_Mean_Sensor); 

cooling_start_index = t_peak;                   
cooling_time = t(cooling_start_index:end);
cooling_temp = A10_Mean_Sensor(cooling_start_index:end);

%% === Exponential Fit Method ===


exp_model = @(b, t) T_amb + (cooling_temp(1) - T_amb) * exp(-t / b);
obj_fun = @(b) sum((exp_model(b, cooling_time - cooling_time(1)) - cooling_temp).^2);
tau_initial = 100;
tau_opt = fminsearch(obj_fun, tau_initial);

% Coil geometry

C = 2 * 8.0128;              % J/K (capacitance of copper in coil)
A = 2 * 863.929106e-6;       % m² (coil surface area, doubled)

h_tau = C / (A * tau_opt);   % W/m²·K based on exponential fit

fitted_curve = exp_model(tau_opt, cooling_time - cooling_time(1));

%% === Per-Step Euler Estimation of h ===


dt = 1;  % time step (assuming 1s between data points; adjust if different)

h_values = zeros(1, length(cooling_temp)-1);

for i = 1:length(cooling_temp)-1
    dT = cooling_temp(i+1) - cooling_temp(i);
    T_diff = cooling_temp(i) - T_amb;
    
    if abs(T_diff) > 1e-4

        h_values(i) = - (C / (A * dt)) * (dT / T_diff);
    else
        h_values(i) = NaN;
    end
end

h_time = cooling_time(1:end-1);

%% === Plotting Results ===


figure;
plot(cooling_time, cooling_temp, 'b', 'DisplayName', 'Experimental',LineWidth=2);
hold on
plot(cooling_time, fitted_curve, 'r--', 'DisplayName', ['Fitted (τ = ' num2str(tau_opt, '%.1f') 's, h = ' num2str(h_tau, '%.1f') ' W/m²K)'],LineWidth=2);
fontsize(24,"points")
xlabel('Time in s');
ylabel('Temperature in C');
legend;
title('Cooling Curve and Exponential Fit (Method 1)');
grid on;

figure;
plot(h_time, h_values, 'k', 'DisplayName', 'Per-step Euler h');
xlabel('Time (s)');
ylabel('h (W/m²·K)');
fontsize(24,"points")
title('Per-Step Convection Coefficient Estimation');
legend;
grid on;
hold off

%% === Finding Average of h Value Evolution ===

windowSize = 20;  % adjust depending on noise level
h_smooth = movmean(h_values, windowSize);

figure
plot(h_time, h_values, 'k:', 'DisplayName', 'Raw h');
hold on;
plot(h_time, h_smooth, 'r', 'LineWidth', 2, 'DisplayName', 'Smoothed h (Moving Avg)',LineWidth=2);
xlabel('Time in s');
ylabel('h in W/m²·K');
fontsize(24,"points")
title('Per-Step Convection Coefficient Estimation (Method 2)');
legend("Per-step determination of h","Smoothed h curve")
grid on;