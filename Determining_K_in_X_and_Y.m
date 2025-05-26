clear all 
close all
clc

%% === Importing Experimental Data ===

Data_510_Press_Run_3 = readmatrix("Formatted_510_Press_Direction_Run_3.xlsx")';

L = length(Data_510_Press_Run_3(1,:));

Heated_Sensor_y = Data_510_Press_Run_3(1,15:322);
Cool_Sensor_y = Data_510_Press_Run_3(2,15:322);

Applied_Voltage_y = Data_510_Press_Run_3(6,15:322);
Active_Power_y = Data_510_Press_Run_3(5,15:322);
Applied_Current_y = Data_510_Press_Run_3(4,15:322);



% Doing the same again for the Perpendicular Direction



Data_510_Perpendicular_Run_1 = readmatrix("Formatted_510_Perpendicular_Direction_Run_1.xlsx")';

L = Data_510_Perpendicular_Run_1(1,:);

Heated_Sensor_x = Data_510_Perpendicular_Run_1(1,7:365);
Cool_Sensor_x = Data_510_Perpendicular_Run_1(2,7:365);

Applied_Voltage_x = Data_510_Perpendicular_Run_1(6,7:365);
Active_Power_x = Data_510_Perpendicular_Run_1(5,7:365);
Applied_Current_x = Data_510_Perpendicular_Run_1(4,7:365);

%% === Defining the System Paramters ===


Cube_Side_Length = 22*10^-3;            % in m
Cube_Area = Cube_Side_Length^2;         % in m^2
Cube_Volume = Cube_Side_Length^3;       % in m^3
Density = 7450;                         % in kg/m^3
Cube_Mass = Density*Cube_Volume;        % in kg
Delta_t = 1;                            % in s, comes from 1 Hz sample frequency of test data

c_p = 408.8;                            % Determined from previous optimisation for SMC
C_th = Cube_Mass*c_p;                   

%% === Calculating Ky ===

K_y = [];

% Rearanging the forward euler differencing algorithm

for i = 1:length(Heated_Sensor_y)-1
    
    term_1 = (Cool_Sensor_y(i+1)-Cool_Sensor_y(i)) / (Heated_Sensor_y(i)-Cool_Sensor_y(i));
    term_2 = (Cube_Side_Length*C_th) / (Cube_Area*Delta_t);
    
    K_y(end+1) = term_1*term_2;
end

K_y_Smoothed = movmean(K_y, 5);                     % Creating a smoothed curve

Final_K_y = mean(K_y_Smoothed(100:300))             % Final mean value

t_y = linspace(1,length(K_y),length(K_y));


%% === Calculating Kx ===

K_x = [];

% Rearanging the forward euler differencing algorithm

for i = 1:length(Heated_Sensor_x)-1
    
    term_1 = (Cool_Sensor_x(i+1)-Cool_Sensor_x(i)) / (Heated_Sensor_x(i)-Cool_Sensor_x(i));
    term_2 = (Cube_Side_Length*C_th) / (Cube_Area*Delta_t);
    
    K_x(end+1) = term_1*term_2;
end

K_x_Smoothed = movmean(K_x, 5);                     % Creating a smoothed curve

Final_K_x = mean(K_x_Smoothed(100:300))             % Final mean value

t_x = linspace(1,length(K_x),length(K_x));

%% === Plotting ===


hold on

plot(t_y,K_y,LineStyle=":")
plot(t_y,K_y_Smoothed, LineWidth=2)
plot(t_x,K_x,'k',LineStyle=":")
plot(t_x,K_x_Smoothed, LineWidth=2)
fontsize(24,"points")
title("Experimental Determination of K in Different Directions")
xlabel("Time in Seconds")
ylabel("Heat Transfer Coefficient in W/mK")
legend("K Press","K Press Smoothed","K Perpendicular","K Perpendicular Smoothed")
ylim([0 25])
grid minor
hold off

%%


figure; clf;
hold on
plot(t_y,Heated_Sensor_y(1:length(t_y)),'r',LineWidth=2);
plot(t_y,Cool_Sensor_y(1:length(t_y)),'b',LineWidth=2);
plot(t_x,Heated_Sensor_x(1:length(t_x)),LineWidth=2,LineStyle="--");
plot(t_x,Cool_Sensor_x(1:length(t_x)),LineWidth=2,LineStyle="--");
fontsize(24,"points")
title("Temperature Data From K Evaluation Experiments")
xlabel("Time in Seconds")
ylabel("Temperature in C")
legend("Heated side press direction","Cool side press direction","Heated side perpendicular direction","Cool side perpendicular direction",Location="southeast")
grid minor
hold off