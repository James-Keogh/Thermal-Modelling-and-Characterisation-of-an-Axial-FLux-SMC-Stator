clear all
close all
clc

%% === Importing Experimental Data ===


Experimental_data_AC = readmatrix("Formatted_AC_100Hz_13A_RMS.xlsx");

AC_Sensor_1to6 = Experimental_data_AC(7:724, 5:10)';                    % Sensor data of stator sensors
AC_T_amb = Experimental_data_AC(7:724, 4)';                             % Ambient sensor data
AC_Sensor_ABC = Experimental_data_AC(7:724, 1:3)';                      % Sensor data of sensors in the coil
AC_Mean_Sensor = [];                                                    % Initialising a mean coil sensor vector


% These values are needed to linearly ramp up the DC Joule loss component
% as the temeprature increases in the coil

A10_V_in = 1.26;                                % Inital voltage measured 
A10_V_fin = 1.42;                               % Final voltage measured

% Creating a Mean value of all three coil sensors

for n = 1:length(AC_Sensor_ABC)
    
    AC_Mean_Sensor(end+1)= mean(AC_Sensor_ABC(:,n));

end

%% === Importing Thermal Conductances ===


Conduction_Conductances = readmatrix("overview sheet AC.xlsx","Range","C27:C47")';
Conduction_Names = readcell("overview sheet AC.xlsx","Range","B27:B47")';

Convection_Conductances = readmatrix("overview sheet AC.xlsx","Range","N27:N42")';
Convection_Names = readcell("overview sheet AC.xlsx","Range","M27:M42")';

Capacitance_Values = readmatrix("overview sheet AC.xlsx","Range","C52:C67")';
Capacitance_Names = readcell("overview sheet AC.xlsx","Range","B52:B67")';

Nodal_Masses = readmatrix("overview sheet AC.xlsx","Range","D52:D67")';

%% === Automatically Allocating Variable Names and Values ===


for i = 1:length(Conduction_Names)
    name = Conduction_Names{i};
    value = Conduction_Conductances(i);
    eval([name ' = ' num2str(value) ';']);
end


for i = 1:length(Convection_Names)
    name = Convection_Names{i};
    value = Convection_Conductances(i);
    eval([name ' = ' num2str(value) ';']);
end


for i = 1:length(Capacitance_Names)
    name = Capacitance_Names{i};
    value = Capacitance_Values(i);
    eval([name ' = ' num2str(value) ';']);
end

%% === Defining Parameters ===


N = 16;                                                                % Number of nodes in the system

T_ambient = mean(AC_T_amb);                                            % Ambient temperature in C taken from experimental data

G_conv = Convection_Conductances;                                      % Manually inserting all nodes that have convection

C_vector = Capacitance_Values';                                        % Heat capacities per node (J/K)

%% === Time Settings ===


t = linspace(0, 718, 100000);                                           % Time vector
dt = t(2) - t(1);                                                       % Time step

%% === Power Input Setup ===


P_d = zeros(1,length(t));                                               % Initialising the power vector

[max_temp, t_peak] = max(AC_Mean_Sensor);                               % Determining the peak time and temperature from experimental data 

Power_increase_percent = A10_V_fin/A10_V_in;                            % Finding the percentual increase of power over the measured heating phase

time_step_factor = 1/dt;                                                % We introduce this facotr because t_peak is in seconds, but our time step in the t vector is much smaller to avoid instability

time_for_power_increase = round(t_peak*time_step_factor);               % we need an integer to replace values in a vector

%% !!!!!!!!!! INPUT !!!!!!!!!!


I_DC = 13;                                                              % Amount of current used in testing (RMS or DC)
f = 100;                                                                % AC frequncy in Hz, set to 0 in DC

%% !!!!!!!!!! INPUT !!!!!!!!!!



% Coil Parameters for power dissipation

resistivity_20C = 1.724*10^-8;                                          % in ohm*m resistivity at 20C
L_cu_tot = 5.54544130;                                                  % in m Total length of copper wire in the coil
D_wire = 1.217*10^-3;                                                   % in m Wire diameter
A_cu = pi*(0.5*D_wire)^2;                                               % in m^2 Cross-sectional area of the wire from coil data

R_el = resistivity_20C*(L_cu_tot/A_cu);                                 % Calculating Joule losses dissipated from coil properties
P_d_initial = ((I_DC^2)*R_el)/2;


% Setting dissipated power as a linear increase matched with the measured results 

P_d(1:time_for_power_increase) = linspace(P_d_initial,P_d_initial*Power_increase_percent,time_for_power_increase);  



%% === Calcuating Iron Losses ===


% Experimentally determined vlaues for the Steinmetz equation

K_h = 1.0977375*10^-1;
K_e = 4.4280188*10^-5;
alpha = 1.75;



B_peak = 1.3;                                                       % in T, the peak magnetic flux density, set to 0 if in DC

Iron_loss = ((K_h*f*B_peak^alpha) + (K_e*(f^2)*(B_peak)^2));        % in W/kg

P_Fe = zeros(N,length(t));

% Empirical Steinmetz equation

for i = 3:12                                                        % We only go from 3-12 becaue we only want to add iron losses to  
                                                                    % the stator nodes along the magnetic circuit with a non 0 volume
    P_Fe(i,1:time_for_power_increase) = Nodal_Masses(i)*Iron_loss;
end

% Adding the iron losses to the general power vector

Power_Matrix = [P_d; P_d; P_Fe(3:end,:)];
%% === Conductance Matrix Setup ===



G_matrix = zeros(N,N);                                               % Initialize G matrix

% Define conduction connections (modify as needed)

% Ensure that each connection is only mentioned once. e.g if [1 3] is present [3 1] should not be added

connections = [1 3; 
    1 6; 
    1 7; 
    2 5; 
    2 6; 
    2 7; 
    3 4;
    4 6;
    4 7;
    5 4;
    6 15;
    7 8;
    7 9;
    8 9;
    8 10;
    9 10;
    10 11;
    11 12;
    12 16;
    7 13;
    10 14];                                                            

% Allocating the cunductive conductances within the matrix

for k = 1:size(connections,1)
    i = connections(k,1);
    j = connections(k,2);
    G_matrix(i,j) = -Conduction_Conductances(k);
    G_matrix(j,i) = -Conduction_Conductances(k);
    G_matrix(i,i) = G_matrix(i,i) + Conduction_Conductances(k);
    G_matrix(j,j) = G_matrix(j,j) + Conduction_Conductances(k);         % We modify the j,j node because our connections matrix only counts each connection once
                                                                        % If we had [1 2] AND [2 1] in the connections matrix then this step would be unnecessary
                                              
end

% Add convection to the diagonal. 

for i = 1:N
    G_matrix(i,i) = G_matrix(i,i) + G_conv(i);                          % Adding the convection correspondant to each node 
end

%% === Forward Euler Differencing ===


% Initializing all nodes to initial sensor temperatures

T = AC_Mean_Sensor(1) * ones(N, length(t));                             

T(15,:) = AC_Sensor_1to6(1,1);
T(3,:) = AC_Sensor_1to6(2,1);
T(7,:) = AC_Sensor_1to6(3,1);
T(9,:) = AC_Sensor_1to6(4,1);
T(10,:) = AC_Sensor_1to6(5,1);
T(16,:) = AC_Sensor_1to6(6,1);


for n = 1:length(t)-1
    current_time = t(n);
    
    % Dynamic power vector

    P_vector_dynamic = Power_Matrix(:,n);

    % Compute net heat flow per node and update temperatures

    Q = P_vector_dynamic - G_matrix * T(:,n) + G_conv' .*T_ambient;
    T(:,n+1) = T(:,n) + (Q ./ C_vector) * dt;
end

%% === Plotting ===

Simulated_sensors = [T(1,:);
T(15,:);
T(3,:);
T(7,:);
T(9,:);
T(10,:);
T(16,:)];

% Plotting only the nodes that also have sensors, to compare to experimental results

figure

hold on
grid on
plot(t, Simulated_sensors,Linewidth=2)

sensor_names = ["AC Simulated ABC mean"];

for i = 1:6
    sensor_names = [sensor_names, ("AC Simulated Sensor " + num2str(i))];
end

%% === AC 13A RMS 100Hz Tests ===


t_exp = linspace(1,length(AC_T_amb),length(AC_T_amb));

plot(t_exp,AC_Mean_Sensor,Linestyle = '--',Linewidth=2)

plot_legend = [sensor_names, "AC ABC Mean"];

for n = 1:length(AC_Sensor_1to6(:,1))

    plot(t_exp,AC_Sensor_1to6(n,:),Linestyle = '--',Linewidth=2 )
    plot_legend = [plot_legend, "AC Sensor " + num2str(n)];

end


plot_legend = [plot_legend, "AC Ambient"];

plot(t_exp, AC_T_amb, Linestyle = '--', Linewidth=2)
fontsize(24,"points")
xlabel('Time in Seconds')
ylabel('Temperature in Celsius')
title('Generalized Nodal Temperature Evolution AC 13A RMS, 100 Hz')
legend(plot_legend)
grid on
hold off

%% === Error Graphing ===

tiledlayout(1,2)
nexttile

% We do a few strange things here just to make the error vectors the same
% lengths. The timestep in the experimental and simulated model are very
% different, thats the reason for all this.

abserror = absolute_error(T(1,1:round(time_step_factor):end-round(2*time_step_factor)),AC_Mean_Sensor(1,:));
pererror = Percentual_error(T(1,1:round(time_step_factor):end-round(2*time_step_factor)),AC_Mean_Sensor(1,:));

plot(t_exp,abserror, LineWidth=2)
fontsize(24,"points")
title("Error of the Coil Temperature 13A RMS AC")
ylabel("Absolute Error in Celsius")
xlabel("Time in Seconds")
yyaxis right
plot(t_exp, pererror, LineWidth=2)
ylabel("Percentual Error w.r.t Measured Values")
legend('Absolute Error','Percentual Error')
grid minor



abserror = absolute_error(T(3,1:round(time_step_factor):end-round(2*time_step_factor)),AC_Sensor_1to6(2,:));
pererror = Percentual_error(T(3,1:round(time_step_factor):end-round(2*time_step_factor)),AC_Sensor_1to6(2,:));

nexttile
plot(t_exp,abserror, LineWidth=2)
fontsize(24,"points")
title("Error of the Sensor 2 Temperature 13A RMS AC")
ylabel("Absolute Error in Celsius")
xlabel("Time in Seconds")
yyaxis right
plot(t_exp, pererror, LineWidth=2)
ylabel("Percentual Error w.r.t Measured Values")
legend('Absolute Error','Percentual Error')
grid minor

%% 
% *===  Error Determining Functions  ===*

function [abserror] = absolute_error(Modelled_Value, Measured_Value)

    abserror = abs(Modelled_Value-Measured_Value);

end


function[pererror] = Percentual_error(Modelled_Value, Measured_Value)

    pererror = (abs(Measured_Value-Modelled_Value))./Measured_Value;

    pererror = pererror.*100;                   % turning fraction into percent values
end