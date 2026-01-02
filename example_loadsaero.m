%% Example Usage of loads() and aero() Functions
% This script demonstrates how to call the aerodynamic analysis functions
% with realistic A320-like wing parameters

clear all
close all
clc

% Change to Q3D directory (required for Q3D_solver to find Storage folder)
cd('C:\Users\Jaime\Desktop\TUD\Master\MDO\MDO-a320\Q3D')

%% Define Wing Geometry Parameters

% Planform dimensions (in meters)
b = 34.1;           % Total wingspan (m) - A320 typical
c_r = 7.2;          % Root chord (m)
c_k = 4.5;          % Kink chord (m)
c_t = 2.0;          % Tip chord (m)
b_k = 5.5;          % Spanwise location of kink (m)

% Wing angles (in degrees)
sweep_te_k = 15;    % Trailing edge sweep angle at kink (deg)
dihedral = 5;       % Dihedral angle (deg)
twist_r = 0;        % Twist at root (deg)
twist_k = -1;       % Twist at kink (deg)
twist_t = -3;       % Twist at tip (deg)

%% Define Flight Conditions

% Flight parameters
h_cr = 10000;       % Cruise altitude (m) - approximately 33,000 ft
V_MO_ref = 250;     % Maximum operating speed (m/s) - about Mach 0.82
n_max = 2.5;        % Maximum load factor (positive)

% Weight breakdown (in Newtons)
W_AminusW = 400000; % Aircraft weight minus wing (N) - approx 40,000 kg
W_wing = 80000;     % Wing weight (N) - approx 8,000 kg
W_fuel = 150000;    % Fuel weight (N) - approx 15,000 kg

%% Define Airfoil Shape using CST Parameters
% CST (Class-Shape Transformation) coefficients for the airfoil
% Format: [upper curve coefficients | lower curve coefficients]
% These represent a typical supercritical airfoil

CST = [0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797;
       0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797];

%% Call the loads() function
% This function computes structural loads (lift and pitching moment distributions)
fprintf('=== Running loads() function (inviscid analysis) ===\n');

[Y_loads, L_loads, M_loads] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                                     n_max, V_MO_ref, W_AminusW, h_cr,...
                                     b, c_r, c_k, c_t, CST, W_wing, W_fuel);

fprintf('Loads analysis complete!\n\n');

%% Call the aero() function  
% This function computes overall aerodynamic coefficients (lift and drag)
fprintf('=== Running aero() function (viscous analysis) ===\n');

[Cl, Cd] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                n_max, V_MO_ref, W_AminusW, h_cr,...
                b, c_r, c_k, c_t, CST, W_wing, W_fuel);

fprintf('Aero analysis complete!\n\n');

%% Display Results Summary

fprintf('=== RESULTS SUMMARY ===\n\n');

% Loads results
fprintf('LOADS (Structural) Analysis:\n');
fprintf('  Number of spanwise stations: %d\n', length(Y_loads));
fprintf('  Max lift per unit span: %.2f N/m\n', max(L_loads));
fprintf('  Max moment per unit span: %.2f Nm/m\n', max(abs(M_loads)));
fprintf('  Total semi-span lift: %.2f kN\n', trapz(Y_loads, L_loads)/1000);
fprintf('\n');

% Aero results
fprintf('AERO (Performance) Analysis:\n');
fprintf('  Wing lift coefficient (CL): %.4f\n', Cl);
fprintf('  Wing drag coefficient (CD): %.4f\n', Cd);
fprintf('  L/D ratio: %.2f\n', Cl / Cd);
fprintf('\n');

%% Plotting Results

figure('Position', [100, 100, 1200, 600]);

% Plot 1: Lift Distribution
subplot(2,2,1)
plot(Y_loads, L_loads, 'b-', 'LineWidth', 2);
grid on;
xlabel('Spanwise Position (m)');
ylabel('Lift per unit span (N/m)');
title('Lift Distribution');
legend('Loads Analysis', 'Location', 'best');

% Plot 2: Moment Distribution
subplot(2,2,2)
plot(Y_loads, M_loads, 'b-', 'LineWidth', 2);
grid on;
xlabel('Spanwise Position (m)');
ylabel('Moment per unit span (Nm/m)');
title('Pitching Moment Distribution (about c/4)');
legend('Loads Analysis', 'Location', 'best');

% Plot 3: Integrated Loads
subplot(2,2,3)
cumulative_lift = cumtrapz(Y_loads, L_loads);
cumulative_moment = cumtrapz(Y_loads, M_loads);
yyaxis left
plot(Y_loads, cumulative_lift/1000, 'b-', 'LineWidth', 2);
ylabel('Cumulative Lift (kN)');
yyaxis right
plot(Y_loads, cumulative_moment/1000, 'r-', 'LineWidth', 2);
ylabel('Cumulative Moment (kNm)');
xlabel('Spanwise Position (m)');
title('Cumulative Loads');
grid on;
legend('Cumulative Lift', 'Cumulative Moment', 'Location', 'best');

% Plot 4: Aerodynamic Coefficients Bar Chart
subplot(2,2,4)
bar_data = [Cl, Cd, Cl/Cd];
bar(bar_data);
set(gca, 'XTickLabel', {'CL', 'CD', 'L/D'});
ylabel('Value');
title('Wing Aerodynamic Coefficients');
grid on;
for i = 1:length(bar_data)
    text(i, bar_data(i), sprintf('  %.4f', bar_data(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

sgtitle('Wing Aerodynamic Analysis Results', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Plots generated successfully!\n');
