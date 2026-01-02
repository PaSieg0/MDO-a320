% runPerformanceAnalysis.m
%
% This script sets up a test case for an A320-like aircraft to demonstrate 
% the use of the calculatePerformance function, which is based on the 
% "AE4-205 MDO for Aerospace Applications 2025/2026" homework assignment.

clear; clc; close all;

%% 1. Define Input Parameters for a Reference Aircraft Case
fprintf('--- Defining Input Parameters for Test Case ---\n');

% --- Operational Conditions (Variables in the MDO problem)
Mcr = 0.78;         % Cruise Mach number
h = 11000;          % Cruise altitude (m)

% --- Atmospheric properties at cruise altitude (ISA model)
T0 = 288.15;        % Sea level temperature (K)
a0 = 340.3;         % Sea level speed of sound (m/s)
L = -0.0065;        % Temperature lapse rate (K/m) in troposphere
T = T0 + L * h;     % Temperature at altitude h (K)
a = a0 * sqrt(T / T0); % Speed of sound at altitude h (m/s)
V = Mcr * a;        % Cruise speed (m/s)

% --- Assumed Aerodynamic and Weight Properties
L_D_ratio = 16;     % Assumed L/D ratio (as suggested on page 5)
W_TO_max_kg = 78000;% Maximum Take-Off Weight (kg)
g = 9.80665;        % Acceleration due to gravity (m/s^2)
W_TO_max = W_TO_max_kg * g; % Convert to Newtons

% --- Estimate cruise weights for this example
% In the full MDO problem, these would be coupled variables. Here we assume
% some realistic values for demonstration.
W_start_cr = W_TO_max * 0.98; % Weight after climb
W_end_cr   = W_TO_max * 0.82; % Weight before descent

% --- Reference Aircraft Parameters (Constants from assignment/datasheet)
V_cr_ref = 0.78 * a; % Reference cruise speed (m/s), assumed to be the same for this case
h_cr_ref = 11000;    % Reference cruise altitude (m)
C_T_ref = 1.8639e-4; % Reference SFC (N/Ns or 1/s), from Page 6

fprintf('Cruise Conditions:\n');
fprintf('  - Mach Number: %.2f\n', Mcr);
fprintf('  - Altitude: %d m\n', h);
fprintf('  - Speed: %.2f m/s\n', V);
fprintf('Weights (N):\n');
fprintf('  - W_TO_max:   %.0f N\n', W_TO_max);
fprintf('  - W_start_cr: %.0f N\n', W_start_cr);
fprintf('  - W_end_cr:   %.0f N\n', W_end_cr);
fprintf('--------------------------------------------------\n\n');


%% 2. Call the Performance Calculation Function
fprintf('--- Calculating Performance ---\n');
[R, W_fuel_calc, C_T] = calculatePerformance(V, h, L_D_ratio, ...
    W_start_cr, W_end_cr, W_TO_max, ...
    V_cr_ref, h_cr_ref, C_T_ref);


%% 3. Display the Results
fprintf('--- Performance Calculation Results ---\n');
fprintf('Calculated Propulsive Efficiency Factor (eta): %.4f\n', C_T_ref / C_T);
fprintf('Calculated Specific Fuel Consumption (C_T): %e 1/s\n', C_T);
fprintf('\n');
fprintf('OBJECTIVE -> Mission Range (R): %.2f km\n', R / 1000);
fprintf('\n');
fprintf('Calculated Total Fuel Weight (W_fuel): %.2f N (%.2f kg)\n', W_fuel_calc, W_fuel_calc / g);
fprintf('-----------------------------------------\n');