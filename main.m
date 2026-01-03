%% Main MDO Analysis Script for A320 Wing Design
clear all;
close all;
clc;

% Setup paths
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
emwetPath = fullfile(currentDir, 'EMWET 1.5');
addpath(q3dPath);
addpath(emwetPath);

%% ============================================================
%  SECTION 1: DATA DEFINITION
%% ============================================================

% Planform dimensions (in meters)
b = 34.1;           % Total wingspan (m)
c_r = 7.2;          % Root chord (m)
c_k = 4.5;          % Kink chord (m)
c_t = 2.0;          % Tip chord (m)
b_k = 5.5;          % Spanwise location of kink (m)

% Wing angles (in degrees) - FIXED FOR CONVERGENCE
% 0.1 deg sweep causes shock issues. Increased to A320 realistic values.
sweep_te_k = 10;    % Trailing edge sweep angle at kink (deg) 
dihedral = 5;       % Dihedral angle (deg)

% Twist (Washout) - CRITICAL FOR CONVERGENCE
% Tips must be twisted down to prevent shock-induced tip stall at Mach 0.78
twist_r = 0;        % Twist at root (deg)
twist_k = -1.0;     % Twist at kink (deg)
twist_t = -3;     % Twist at tip (deg)

% Flight parameters
h_cr = 10000;       % Cruise altitude (m)
V_MO_ref = 250;     % Maximum operating speed (m/s)
n_max = 2.5;        % Maximum load factor (Structural)
n_cruise = 1.0;     % Cruise load factor (Performance)

% Weight breakdown (in Newtons)
W_AminusW = 400000; % Aircraft weight minus wing (N)
W_wing = 50000;     % Wing weight (N) - initial guess (~5100 kg, reasonable for A320)
W_fuel = 150000;    % Fuel weight (N)

% Structures parameters
MTOW = (W_AminusW + W_wing + W_fuel) / 9.81; % kg
ZFW = (W_AminusW + W_wing) / 9.81;           % kg
S_ref = 122.6;      % Wing reference area (m^2)

% Airfoil shape (CST parameters)
CST = [ 0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797 ];

% Spar and tank
spar_locs = [0.15, 0.65];
tank_limits = [0.1, 0.9];

% Engine data
engine_data.count = 1;
engine_data.y_location = 4.7;
engine_data.weight = 1969;

% Airfoils
airfoils.root = 'b737a';
airfoils.kink = 'b737a';
airfoils.tip  = 'b737a';

% Materials
m_up = [7.1e10, 2795, 4.8e8, 4.6e8];
m_fr = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.upper = m_up;
mat_props.lower = m_up;
mat_props.front = m_fr;
mat_props.rear  = m_up;

% Atmospheric constants
T0 = 288.15;        % Sea level temperature (K)
a0 = 340.3;         % Sea level speed of sound (m/s)
L_atm = -0.0065;    % Temperature lapse rate (K/m)

% Cruise conditions
Mcr = 0.78;         % Cruise Mach number

% Calculate cruise speed and atmospheric conditions
T = T0 + L_atm * h_cr;  % Temperature at cruise altitude (K)
a = a0 * sqrt(T / T0);  % Speed of sound at cruise altitude (m/s)
V_cr = Mcr * a;         % Cruise speed (m/s)

% Engine reference conditions NEEDS CHECK

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

V_cr_ref = V_cr;        % Reference cruise speed (m/s) - matches actual cruise
h_cr_ref = h_cr;        % Reference cruise altitude (m)
C_T_ref = 1.8639e-4;    % Reference specific fuel consumption (1/s)


%% ============================================================
%  SECTION 2: FUNCTION CALLS
%% ============================================================

% Plot wing geometry for verification
% plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, dihedral);
% plot_airfoil(CST);

% Evaluate fuel tank volume (performance)
W_fuel = performance(b, c_r, c_k, c_t, b_k, tank_limits);
fprintf('Initial fuel capacity: %.2f kg\n', W_fuel / 9.81);

%% ============================================================
%  CONVERGENCE LOOP
%% ============================================================

% Convergence parameters
max_iter = 20;
tol = 0.001;  % 0.1% convergence tolerance
converged = false;
iter = 0;

fprintf('\n========== STARTING MDO CONVERGENCE LOOP ==========\n');

while ~converged && iter < max_iter
    iter = iter + 1;
    fprintf('\n--- Iteration %d ---\n', iter);
    
    W_wing_old = W_wing;
    
    % Loads analysis - Use actual loads() function
    [Y, L, M] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                      n_max, V_MO_ref, W_AminusW, h_cr,...
                      b, c_r, c_k, c_t, CST, W_wing, W_fuel);

    % Structures analysis
    W_wing_new = 9.81 * structures(Y, L, M, ...
                            sweep_te_k, b_k, dihedral, ...
                            b, c_r, c_k, c_t, ...
                            MTOW, ZFW, 1, S_ref, ...
                            mat_props, spar_locs, ...
                            tank_limits, engine_data, airfoils, ...
                            'a320_main');
    
    err_wing = abs(W_wing_new - W_wing_old) / W_wing_old;
    
    fprintf('  W_wing: %.2f kg (change: %.2f%%)\n', W_wing_new/9.81, err_wing*100);
    
    if err_wing < tol
        converged = true;
        fprintf('>>> CONVERGED after %d iterations <<<\n', iter);
    end
    
    % Update for next iteration
    W_wing = W_wing_new;
    MTOW = (W_AminusW + W_wing + W_fuel) / 9.81;
    ZFW = (W_AminusW + W_wing) / 9.81;
end

if ~converged
    fprintf('>>> WARNING: Did not converge after %d iterations <<<\n', max_iter);
end

fprintf('===================================================\n\n');

%% ============================================================
%  CALCULATE RANGE (POST-CONVERGENCE)
%% ============================================================

% Aerodynamic analysis
[Cl, Cd] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                n_cruise, V_cr, W_AminusW, h_cr,...
                b, c_r, c_k, c_t, CST, W_wing_new, W_fuel);

% Calculate L/D ratio from aerodynamic coefficients
L_D_ratio = Cl / Cd;

% Define weights for cruise segment
W_TO_max = W_AminusW + 2 * W_wing_new + W_fuel;
W_start_cr = W_TO_max;  % Weight at start of cruise (after taxi, takeoff, climb)
W_end_cr = (1 - W_fuel / W_TO_max) * W_start_cr / (0.938);

% Calculate performance using Breguet equations
[R, W_fuel_calc, C_T] = calculatePerformance(V_cr, h_cr, L_D_ratio, ...
    W_start_cr, W_end_cr, W_TO_max, ...
    V_cr_ref, h_cr_ref, C_T_ref);

%% ============================================================
%  SECTION 4: DISPLAY RESULTS
%% ============================================================

fprintf('\n========== MDO ANALYSIS RESULTS ==========\n');
fprintf('Wing Weight:        %.2f kg\n', W_wing_new / 9.81);
fprintf('Fuel Capacity:      %.2f kg\n', W_fuel / 9.81);
fprintf('Lift Coefficient:   %.4f\n', Cl);
fprintf('Drag Coefficient:   %.4f\n', Cd);
fprintf('L/D Ratio:          %.2f\n', L_D_ratio);
fprintf('Mission Range:      %.2f km\n', R / 1000);
fprintf('Fuel Consumption:   %.2f kg\n', W_fuel_calc / 9.81);
fprintf('Iterations:         %d\n', iter);
fprintf('==========================================\n\n');