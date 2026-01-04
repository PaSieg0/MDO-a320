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

% Planform dimensions (in meters) - Source: Assignment specification
b = 34.0;          % Total wingspan (m) [Specification]
c_r = 7.0;          % Root chord (m) [Drawing]
c_k = 3.7;         % Kink chord (m) [Drawing]
c_t = 1.6;          % Tip chord (m) [Specification]
b_k = 4.36 + 3.95/2;          % Spanwise location of kink (m) [Estimated, drawing]

% Wing angles (in degrees)
sweep_te_k = 0.01;    % Trailing edge sweep angle at kink (deg) [Estimated for aerodynamic efficiency]
dihedral = 5;       % Dihedral angle (deg) [Typical A320 value]

% Twist (Washout) - CRITICAL FOR CONVERGENCE
% Tips must be twisted down to prevent shock-induced tip stall at Mach 0.78
twist_r = 0;        % Twist at root (deg) [Standard]
twist_k = -1.0;     % Twist at kink (deg) [Typical washout]
twist_t = -3;       % Twist at tip (deg) [Typical washout]

% Flight parameters - Source: Assignment specification
M_cr = 0.78;        % Cruise Mach number [-] [Specification]
M_MO_ref = 0.82;   % Maximum operating Mach number [-] [Specification]
h_cr = 11278.4;     % Cruise altitude (m) = 37,000 ft [Specification: 37000 ft]
n_max = 2.5;        % Maximum load factor (Structural) [CS-25 requirement]
n_cruise = 1.0;     % Cruise load factor (Performance)

% Calculate cruise speed from Mach number using ISA model
[~, ~, T_cr] = stdatm(h_cr);  % Get temperature at cruise altitude
gamma = 1.4;        % Ratio of specific heats [-] [ISA]
R = 287.058;        % Specific gas constant (J/kg-K) [ISA]
a_cr = sqrt(gamma * R * T_cr);  % Speed of sound at cruise
V_MO_ref = M_MO_ref * a_cr;  % Maximum operating speed (m/s)

% Weight breakdown (in Newtons)
W_AminusW = 400000; % Aircraft weight minus wing (N) [Estimated for A320 class]
W_wing = 50000;     % Wing weight (N) - initial guess (~5100 kg) [Typical for A320]
W_fuel = 150000;    % Fuel weight (N) [Estimated for medium range]

% Structures parameters
MTOW = (W_AminusW + 2*W_wing + W_fuel) / 9.81; % kg
ZFW = (W_AminusW + 2*W_wing) / 9.81;           % kg
S_ref = (c_r + c_k) * b_k + (c_k + c_t) * (b/2 - b_k); % Area

% Airfoil shape (CST parameters) - Source: Baseline airfoil from assignment
% Thickness-to-chord ratio: 0.13 [Specification]
CST = [0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797];

% Spar and tank
spar_locs = [0.2, 0.6];  % Front and rear spar locations [% chord, typical wing box]
tank_limits = [0, 0.85];  % Fuel tank spanwise limits [% half-span, typical]

% Engine data
engine_data.count = 1;      % Engines per half-wing [Standard twin-engine config]
engine_data.y_location = 4.7;  % Engine spanwise position (m) [~35% span, typical]
engine_data.weight = 1969;  % Engine weight (kg) [CFM56-5B class turbofan]

% Airfoils
airfoils.root = 'b737a';  % Root airfoil profile [Placeholder, similar t/c]
airfoils.kink = 'b737a';  % Kink airfoil profile
airfoils.tip  = 'b737a';

% Materials
m_up = [7.1e10, 2795, 4.8e8, 4.6e8];
m_fr = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.upper = m_up;
mat_props.lower = m_up;
mat_props.front = m_fr;
mat_props.rear  = m_up;

% Calculate cruise speed and atmospheric conditions using ISA
[~, ~, T] = stdatm(h_cr);  % Get temperature at cruise altitude (K)
gamma = 1.4;        % Ratio of specific heats [-]
R = 287.058;        % Specific gas constant (J/kg-K)
a = sqrt(gamma * R * T);  % Speed of sound at cruise altitude (m/s)
V_cr = M_cr * a;    % Cruise speed (m/s)

% Engine reference conditions NEEDS CHECK

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

V_cr_ref = V_cr;        % Reference cruise speed (m/s) - matches actual cruise
h_cr_ref = h_cr;        % Reference cruise altitude (m)
C_T_ref = 1.8639e-4;    % Reference specific fuel consumption (1/s)


%% ============================================================
%  SECTION 2: FUNCTION CALLS
%% ============================================================

% Plot wing geometry for verification
plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, dihedral);
plot_airfoil(CST);

% Evaluate fuel tank volume (performance)
W_fuel = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);
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
    MTOW = (W_AminusW + 2*W_wing + W_fuel) / 9.81;
    ZFW = (W_AminusW + 2*W_wing) / 9.81;
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
                M_cr, W_AminusW, h_cr,...
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
fprintf('Total Aircraft Weight: %.2f kg\n', MTOW);
fprintf('Fuel Capacity:      %.2f kg\n', W_fuel / 9.81);
fprintf('Lift Coefficient:   %.4f\n', Cl);
fprintf('Drag Coefficient:   %.4f\n', Cd);
fprintf('L/D Ratio:          %.2f\n', L_D_ratio);
fprintf('Mission Range:      %.2f km\n', R / 1000);
fprintf('Fuel Consumption:   %.2f kg\n', W_fuel_calc / 9.81);
fprintf('Iterations:         %d\n', iter);
fprintf('==========================================\n\n');