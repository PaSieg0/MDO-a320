%% Example Usage of loads() and aero() Functions
clear all; close all; clc;

% Add Q3D directory to path
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
addpath(q3dPath);

%% Define Wing Geometry Parameters
b = 34.1;           % Total wingspan (m)
c_r = 7.2;          % Root chord (m)
c_k = 4.5;          % Kink chord (m)
c_t = 2.0;          % Tip chord (m)
b_k = 5.5;          % Spanwise location of kink (m)
sweep_te_k = 15;    % Trailing edge sweep angle at kink (deg)
dihedral = 5;       % Dihedral angle (deg)
twist_r = 0;        % Twist at root (deg)
twist_k = -1;       % Twist at kink (deg)
twist_t = -3;       % Twist at tip (deg)

%% Define Flight Conditions
h_cr = 10000;       % Cruise altitude (m)
V_MO_ref = 250;     % Speed (m/s) ~ Mach 0.82
V_cr = 233;    % Cruise speed (m/s)
n_max = 2.5;        % Max load factor (for STRUCTURAL analysis)
n_cruise = 1.0;     % Cruise load factor (for AERO analysis)

% Weight breakdown (in Newtons)
W_AminusW = 400000; 
W_wing = 80000;     
W_fuel = 150000;    

CST = [0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797];

%% 1. Run loads() (Inviscid, 2.5g Structural Limit)
fprintf('=== Running loads() function (Inviscid, 2.5g) ===\n');
% We use n_max here because we need to know if the wing breaks at max load
[Y_loads, L_loads, M_loads] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                                     n_max, V_MO_ref, W_AminusW, h_cr,...
                                     b, c_r, c_k, c_t, CST, W_wing, W_fuel);
fprintf('Loads analysis complete!\n\n');

%% 2. Run aero() (Viscous, 1.0g Cruise Performance)
fprintf('=== Running aero() function (Viscous, 1.0g Cruise) ===\n');

% CRITICAL FIX: Pass 'n_cruise' (1.0), not 'n_max'. 
% Viscous solvers often diverge at 2.5g transonic conditions due to stall.
[Cl, Cd] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                n_cruise, V_cr, W_AminusW, h_cr,...
                b, c_r, c_k, c_t, CST, W_wing, W_fuel);

fprintf('Aero analysis complete!\n\n');

%% Display Results Summary
fprintf('=== RESULTS SUMMARY ===\n');
fprintf('AERO (Performance at 1g Cruise):\n');
fprintf('  CL: %.4f\n', Cl);
fprintf('  CD: %.4f\n', Cd);
fprintf('  L/D: %.2f\n', Cl / Cd);