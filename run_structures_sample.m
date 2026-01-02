% Script to run structures.m with sample data
clc; clear; close all;

%% 1. Define Sample Inputs
% Note: caseName will default to 'a320' inside the function if we don't pass it.
myCaseName = 'test_run'; 

% Weights
MTOW  = 52390;
ZFW   = 46720;
n_max = 2.5;

% Geometry Inputs (Naming convention from loads/aero)
S_ref = 91.04;
b     = 28.35;
b_k   = 4.8;
c_r   = 7.38;
c_k   = 5.00;
c_t   = 1.51;

% Angles
% Converting Leading Edge Sweep (25 deg) to Trailing Edge Sweep for inputs
Lambda_LE_in = 25; 
x_le_k_equiv = b_k * tan(deg2rad(Lambda_LE_in));
x_te_k_equiv = x_le_k_equiv + c_k;
x_te_r       = c_r;
sweep_te_k   = rad2deg(atan((x_te_k_equiv - x_te_r) / b_k));
dihedral     = 0;

% Structs for organized data passing
spar_locs = [0.15, 0.65]; % Front, Rear

tank_limits = [0.1, 0.9];

engine_data.count = 1;
engine_data.y_location = 4.7;
engine_data.weight = 1969;

airfoils.root = 'a320'; % Ensure .dat files exist
airfoils.kink = 'a320';
airfoils.tip  = 'a320';

% Materials
m_up = [7.1e10, 2795, 4.8e8, 4.6e8];
m_fr = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.upper = m_up;
mat_props.lower = m_up;
mat_props.front = m_fr;
mat_props.rear  = m_up;

%% 2. Create Dummy Load Distributions
num_sections = 15;
Y = linspace(0, b/2, num_sections)';
L_total = MTOW * 9.81 * n_max; 
eta = Y ./ (b/2);
dist_shape = sqrt(1 - eta.^2);
L_avg = L_total / b; 
L = dist_shape * L_avg * (4/pi); 
M = -0.1 * L .* c_r; 

%% 3. Run Structures Module

% Example: Run with specific name and VERBOSE = true
W_wing = structures(Y, L, M, ...
                    sweep_te_k, b_k, dihedral, ...
                    b, c_r, c_k, c_t, ...
                    MTOW, ZFW, n_max, S_ref, ...
                    mat_props, spar_locs, ...
                    tank_limits, engine_data, airfoils, ...
                    myCaseName, true); % <-- Verbose ON

fprintf('Final Result returned to script: W_wing = %.2f kg\n', W_wing);