% Test EMWET sensitivity to load variations
% This script tests whether EMWET Student v1.5 responds to different load magnitudes

clear; clc;
fprintf('========== TESTING EMWET LOAD SENSITIVITY ==========\n\n');

%% Define geometry (same as main.m)
b = 34.1;           % Total wingspan (m)
b_k = 5.5;          % Kink station spanwise position (m)
c_r = 7.2;          % Root chord (m)
c_k = 4.5;          % Kink chord (m)
c_t = 2.0;          % Tip chord (m)
sweep_te_k = 3.67;  % Trailing edge sweep at kink (deg)
dihedral = 0.48;    % Dihedral angle (deg)

% Flight parameters
n_max = 2.5;
V_MO_ref = 250;
h_cr = 10000;
S_ref = 122.6;

% Weight parameters
W_AminusW = 400000;  % N
W_fuel = 150000;     % N
MTOW = 64220;  % kg
ZFW = 48930;   % kg

% Twist
twist_r = 0;
twist_k = -1.0;
twist_t = -3;

% Materials (FIXED: use different material for lower panel like B737)
m_up = [7.1e10, 2795, 4.8e8, 4.6e8];
m_lo = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];  % Stronger material for lower panel
m_fr = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.upper = m_up;
mat_props.lower = m_lo;  % FIXED: was m_up
mat_props.front = m_fr;
mat_props.rear  = m_up;

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

%% Generate simple baseline load distribution (matching B737 pattern)
% Create simple elliptical-like distribution similar to B737.load
% B737 loads: 5.99e4 N at root down to 0.27e4 N at tip
% B737 moments: -4.34e4 to +4.00e4 Nm (negative at root, positive at tip)

Y = linspace(0, b/2, 15)';  % 15 stations like B737

% Elliptical lift distribution (high at root, low at tip)
L_baseline = 6.0e4 * sqrt(1 - (2*Y/b).^2);  % N/m, elliptical shape

% Pitching moment (negative at root, transitions to positive at tip)
M_baseline = -4.5e4 + 5.0e4 * (Y/(b/2));  % Nm/m, linear transition

fprintf('Baseline loads generated (simple elliptical pattern):\n');
fprintf('  Mean L = %.2e N/m, Max L = %.2e N/m\n', mean(L_baseline), max(L_baseline));
fprintf('  Mean M = %.2e Nm/m, Max M = %.2e Nm/m\n\n', mean(abs(M_baseline)), max(abs(M_baseline)));

%% Test with different load multipliers
multipliers = [0.5, 1.0, 2.0, 5.0, 10.0];
results = zeros(length(multipliers), 1);

for i = 1:length(multipliers)
    mult = multipliers(i);
    
    % Scale loads
    L_test = L_baseline * mult;
    M_test = M_baseline * mult;
    
    fprintf('========================================\n');
    fprintf('TEST %d: Load multiplier = %.1fx\n', i, mult);
    fprintf('========================================\n');
    fprintf('  Scaled Mean L = %.2e N/m\n', mean(L_test));
    fprintf('  Scaled Mean M = %.2e Nm/m\n', mean(abs(M_test)));
    
    % Call structures with scaled loads
    caseName = sprintf('emwet_test_%dx', round(mult*10));
    W_wing_result = 9.81 * structures(Y, L_test, M_test, ...
                                      sweep_te_k, b_k, dihedral, ...
                                      b, c_r, c_k, c_t, ...
                                      MTOW, ZFW, n_max, S_ref, ...
                                      mat_props, spar_locs, ...
                                      tank_limits, engine_data, airfoils, ...
                                      caseName, false);
    
    results(i) = W_wing_result;
    fprintf('  >>> RESULT: W_wing = %.2f kg\n\n', W_wing_result/9.81);
end

%% Summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('Load Multiplier | W_wing (kg)\n');
fprintf('----------------+------------\n');
for i = 1:length(multipliers)
    fprintf('     %.1fx       |   %.2f\n', multipliers(i), results(i)/9.81);
end
fprintf('============================\n\n');

% Check if all results are identical
if all(abs(results - results(1)) < 0.01)
    fprintf('CONCLUSION: EMWET returns IDENTICAL wing weight regardless of load magnitude.\n');
    fprintf('            This confirms EMWET Student v1.5 is insensitive to loads.\n');
    fprintf('            Wing weight is driven to minimum gauge (0.8mm) in all cases.\n\n');
else
    fprintf('CONCLUSION: EMWET DOES respond to load variations.\n');
    fprintf('            The MDO coupling should work properly with sufficient load magnitudes.\n\n');
    
    % Show variation
    weight_range = (max(results) - min(results)) / 9.81;
    fprintf('            Weight variation: %.2f kg (%.1f%% of mean)\n', ...
            weight_range, 100*weight_range/mean(results)*9.81);
end
