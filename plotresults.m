% plotresults.m
% Script to plot the optimized wing planform and airfoil shape

%% Optimized Design Variables
b = 34.21;          % Wingspan (m)
c_r = 7.00;         % Root chord (m)
c_k = 3.61;         % Kink chord (m)
c_t = 1.69;         % Tip chord (m)

% Wing geometry parameters (from A320 baseline)
b_k = 6.0;          % Kink spanwise location (m) - typical for A320
sweep_te_k = 7.0;   % Trailing edge sweep at kink (deg) - typical value
dihedral = 5.0;     % Dihedral angle (deg)

% CST Parameters
CST = [0.1596, -0.0124, 0.2591, 0.2464, 0.3512, ...
    -0.0289, -0.2925, -0.0904, -0.0001, 0.0601];

%% Plot Wing Planform
plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, dihedral);
title('Optimized A320 Wing Planform - Top View', 'FontSize', 14, 'FontWeight', 'bold');

% Add optimization results annotation
annotation('textbox', [0.15, 0.02, 0.7, 0.1], ...
    'String', sprintf('Range: 6123.93 km | L/D: 26.64 | MTOW: 73659 kg'), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'EdgeColor', 'none');

%% Plot Airfoil Shape
plot_airfoil(CST);
title('Optimized Airfoil Shape (CST Parameterization)', 'FontSize', 14, 'FontWeight', 'bold');

% Add CST values annotation
annotation('textbox', [0.15, 0.02, 0.7, 0.08], ...
    'String', sprintf('Upper: [%.3f, %.3f, %.3f, %.3f, %.3f]\nLower: [%.3f, %.3f, %.3f, %.3f, %.3f]', ...
    CST(1:5), CST(6:10)), ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'EdgeColor', 'none');

%% Display Summary
fprintf('\n========== PLOTTED OPTIMIZATION RESULTS ==========\n');
fprintf('Wing Planform:\n');
fprintf('  Wingspan:    %.2f m\n', b);
fprintf('  Root Chord:  %.2f m\n', c_r);
fprintf('  Kink Chord:  %.2f m\n', c_k);
fprintf('  Tip Chord:   %.2f m\n', c_t);
fprintf('\nPerformance:\n');
fprintf('  Range:       6123.93 km\n');
fprintf('  L/D:         26.64\n');
fprintf('  MTOW:        73659 kg\n');
fprintf('===================================================\n');