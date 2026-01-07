% plotresults.m
% Script to plot the optimized wing planform and airfoil shape

%% Optimized Design Variables
b = 33.21;          % Wingspan (m)
c_r = 7.10;         % Root chord (m)
c_k = 3.80;         % Kink chord (m)
c_t = 1.60;         % Tip chord (m)

% Flight conditions
M_cr = 0.782;       % Cruise Mach number
h_cr = 11267;       % Cruise altitude (m)
W_fuel = 16275.69;  % Fuel weight (kg)

% Wing geometry parameters (from A320 baseline)
b_k = 6.0;          % Kink spanwise location (m) - typical for A320
sweep_te_k = 7.0;   % Trailing edge sweep at kink (deg) - typical value
dihedral = 5.0;     % Dihedral angle (deg)

% CST Parameters (6 upper + 6 lower)
CST = [0.1731, 0.0942, 0.2982, 0.1242, 0.2848, 0.3828, ...
      -0.2023, -0.1676, -0.0620, -0.3735, 0.0652, 0.2221];

% Performance results
Range_km = 6363.33;
L_D = 27.86;
MTOW_kg = 73971;
W_wing_kg = 5912;

%% Plot Wing Planform
plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, dihedral);
title('Optimized A320 Wing Planform - Top View', 'FontSize', 14, 'FontWeight', 'bold');

% Add optimization results annotation
annotation('textbox', [0.15, 0.02, 0.7, 0.1], ...
    'String', sprintf('Range: %.2f km | L/D: %.2f | MTOW: %d kg | W_{wing}: %d kg', ...
    Range_km, L_D, MTOW_kg, W_wing_kg), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'EdgeColor', 'none');

%% Plot Airfoil Shape
plot_airfoil(CST);
title('Optimized Airfoil Shape (CST Parameterization)', 'FontSize', 14, 'FontWeight', 'bold');

% Add CST values annotation
annotation('textbox', [0.15, 0.02, 0.7, 0.08], ...
    'String', sprintf('Upper: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\nLower: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]', ...
    CST(1:6), CST(7:12)), ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'EdgeColor', 'none');

%% Display Summary
fprintf('\n========== PLOTTED OPTIMIZATION RESULTS ==========\n');
fprintf('Wing Planform:\n');
fprintf('  Wingspan:    %.2f m\n', b);
fprintf('  Root Chord:  %.2f m\n', c_r);
fprintf('  Kink Chord:  %.2f m\n', c_k);
fprintf('  Tip Chord:   %.2f m\n', c_t);
fprintf('Flight Conditions:\n');
fprintf('  Cruise Mach: %.3f\n', M_cr);
fprintf('  Cruise Alt:  %.0f m\n', h_cr);
fprintf('Performance:\n');
fprintf('  Range:       %.2f km\n', Range_km);
fprintf('  L/D:         %.2f\n', L_D);
fprintf('  MTOW:        %d kg\n', MTOW_kg);
fprintf('  Wing Weight: %d kg\n', W_wing_kg);
fprintf('===================================================\n');