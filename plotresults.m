% plotresults.m
% Script to compare original and optimized wing designs

clear; close all; clc;

% Setup paths
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
emwetPath = fullfile(currentDir, 'EMWET');
addpath(q3dPath);
addpath(emwetPath);

%% ORIGINAL DESIGN (from realmain.m baseline)
orig.b = 34.0;          % Wingspan (m)
orig.c_r = 7.0;         % Root chord (m)
orig.c_k = 3.7;         % Kink chord (m)
orig.c_t = 1.6;         % Tip chord (m)
orig.M_cr = 0.78;       % Cruise Mach number
orig.h_cr = 11278.4;    % Cruise altitude (m)
orig.W_fuel = 15981.29 * 9.81;  % Fuel weight (N)
orig.CST = [0.2337, 0.0796, 0.2683, 0.0887, 0.2789, 0.3811, ...
           -0.2254, -0.1634, -0.0470, -0.4771, 0.0735, 0.3255];

%% OPTIMIZED DESIGN
opt.b = 33.21;          % Wingspan (m)
opt.c_r = 7.10;         % Root chord (m)
opt.c_k = 3.80;         % Kink chord (m)
opt.c_t = 1.60;         % Tip chord (m)
opt.M_cr = 0.782;       % Cruise Mach number
opt.h_cr = 11267;       % Cruise altitude (m)
opt.W_fuel = 16275.69 * 9.81;  % Fuel weight (N)
opt.CST = [0.1731, 0.0942, 0.2982, 0.1242, 0.2848, 0.3828, ...
          -0.2023, -0.1676, -0.0620, -0.3735, 0.0652, 0.2221];

%% Common geometry parameters
b_k = 4.36 + 3.95/2;    % Kink spanwise location (m)
sweep_te_k = 0.01;      % Trailing edge sweep at kink (deg)
dihedral = 5.0;         % Dihedral angle (deg)

%% Calculate Range for both designs
fprintf('Calculating Original Design Performance...\n');
x_orig = [orig.b, orig.c_r, orig.c_k, orig.c_t, orig.M_cr, orig.h_cr, orig.W_fuel, orig.CST];
orig.Range = optimize(x_orig);
orig.Range_km = orig.Range / 1000;

fprintf('Calculating Optimized Design Performance...\n');
x_opt = [opt.b, opt.c_r, opt.c_k, opt.c_t, opt.M_cr, opt.h_cr, opt.W_fuel, opt.CST];
opt.Range = optimize(x_opt);
opt.Range_km = opt.Range / 1000;

%% Calculate wing areas
orig.S = (orig.c_r + orig.c_k) * b_k + (orig.c_k + orig.c_t) * (orig.b/2 - b_k);
opt.S = (opt.c_r + opt.c_k) * b_k + (opt.c_k + opt.c_t) * (opt.b/2 - b_k);

%% PLOT 1: Wing Planforms Comparison
figure('Name', 'Wing Planform Comparison', 'NumberTitle', 'off', 'Position', [100 100 900 600]);
hold on; grid on; axis equal;

% Function to compute wing outline
compute_wing = @(b, c_r, c_k, c_t) deal(...
    [0, c_r + b_k*tand(sweep_te_k) - c_k, 0 + (b/2)*((c_r + b_k*tand(sweep_te_k) - c_k)/b_k), ...
     0 + (b/2)*((c_r + b_k*tand(sweep_te_k) - c_k)/b_k) + c_t, c_r + b_k*tand(sweep_te_k), c_r, 0], ...
    [0, b_k, b/2, b/2, b_k, 0, 0]);

% Original wing coordinates
x_le_r_orig = 0; x_te_r_orig = orig.c_r;
x_te_k_orig = x_te_r_orig + b_k * tand(sweep_te_k);
x_le_k_orig = x_te_k_orig - orig.c_k;
tan_sweep_le_orig = (x_le_k_orig - x_le_r_orig) / b_k;
x_le_t_orig = x_le_r_orig + (orig.b/2) * tan_sweep_le_orig;
x_te_t_orig = x_le_t_orig + orig.c_t;
x_orig_wing = [x_le_r_orig, x_le_k_orig, x_le_t_orig, x_te_t_orig, x_te_k_orig, x_te_r_orig, x_le_r_orig];
y_orig_wing = [0, b_k, orig.b/2, orig.b/2, b_k, 0, 0];

% Optimized wing coordinates
x_le_r_opt = 0; x_te_r_opt = opt.c_r;
x_te_k_opt = x_te_r_opt + b_k * tand(sweep_te_k);
x_le_k_opt = x_te_k_opt - opt.c_k;
tan_sweep_le_opt = (x_le_k_opt - x_le_r_opt) / b_k;
x_le_t_opt = x_le_r_opt + (opt.b/2) * tan_sweep_le_opt;
x_te_t_opt = x_le_t_opt + opt.c_t;
x_opt_wing = [x_le_r_opt, x_le_k_opt, x_le_t_opt, x_te_t_opt, x_te_k_opt, x_te_r_opt, x_le_r_opt];
y_opt_wing = [0, b_k, opt.b/2, opt.b/2, b_k, 0, 0];

% Plot both wings
fill(x_orig_wing, y_orig_wing, [0.9 0.9 0.9], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2, 'FaceAlpha', 0.5, 'DisplayName', 'Original');
fill(x_opt_wing, y_opt_wing, [0.7 0.85 1], 'EdgeColor', [0 0.4 0.7], 'LineWidth', 2, 'FaceAlpha', 0.7, 'DisplayName', 'Optimized');

xlabel('Chordwise Position (m)', 'FontSize', 12);
ylabel('Spanwise Position (m)', 'FontSize', 12);
title('Wing Planform Comparison: Original vs Optimized', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);

% Add annotation box
annotation('textbox', [0.15, 0.75, 0.35, 0.15], ...
    'String', sprintf('Original:\n  b=%.2fm, S=%.1fm²\n  Range=%.0fkm', orig.b, orig.S, orig.Range_km), ...
    'FontSize', 10, 'BackgroundColor', [0.95 0.95 0.95], 'EdgeColor', [0.5 0.5 0.5]);
annotation('textbox', [0.55, 0.75, 0.35, 0.15], ...
    'String', sprintf('Optimized:\n  b=%.2fm, S=%.1fm²\n  Range=%.0fkm (+%.1f%%)', opt.b, opt.S, opt.Range_km, (opt.Range-orig.Range)/orig.Range*100), ...
    'FontSize', 10, 'BackgroundColor', [0.85 0.92 1], 'EdgeColor', [0 0.4 0.7]);

%% PLOT 2: Airfoil Shape Comparison
figure('Name', 'Airfoil Shape Comparison', 'NumberTitle', 'off', 'Position', [150 150 900 500]);
hold on; grid on; axis equal;

% Generate airfoil coordinates using CST
x_af = linspace(0, 1, 200);

% CST shape function
cst_shape = @(x, A) sqrt(x) .* (1-x) .* sum(cell2mat(arrayfun(@(i) A(i) .* nchoosek(length(A)-1, i-1) .* x.^(i-1) .* (1-x).^(length(A)-i), 1:length(A), 'UniformOutput', false)'), 1);

% Original airfoil
y_upper_orig = cst_shape(x_af, orig.CST(1:6));
y_lower_orig = cst_shape(x_af, orig.CST(7:12));

% Optimized airfoil
y_upper_opt = cst_shape(x_af, opt.CST(1:6));
y_lower_opt = cst_shape(x_af, opt.CST(7:12));

% Plot airfoils
plot(x_af, y_upper_orig, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', 'Original Upper');
plot(x_af, y_lower_orig, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', 'Original Lower');
plot(x_af, y_upper_opt, '-', 'Color', [0 0.4 0.7], 'LineWidth', 2, 'DisplayName', 'Optimized Upper');
plot(x_af, y_lower_opt, '-', 'Color', [0 0.4 0.7], 'LineWidth', 2, 'DisplayName', 'Optimized Lower');

% Fill areas for visualization
fill([x_af, fliplr(x_af)], [y_upper_orig, fliplr(y_lower_orig)], [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([x_af, fliplr(x_af)], [y_upper_opt, fliplr(y_lower_opt)], [0.7 0.85 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('x/c', 'FontSize', 12);
ylabel('y/c', 'FontSize', 12);
title('Airfoil Shape Comparison: Original vs Optimized', 'FontSize', 14, 'FontWeight', 'bold');
legend('Original Upper', 'Original Lower', 'Optimized Upper', 'Optimized Lower', 'Location', 'best', 'FontSize', 10);
xlim([-0.05 1.05]);

%% Check Constraints for both designs
fprintf('\nChecking constraints...\n');
[c_orig, ceq_orig] = constraints(x_orig);
[c_opt, ceq_opt] = constraints(x_opt);

%% Display Summary
fprintf('\n==================== DESIGN COMPARISON ====================\n');
fprintf('                        Original      Optimized      Change\n');
fprintf('------------------------------------------------------------\n');
fprintf('Wingspan (m):           %8.2f      %8.2f      %+6.2f%%\n', orig.b, opt.b, (opt.b-orig.b)/orig.b*100);
fprintf('Root Chord (m):         %8.2f      %8.2f      %+6.2f%%\n', orig.c_r, opt.c_r, (opt.c_r-orig.c_r)/orig.c_r*100);
fprintf('Kink Chord (m):         %8.2f      %8.2f      %+6.2f%%\n', orig.c_k, opt.c_k, (opt.c_k-orig.c_k)/orig.c_k*100);
fprintf('Tip Chord (m):          %8.2f      %8.2f      %+6.2f%%\n', orig.c_t, opt.c_t, (opt.c_t-orig.c_t)/orig.c_t*100);
fprintf('Wing Area (m²):         %8.2f      %8.2f      %+6.2f%%\n', orig.S, opt.S, (opt.S-orig.S)/orig.S*100);
fprintf('Cruise Mach:            %8.3f      %8.3f      %+6.2f%%\n', orig.M_cr, opt.M_cr, (opt.M_cr-orig.M_cr)/orig.M_cr*100);
fprintf('Cruise Alt (m):         %8.0f      %8.0f      %+6.2f%%\n', orig.h_cr, opt.h_cr, (opt.h_cr-orig.h_cr)/orig.h_cr*100);
fprintf('Fuel Weight (kg):       %8.2f      %8.2f      %+6.2f%%\n', orig.W_fuel/9.81, opt.W_fuel/9.81, (opt.W_fuel-orig.W_fuel)/orig.W_fuel*100);
fprintf('------------------------------------------------------------\n');
fprintf('RANGE (km):             %8.2f      %8.2f      %+6.2f%%\n', orig.Range_km, opt.Range_km, (opt.Range-orig.Range)/orig.Range*100);
fprintf('============================================================\n');

%% Display Constraint Verification
fprintf('\n================ CONSTRAINT VERIFICATION ===================\n');
fprintf('Inequality Constraints (c <= 0 means SATISFIED):\n');
fprintf('------------------------------------------------------------\n');
fprintf('                        Original      Optimized\n');
fprintf('c1 (Wing Loading):    %10.2f %s  %10.2f %s\n', ...
    c_orig(1), status_str(c_orig(1) <= 0), c_opt(1), status_str(c_opt(1) <= 0));
fprintf('c2 (Fuel Capacity):   %10.2f %s  %10.2f %s\n', ...
    c_orig(2), status_str(c_orig(2) <= 0), c_opt(2), status_str(c_opt(2) <= 0));
fprintf('------------------------------------------------------------\n');
if isempty(ceq_orig) && isempty(ceq_opt)
    fprintf('Equality Constraints: None\n');
else
    fprintf('Equality Constraints (ceq = 0 means SATISFIED):\n');
    for i = 1:length(ceq_orig)
        fprintf('ceq%d:                 %10.4f %s  %10.4f %s\n', i, ...
            ceq_orig(i), status_str(abs(ceq_orig(i)) < 1e-3), ...
            ceq_opt(i), status_str(abs(ceq_opt(i)) < 1e-3));
    end
end
fprintf('============================================================\n');

%% Helper function for constraint status
function s = status_str(satisfied)
    if satisfied
        s = '[OK]';
    else
        s = '[X] ';
    end
end