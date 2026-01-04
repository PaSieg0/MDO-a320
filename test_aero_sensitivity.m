%% Test script for aerodynamic discipline (aero.m)
% Tests sensitivity of Cl and Cd to various input parameters
% Isolates aero discipline from loads, structures, and performance

clear all;
close all;
clc;

% Setup paths
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
addpath(q3dPath);

fprintf('========== AERODYNAMIC DISCIPLINE SENSITIVITY TEST ==========\n');
fprintf('NOTE: Testing near main.m baseline conditions (M=0.78, h=11278m)\n');
fprintf('Goal: Identify where Q3D transonic analysis diverges\n\n');

%% Baseline parameters (EXACT from main.m)

% Geometry - EXACT from main.m
b = 34.0;           % Wingspan (m)
c_r = 7.0;          % Root chord (m)
c_k = 3.7;          % Kink chord (m)
c_t = 1.6;          % Tip chord (m)
b_k = 4.36 + 3.95/2;  % Kink location (m)

% Wing angles - EXACT from main.m
sweep_te_k = 0.01;  % TE sweep (deg)
dihedral = 5;       % Dihedral (deg)
twist_r = 0;        % Root twist (deg)
twist_k = -1.0;     % Kink twist (deg)
twist_t = -3;       % Tip twist (deg)

% Airfoil - EXACT from main.m
CST = [0.1800 0.2800 0.2400 0.2200 0.2400 -0.1000 -0.1800 -0.1200 -0.0300 0.0600];

% Flight condition - EXACT from main.m (the problematic one)
M_cr = 0.78;        % Cruise Mach (DIVERGES in main.m)
h_cr = 11278.4;     % Cruise altitude (37,000 ft)

% Weights - EXACT from main.m
W_AminusW = 400000; % Aircraft minus wing (N)
W_wing = 50000;     % Wing weight (N)
W_fuel = 150000;    % Fuel weight (N)

% Constants
gamma = 1.4;
R = 287.058;

%% Test 1: Baseline evaluation (EXACT main.m conditions)

fprintf('TEST 1: Baseline Configuration (EXACT from main.m)\n');
fprintf('---------------------------------------------------\n');

[Cl_base, Cd_base] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                           M_cr, W_AminusW, h_cr,...
                           b, c_r, c_k, c_t, CST, W_wing, W_fuel);

fprintf('Mach: %.2f, Altitude: %.0f m\n', M_cr, h_cr);
if isnan(Cd_base)
    fprintf('*** DIVERGED *** Cl = %.4f, Cd = NaN\n\n', Cl_base);
else
    fprintf('Cl = %.4f, Cd = %.4f, L/D = %.2f\n\n', Cl_base, Cd_base, Cl_base/Cd_base);
end

%% Test 2: Mach number variations (narrow range around 0.78)

fprintf('TEST 2: Mach Number Sensitivity (around M=0.78)\n');
fprintf('------------------------------------------------\n');

M_test = [0.72, 0.74, 0.76, 0.77, 0.78, 0.79, 0.80];  % Narrow range to find divergence point
n_M = length(M_test);
Cl_M = zeros(1, n_M);
Cd_M = zeros(1, n_M);
LD_M = zeros(1, n_M);

for i = 1:n_M
    [Cl_M(i), Cd_M(i)] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                               M_test(i), W_AminusW, h_cr,...
                               b, c_r, c_k, c_t, CST, W_wing, W_fuel);
    if isnan(Cd_M(i))
        fprintf('M = %.2f: DIVERGED (Cd = NaN)\n', M_test(i));
        LD_M(i) = NaN;
    else
        LD_M(i) = Cl_M(i) / Cd_M(i);
        fprintf('M = %.2f: Cl = %.4f, Cd = %.4f, L/D = %.2f\n', M_test(i), Cl_M(i), Cd_M(i), LD_M(i));
    end
end
fprintf('\n');

%% Test 3: Altitude variations (around 11278m)

fprintf('TEST 3: Altitude Sensitivity (around h=11278m)\n');
fprintf('-----------------------------------------------\n');

h_test = [10000, 10500, 11000, 11278.4, 11500, 12000]; % Around cruise altitude
n_h = length(h_test);
Cl_h = zeros(1, n_h);
Cd_h = zeros(1, n_h);
LD_h = zeros(1, n_h);

for i = 1:n_h
    [Cl_h(i), Cd_h(i)] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                               M_cr, W_AminusW, h_test(i),...
                               b, c_r, c_k, c_t, CST, W_wing, W_fuel);
    if isnan(Cd_h(i))
        fprintf('h = %.0f m: DIVERGED (Cd = NaN)\n', h_test(i));
        LD_h(i) = NaN;
    else
        LD_h(i) = Cl_h(i) / Cd_h(i);
        fprintf('h = %.0f m: Cl = %.4f, Cd = %.4f, L/D = %.2f\n', h_test(i), Cl_h(i), Cd_h(i), LD_h(i));
    end
end
fprintf('\n');

%% Test 4: Combined Mach-Altitude matrix (find safe operating region)

fprintf('TEST 4: Combined Mach-Altitude Matrix\n');
fprintf('--------------------------------------\n');

M_matrix = [0.74, 0.76, 0.78];
h_matrix = [10000, 11278.4, 12000];
n_M_mat = length(M_matrix);
n_h_mat = length(h_matrix);

fprintf('      ');
for j = 1:n_h_mat
    fprintf('h=%.0fm   ', h_matrix(j));
end
fprintf('\n');

for i = 1:n_M_mat
    fprintf('M=%.2f: ', M_matrix(i));
    for j = 1:n_h_mat
        [Cl_temp, Cd_temp] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                                   M_matrix(i), W_AminusW, h_matrix(j),...
                                   b, c_r, c_k, c_t, CST, W_wing, W_fuel);
        if isnan(Cd_temp)
            fprintf('DIVERG  ');
        else
            fprintf('L/D=%.1f  ', Cl_temp/Cd_temp);
        end
    end
    fprintf('\n');
end
fprintf('\n');

%% Convergence Summary

fprintf('========== CONVERGENCE SUMMARY ==========\n');
fprintf('Test 1 (Baseline M=%.2f, h=%.0fm):  %s\n', M_cr, h_cr, iif(isnan(Cd_base), '*** DIVERGED ***', 'CONVERGED'));
fprintf('Test 2 (Mach variations):            %d/%d converged\n', sum(~isnan(Cd_M)), n_M);
fprintf('Test 3 (Altitude variations):        %d/%d converged\n', sum(~isnan(Cd_h)), n_h);
fprintf('\n');

% Identify convergence boundaries
M_converged = M_test(~isnan(Cd_M));
M_diverged = M_test(isnan(Cd_M));
h_converged = h_test(~isnan(Cd_h));
h_diverged = h_test(isnan(Cd_h));

if ~isempty(M_converged) && ~isempty(M_diverged)
    fprintf('Mach convergence boundary: between M=%.2f and M=%.2f\n', max(M_converged), min(M_diverged));
end
if ~isempty(h_converged) && ~isempty(h_diverged)
    fprintf('Altitude convergence boundary: between h=%.0f m and h=%.0f m\n', max(h_converged), min(h_diverged));
end

fprintf('=========================================\n\n');

if isnan(Cd_base)
    fprintf('*** CRITICAL: Baseline condition (M=%.2f, h=%.0f m) DIVERGES ***\n', M_cr, h_cr);
    fprintf('Root cause likely:\n');
    fprintf('  1. Transonic shocks on airfoil (M > M_crit)\n');
    fprintf('  2. CST airfoil not suitable for M=0.78 operation\n');
    fprintf('  3. Q3D transonic solver convergence limits\n');
    fprintf('Recommendation: Use lower Mach (M<0.76) or different airfoil\n\n');
end

%% Plotting results (only if some cases converged)

if sum(~isnan(Cd_M)) > 1 || sum(~isnan(Cd_h)) > 1
    figure('Position', [100 100 1200 500]);
    
    % Mach sweep
    subplot(1,2,1);
    yyaxis left
    plot(M_test, Cl_M, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_l', 'FontSize', 12);
    ylim([0 max(Cl_M(~isnan(Cl_M)))*1.2]);
    yyaxis right
    plot(M_test, Cd_M, 's-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_d', 'FontSize', 12);
    xlabel('Mach Number', 'FontSize', 12);
    title('Mach Sensitivity (Converged Points Only)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Altitude sweep
    subplot(1,2,2);
    yyaxis left
    plot(h_test/1000, Cl_h, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_l', 'FontSize', 12);
    ylim([0 max(Cl_h(~isnan(Cl_h)))*1.2]);
    yyaxis right
    plot(h_test/1000, Cd_h, 's-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_d', 'FontSize', 12);
    xlabel('Altitude (km)', 'FontSize', 12);
    title('Altitude Sensitivity (Converged Points Only)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    sgtitle(sprintf('Aerodynamic Sensitivity Near M=%.2f, h=%.0fm', M_cr, h_cr), 'FontSize', 16, 'FontWeight', 'bold');
else
    fprintf('WARNING: Insufficient converged cases for plotting\n');
end

fprintf('\n========== TEST COMPLETE ==========\n');
fprintf('------------------------------\n');

twist_t_test = [-5, -4, -3, -2, -1, 0];
n_twist = length(twist_t_test);
Cl_twist = zeros(1, n_twist);
Cd_twist = zeros(1, n_twist);
LD_twist = zeros(1, n_twist);

for i = 1:n_twist
    [Cl_twist(i), Cd_twist(i)] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t_test(i),... 
                                       M_cr, W_AminusW, h_cr,...
                                       b, c_r, c_k, c_t, CST, W_wing, W_fuel);
    if isnan(Cd_twist(i))
        fprintf('Twist_tip = %+.1f deg: DIVERGED (Cd = NaN)\n', twist_t_test(i));
        LD_twist(i) = NaN;
    else
        LD_twist(i) = Cl_twist(i) / Cd_twist(i);
        fprintf('Twist_tip = %+.1f deg: Cl = %.4f, Cd = %.4f, L/D = %.2f\n', ...
                twist_t_test(i), Cl_twist(i), Cd_twist(i), LD_twist(i));
    end
end
fprintf('\n');

%% Plotting results (only if some cases converged)

if sum(~isnan(Cd_M)) > 1 || sum(~isnan(Cd_h)) > 1
    figure('Position', [100 100 1200 500]);
    
    % Mach sweep
    subplot(1,2,1);
    yyaxis left
    plot(M_test, Cl_M, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_l', 'FontSize', 12);
    ylim([0 max(Cl_M(~isnan(Cl_M)))*1.2]);
    yyaxis right
    plot(M_test, Cd_M, 's-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_d', 'FontSize', 12);
    xlabel('Mach Number', 'FontSize', 12);
    title('Mach Sensitivity (Converged Points Only)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    % Altitude sweep
    subplot(1,2,2);
    yyaxis left
    plot(h_test/1000, Cl_h, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_l', 'FontSize', 12);
    ylim([0 max(Cl_h(~isnan(Cl_h)))*1.2]);
    yyaxis right
    plot(h_test/1000, Cd_h, 's-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('C_d', 'FontSize', 12);
    xlabel('Altitude (km)', 'FontSize', 12);
    title('Altitude Sensitivity (Converged Points Only)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    sgtitle(sprintf('Aerodynamic Sensitivity Near M=%.2f, h=%.0fm', M_cr, h_cr), 'FontSize', 16, 'FontWeight', 'bold');
else
    fprintf('WARNING: Insufficient converged cases for plotting\n');
end

fprintf('\n========== TEST COMPLETE ==========\n');

%% Helper function
function result = iif(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end
