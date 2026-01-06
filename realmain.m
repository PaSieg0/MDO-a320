%% Main MDO Optimization Script for A320 Wing Design
% This script sets up and runs the optimization using the optimize.m function
clear all;
close all;
clc;

% Setup paths
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
emwetPath = fullfile(currentDir, 'EMWET');
addpath(q3dPath);
addpath(emwetPath);

%%  SECTION 1: DEFINE DESIGN VECTOR INITIAL VALUES

% Planform dimensions (in meters) - Source: Assignment specification
b = 34.0;          % Total wingspan (m) [Drawing]
c_r = 7.0;          % Root chord (m) [Drawing]
c_k = 3.7;         % Kink chord (m) [Drawing]
c_t = 1.6;          % Tip chord (m) [Drawing]

% Flight parameters - Source: Assignment specification
M_cr = 0.78;        % Cruise Mach number [-] [Specification]
h_cr = 11278.4;     % Cruise altitude (m) = 37,000 ft [Specification: 37000 ft]

% Airfoil shape (CST parameters) - Source: Baseline airfoil from assignment
CST = [0.1800 0.2800 0.2400 0.2200 0.2400 -0.1000 -0.1800 -0.1200 -0.0300 0.0600];

% Fuel weight (in Newtons) - Initial guess
W_fuel = 16275.44 * 9.81;  % Fuel weight (N) [Estimated for medium range]

% Construct design vector for optimization
% x = [b, c_r, c_k, c_t, M_cr, h_cr, W_fuel, CST(1:10)]
x0 = [b, c_r, c_k, c_t, M_cr, h_cr, W_fuel, CST];

%%  SECTION 2: DEFINE BOUNDS FOR OPTIMIZATION

% Lower bounds for design variables
b_lb = 28.0;        % Minimum wingspan (m)
c_r_lb = 5.0;       % Minimum root chord (m)
c_k_lb = 2.5;       % Minimum kink chord (m)
c_t_lb = 1.0;       % Minimum tip chord (m)
M_cr_lb = 0.70;     % Minimum cruise Mach
h_cr_lb = 9000;     % Minimum cruise altitude (m)
W_fuel_lb = 5000 * 9.81;  % Minimum fuel weight (N) ~5000 kg
CST_lb = [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5];

lb = [b_lb, c_r_lb, c_k_lb, c_t_lb, M_cr_lb, h_cr_lb, W_fuel_lb, CST_lb];

% Upper bounds for design variables
b_ub = 40.0;        % Maximum wingspan (m)
c_r_ub = 9.0;       % Maximum root chord (m)
c_k_ub = 5.0;       % Maximum kink chord (m)
c_t_ub = 2.5;       % Maximum tip chord (m)
M_cr_ub = 0.82;     % Maximum cruise Mach
h_cr_ub = 13000;    % Maximum cruise altitude (m)
W_fuel_ub = 25000 * 9.81;  % Maximum fuel weight (N) ~25000 kg
CST_ub = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

ub = [b_ub, c_r_ub, c_k_ub, c_t_ub, M_cr_ub, h_cr_ub, W_fuel_ub, CST_ub];

%%  SECTION 3: RUN INITIAL EVALUATION

fprintf('\n========== INITIAL DESIGN EVALUATION ==========\n');
fprintf('Running initial design point...\n');

Range_initial = optimize(x0);

fprintf('\nInitial Design Results:\n');
fprintf('  Wingspan:        %.2f m\n', x0(1));
fprintf('  Root Chord:      %.2f m\n', x0(2));
fprintf('  Kink Chord:      %.2f m\n', x0(3));
fprintf('  Tip Chord:       %.2f m\n', x0(4));
fprintf('  Cruise Mach:     %.2f\n', x0(5));
fprintf('  Cruise Altitude: %.0f m\n', x0(6));
fprintf('  Fuel Weight:     %.2f kg\n', x0(7) / 9.81);
fprintf('  Range:           %.2f km\n', Range_initial / 1000);
fprintf('================================================\n\n');

%%  SECTION 4: SETUP AND RUN OPTIMIZER

fprintf('========== STARTING OPTIMIZATION ==========\n');

% Define objective function (negative range for maximization)
objFun = @(x) -optimize(x);

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 50, ...
    'MaxFunctionEvaluations', 500, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-6, ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'PlotFcn', @optimplotfval);

% Run optimization (unconstrained except for bounds)
[x_opt, fval, exitflag, output] = fmincon(objFun, x0, [], [], [], [], lb, ub, @constraints, options);

%%  SECTION 5: DISPLAY OPTIMIZATION RESULTS

fprintf('\n========== OPTIMIZATION RESULTS ==========\n');
fprintf('Exit Flag: %d\n', exitflag);
fprintf('Iterations: %d\n', output.iterations);
fprintf('Function Evaluations: %d\n', output.funcCount);
fprintf('\nOptimal Design Variables:\n');
fprintf('  Wingspan:        %.2f m (change: %+.2f%%)\n', x_opt(1), (x_opt(1)-x0(1))/x0(1)*100);
fprintf('  Root Chord:      %.2f m (change: %+.2f%%)\n', x_opt(2), (x_opt(2)-x0(2))/x0(2)*100);
fprintf('  Kink Chord:      %.2f m (change: %+.2f%%)\n', x_opt(3), (x_opt(3)-x0(3))/x0(3)*100);
fprintf('  Tip Chord:       %.2f m (change: %+.2f%%)\n', x_opt(4), (x_opt(4)-x0(4))/x0(4)*100);
fprintf('  Cruise Mach:     %.3f (change: %+.2f%%)\n', x_opt(5), (x_opt(5)-x0(5))/x0(5)*100);
fprintf('  Cruise Altitude: %.0f m (change: %+.2f%%)\n', x_opt(6), (x_opt(6)-x0(6))/x0(6)*100);
fprintf('  Fuel Weight:     %.2f kg (change: %+.2f%%)\n', x_opt(7)/9.81, (x_opt(7)-x0(7))/x0(7)*100);
fprintf('\nCST Parameters (Upper Surface):\n');
fprintf('  [%.4f, %.4f, %.4f, %.4f, %.4f]\n', x_opt(8:12));
fprintf('CST Parameters (Lower Surface):\n');
fprintf('  [%.4f, %.4f, %.4f, %.4f, %.4f]\n', x_opt(13:17));
fprintf('\nPerformance:\n');
fprintf('  Initial Range: %.2f km\n', Range_initial / 1000);
fprintf('  Optimal Range: %.2f km\n', -fval / 1000);
fprintf('  Improvement:   %+.2f%%\n', (-fval - Range_initial) / Range_initial * 100);
fprintf('==========================================\n\n');

%%  SECTION 6: FINAL VERIFICATION

fprintf('========== FINAL VERIFICATION ==========\n');
Range_final = optimize(x_opt);
fprintf('Verified Optimal Range: %.2f km\n', Range_final / 1000);
fprintf('========================================\n');