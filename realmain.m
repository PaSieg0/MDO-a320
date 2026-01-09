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

global W_wing; % <<< ADD THIS LINE
W_wing = 6344*9.81; % Initialize the guess ONCE here.

%%  SECTION 1: DEFINE DESIGN VECTOR INITIAL VALUES (in physical units)

% Planform dimensions (in meters) - Source: Assignment specification
b = 34.0;          % Total wingspan (m) [Drawing]
c_r = 7.0;          % Root chord (m) [Drawing]
c_k = 3.7;         % Kink chord (m) [Drawing]
c_t = 1.6;          % Tip chord (m) [Drawing]

% Flight parameters - Source: Assignment specification
M_cr = 0.78;        % Cruise Mach number [-] [Specification]
h_cr = 11278.4;     % Cruise altitude (m) = 37,000 ft [Specification: 37000 ft]

% Airfoil shape (CST parameters) - Source: Baseline airfoil from assignment
CST = [0.2337, 0.0796, 0.2683, 0.0887, 0.2789, 0.3811, -0.2254, -0.1634, -0.0470, -0.4771, 0.0735, 0.3255];    % CST parameters (Upper and Lower)
% Fuel weight (in Newtons) - Initial guess
spar_locs = [0.2, 0.6];  % Front and rear spar locations
tank_limits = [0, 0.85]; % Fuel tank spanwise limits
b_k = 4.36 + 3.95/2;          % Spanwise location of kink (m) [Estimated, drawing]
W_fuel = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);

% Construct design vector for optimization (in physical units)
x0 = [b, c_r, c_k, c_t, M_cr, h_cr, W_fuel, CST];

%%  SECTION 2: DEFINE BOUNDS FOR OPTIMIZATION (in physical units)

% Lower bounds for design variables
b_lb = 28.0;       % Minimum wingspan (m)
c_r_lb = c_k;       % Minimum root chord (m)
c_k_lb = c_t;       % Minimum kink chord (m)
c_t_lb = 0.5;       % Minimum tip chord (m)
M_cr_lb = 0.9*M_cr;   % Minimum cruise Mach number
h_cr_lb = 0.9*h_cr;    % Minimum cruise altitude (m)
W_fuel_lb = 0.5 * W_fuel;  % Minimum fuel weight (N) ~5000 kg
CST_lb = -0.5 * ones(1,12); % Lower bounds for CST parameters

lb = [b_lb, c_r_lb, c_k_lb, c_t_lb, M_cr_lb, h_cr_lb, W_fuel_lb, CST_lb];

% Upper bounds for design variables
b_ub = 36.0;       % Maximum wingspan (m)
c_r_ub = 9.39;       % Maximum root chord (m)
c_k_ub = c_r;       % Maximum kink chord (m)
c_t_ub = c_k;       % Maximum tip chord (m)
M_cr_ub = 0.82;   % Maximum cruise Mach number
h_cr_ub = 1.1*h_cr;    % Maximum cruise altitude (m)
W_fuel_ub = 1.5 * W_fuel;  % Maximum fuel weight (N) ~25000 kg
CST_ub = 0.5 * ones(1,12); % Upper bounds for CST parameters

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
fprintf('  Cruise Mach:     %.3f\n', x0(5));
fprintf('  Cruise Altitude: %.0f m\n', x0(6));
fprintf('  Fuel Weight:     %.2f kg\n', x0(7) / 9.81);
fprintf('  CST Parameters:  '); fprintf('%.4f ', x0(8:end)); fprintf('\n');
fprintf('  Range:           %.2f km\n', Range_initial / 1000);
fprintf('================================================\n\n');

%%  SECTION 4: SETUP AND RUN OPTIMIZER

fprintf('========== STARTING OPTIMIZATION ==========\n');

% <<< MODIFICATION: NORMALIZATION SETUP >>>
% Create functions to map variables between physical space and normalized [-1, 1] space.
% Denormalize: Maps a normalized vector x_norm from [-1, 1] to the physical range [lb, ub]
denormalize = @(x_norm) lb + (ub - lb) .* (x_norm + 1) / 2;
% Normalize: Maps a physical vector x_phys from [lb, ub] to the normalized range [-1, 1]
normalize = @(x_phys) 2 * (x_phys - lb) ./ (ub - lb) - 1;

% <<< MODIFICATION: NORMALIZE INITIAL GUESS AND SET NORMALIZED BOUNDS >>>
x0_norm = normalize(x0);
lb_norm = -ones(size(x0));
ub_norm = ones(size(x0));

% <<< MODIFICATION: CREATE WRAPPER FUNCTIONS FOR OPTIMIZER >>>
% These wrappers take a normalized vector, denormalize it, and then call your original functions.
objFun_norm = @(x_norm) - log(optimize(denormalize(x_norm)));
constrFun_norm = @(x_norm) constraints(denormalize(x_norm));

% Optimization options - Now that variables are scaled, these should be more effective.
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ... % 'sqp' is often a good choice for scaled problems
    'Display', 'iter-detailed', ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 1000, ...
    'OptimalityTolerance', 1e-3, ... % Can be a bit tighter now
    'StepTolerance', 1e-4, ...       % Tighter step tolerance
    'ConstraintTolerance', 1e-6, ... % Much tighter constraint tolerance
    'FiniteDifferenceStepSize', 1e-4, ... % A single small value is now effective
    'FiniteDifferenceType', 'forward', ...
    'PlotFcn', {@optimplotfval, @optimplotconstrviolation, @optimplotstepsize}, ...
    'UseParallel', false);

% <<< MODIFICATION: RUN OPTIMIZER WITH NORMALIZED VALUES >>>
tic
[x_opt_norm, fval, exitflag, output] = fmincon(objFun_norm, x0_norm, [], [], [], [], lb_norm, ub_norm, constrFun_norm, options);
toc

% <<< MODIFICATION: DENORMALIZE THE FINAL RESULT FOR ANALYSIS >>>
x_opt = denormalize(x_opt_norm);

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
fprintf('  CST Parameters:  '); fprintf('%.4f ', x_opt(8:end)); fprintf('\n');
fprintf('\nPerformance:\n');
fprintf('  Initial Range: %.2f km\n', Range_initial / 1000);
fprintf('  Optimal Range: %.2f km\n', -fval / 1000);
fprintf('  Improvement:   %+.2f%%\n', (-fval - Range_initial) / Range_initial * 100);
fprintf('==========================================\n\n');

%%  SECTION 6: FINAL VERIFICATION

fprintf('========== FINAL VERIFICATION ==========\n');
Range_final = optimize(x_opt);
fprintf('Verified Optimal Range: %.2f km\n', Range_final / 1000);
% Check constraints at optimal solution
[c_final, ceq_final] = constraints(x_opt);
fprintf('\nConstraint Verification:\n');
fprintf('  Inequality constraints (c <= 0):\n');
for i = 1:length(c_final)
    if c_final(i) <= 0
        fprintf('    c(%d) = %.4f [SATISFIED]\n', i, c_final(i));
    else
        fprintf('    c(%d) = %.4f [VIOLATED]\n', i, c_final(i));
    end
end
fprintf('  Equality constraints (ceq = 0):\n');
for i = 1:length(ceq_final)
    if abs(ceq_final(i)) < 1e-3
        fprintf('    ceq(%d) = %.4f [SATISFIED]\n', i, ceq_final(i));
    else
        fprintf('    ceq(%d) = %.4f [VIOLATED]\n', i, ceq_final(i));
    end
end
fprintf('========================================\n');
