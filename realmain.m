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
CST = [0.2337, 0.0796, 0.2683, 0.0887, 0.2789, 0.3811, -0.2254, -0.1634, -0.0470, -0.4771, 0.0735, 0.3255];    % CST parameters (Upper and Lower)
% Fuel weight (in Newtons) - Initial guess
spar_locs = [0.2, 0.6];  % Front and rear spar locations
tank_limits = [0, 0.85]; % Fuel tank spanwise limits
b_k = 4.36 + 3.95/2;          % Spanwise location of kink (m) [Estimated, drawing]
W_fuel = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);

% Construct design vector for optimization (ONLY c_r, c_k, c_t, W_fuel)
% x = [c_r, c_k, c_t, W_fuel]
x0 = [c_r, c_k, c_t, W_fuel];

%%  SECTION 2: DEFINE BOUNDS FOR OPTIMIZATION

% Lower bounds for design variables
c_r_lb = c_k;       % Minimum root chord (m)
c_k_lb = c_t;       % Minimum kink chord (m)
c_t_lb = 0.5;       % Minimum tip chord (m)
W_fuel_lb = 0.5 * W_fuel;  % Minimum fuel weight (N) ~5000 kg

lb = [c_r_lb, c_k_lb, c_t_lb, W_fuel_lb];

% Upper bounds for design variables
c_r_ub = 9.39;       % Maximum root chord (m)
c_k_ub = c_r;       % Maximum kink chord (m)
c_t_ub = c_k;       % Maximum tip chord (m)
W_fuel_ub = 1.5 * W_fuel;  % Maximum fuel weight (N) ~25000 kg

ub = [c_r_ub, c_k_ub, c_t_ub, W_fuel_ub];

%%  SECTION 3: RUN INITIAL EVALUATION

fprintf('\n========== INITIAL DESIGN EVALUATION ==========\n');
fprintf('Running initial design point...\n');

Range_initial = optimize(x0);

fprintf('\nInitial Design Results:\n');
fprintf('  Root Chord:      %.2f m\n', x0(1));
fprintf('  Kink Chord:      %.2f m\n', x0(2));
fprintf('  Tip Chord:       %.2f m\n', x0(3));
fprintf('  Fuel Weight:     %.2f kg\n', x0(4) / 9.81);
fprintf('  Range:           %.2f km\n', Range_initial / 1000);
fprintf('================================================\n\n');

%%  SECTION 4: SETUP AND RUN OPTIMIZER

fprintf('========== STARTING OPTIMIZATION ==========\n');

% Define objective function (negative range for maximization)
objFun = @(x) -optimize(x);

% Optimization options - tuned for speed
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...                      % Try 'sqp' to potentially reduce function evaluations
    'Display', 'iter-detailed', ...
    'MaxIterations', 100, ...
    'MaxFunctionEvaluations', 1000, ...
    'OptimalityTolerance', 1e-2, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-2, ...
    'FiniteDifferenceStepSize', 1e-6, ...  % Let MATLAB choose the step size, which can be better for scaled problems
    'FiniteDifferenceType', 'forward', ...
    'PlotFcn', {@optimplotfval, @optimplotconstrviolation, @optimplotstepsize}, ... % Added step size plot
    'UseParallel', true);                        % Use if you have the Parallel Computing Toolbox

% Run optimization (unconstrained except for bounds)
tic
[x_opt, fval, exitflag, output] = fmincon(objFun, x0, [], [], [], [], lb, ub, @constraints, options);
toc

%%  SECTION 5: DISPLAY OPTIMIZATION RESULTS

fprintf('\n========== OPTIMIZATION RESULTS ==========\n');
fprintf('Exit Flag: %d\n', exitflag);
fprintf('Iterations: %d\n', output.iterations);
fprintf('Function Evaluations: %d\n', output.funcCount);
fprintf('\nOptimal Design Variables:\n');
fprintf('  Root Chord:      %.2f m (change: %+.2f%%)\n', x_opt(1), (x_opt(1)-x0(1))/x0(1)*100);
fprintf('  Kink Chord:      %.2f m (change: %+.2f%%)\n', x_opt(2), (x_opt(2)-x0(2))/x0(2)*100);
fprintf('  Tip Chord:       %.2f m (change: %+.2f%%)\n', x_opt(3), (x_opt(3)-x0(3))/x0(3)*100);
fprintf('  Fuel Weight:     %.2f kg (change: %+.2f%%)\n', x_opt(4)/9.81, (x_opt(4)-x0(4))/x0(4)*100);
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

%%  SECTION 7: CUSTOM CONSTRAINT PLOT (two points per iteration, one per constraint)
if isfield(output, 'iterations') && isfield(output, 'firstorderopt')
    % If output structure contains iteration info, try to plot constraints
    % (fmincon does not store all iterates by default, so this is a custom plot for initial/final only)
    % Evaluate constraints at initial and optimal points
    c_hist = zeros(2,2);
    c_hist(1,:) = constraints(x0);    % Initial
    c_hist(2,:) = constraints(x_opt); % Final
    figure('Name','Constraint Values: Initial vs Optimal','NumberTitle','off');
    hold on; grid on;
    plot([1 2], c_hist(:,1), 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName','Constraint 1');
    plot([1 2], c_hist(:,2), 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName','Constraint 2');
    set(gca,'XTick',[1 2],'XTickLabel',{'Initial','Optimal'});
    ylabel('Constraint Value');
    title('Constraint Values at Initial and Optimal Points');
    legend('show');
    yline(0,'--k','Constraint Boundary');
    hold off;
end