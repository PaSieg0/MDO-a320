%% MDO Assignment - Part 2: Deliverables Generator (Final Master)
% Generates tables and 6 subplots in a single figure window.
% Fixes: Planform Legend, 3D Transparency/Lighting, Robust Drag Calc.

clear; close all; clc;

% --- Setup Paths ---
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
emwetPath = fullfile(currentDir, 'EMWET');

% Validate Paths
if ~exist(q3dPath, 'dir')
    error('Q3D folder not found at: %s. Please run from the correct directory.', q3dPath);
end

addpath(currentDir); 
addpath(q3dPath);
addpath(emwetPath);

global W_wing
W_wing = 6344*9.81; % Initialize global weight

%% =========================================================================
%% 1. DEFINE DESIGN VECTORS
%% =========================================================================

% --- INITIAL DESIGN (x0) ---
b_init = 34.0;
c_r_init = 7.0;
c_k_init = 3.7;
c_t_init = 1.6;
M_cr_init = 0.78;
h_cr_init = 11278.4;
CST_init = [0.2337, 0.0796, 0.2683, 0.0887, 0.2789, 0.3811, -0.2254, -0.1634, -0.0470, -0.4771, 0.0735, 0.3255];
spar_locs = [0.2, 0.6];
tank_limits = [0, 0.85];
b_k_static = 4.36 + 3.95/2; 

% Calculate initial fuel guess
W_fuel_init = performance(b_init, c_r_init, c_k_init, c_t_init, b_k_static, spar_locs, tank_limits);

x_initial = [b_init, c_r_init, c_k_init, c_t_init, M_cr_init, h_cr_init, W_fuel_init, CST_init];

% --- USER OPTIMIZED DESIGN (x_opt) ---
b_optim = 33.48;
c_r_optim = 6.98;
c_k_optim = 4.21;
c_t_optim = 1.46;
M_cr_optim = 0.775;
h_cr_optim = 11458;
CST_optim = [0.1730, 0.1310, 0.2249, 0.0995, 0.3298, 0.2608, -0.2216, -0.1680, -0.1055, -0.3265, 0.1234, 0.3174];

% Calculate fuel for optimized design
W_fuel_optim = performance(b_optim, c_r_optim, c_k_optim, c_t_optim, b_k_static, spar_locs, tank_limits);

x_optimal = [b_optim, c_r_optim, c_k_optim, c_t_optim, M_cr_optim, h_cr_optim, W_fuel_optim, CST_optim];

%% =========================================================================
%% 2. RUN ANALYSIS
%% =========================================================================
fprintf('Analyzing Initial Design... (This may take a minute)\n');
Data_Init = analyze_design_full(x_initial, 'Initial');

fprintf('Analyzing Optimized Design... (This may take a minute)\n');
Data_Opt = analyze_design_full(x_optimal, 'Optimized');

%% =========================================================================
%% 3. GENERATE TABLES (Strict Assignment Order)
%% =========================================================================

% --- TABLE 1: OPTIMIZATION VARIABLES, CONSTRAINTS & WEIGHTS ---
fprintf('\n\n=======================================================================================\n');
fprintf('TABLE 1: OPTIMIZATION RESULTS (Variables, Constraints, Weights)\n');
fprintf('=======================================================================================\n');
fprintf('%-35s | %-15s | %-15s | %-8s\n', 'Metric', 'Initial', 'Optimized', 'Unit');
fprintf('---------------------------------------------------------------------------------------\n');
% Objective
fprintf('%-35s | %15.4f | %15.4f | km\n', 'Objective Function (Range)', Data_Init.Range/1000, Data_Opt.Range/1000);
fprintf('---------------------------------------------------------------------------------------\n');
% Design Variables
fprintf('%-35s | %15.4f | %15.4f | m\n', 'Wingspan (b)', Data_Init.b, Data_Opt.b);
fprintf('%-35s | %15.4f | %15.4f | m\n', 'Root Chord (c_r)', Data_Init.c_r, Data_Opt.c_r);
fprintf('%-35s | %15.4f | %15.4f | m\n', 'Kink Chord (c_k)', Data_Init.c_k, Data_Opt.c_k);
fprintf('%-35s | %15.4f | %15.4f | m\n', 'Tip Chord (c_t)', Data_Init.c_t, Data_Opt.c_t);
fprintf('%-35s | %15.4f | %15.4f | -', 'Cruise Mach (M_cr)', Data_Init.M_cr, Data_Opt.M_cr);
fprintf('%-35s | %15.0f | %15.0f | m\n', 'Cruise Altitude (h_cr)', Data_Init.h_cr, Data_Opt.h_cr);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'Fuel Weight (W_fuel)', Data_Init.W_fuel/9.81, Data_Opt.W_fuel/9.81);
fprintf('---------------------------------------------------------------------------------------\n');
% Constraints
fprintf('%-35s | %15.4f | %15.4f | -\n', 'Constr: Wing Loading (c<=0)', Data_Init.Constr_WL, Data_Opt.Constr_WL);
fprintf('%-35s | %15.4f | %15.4f | -\n', 'Constr: Tank Volume (c<=0)', Data_Init.Constr_Vol, Data_Opt.Constr_Vol);
fprintf('---------------------------------------------------------------------------------------\n');
% Weights & CO2
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'MTOW (W_TO_max)', Data_Init.W_TO_max/9.81, Data_Opt.W_TO_max/9.81);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'Wing Structure (W_str_wing)', Data_Init.W_wing/9.81, Data_Opt.W_wing/9.81);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'CO2 Emitted (W_CO2)', Data_Init.W_CO2, Data_Opt.W_CO2);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'Weight A-W-Fuel (W_AW)', Data_Init.W_AW/9.81, Data_Opt.W_AW/9.81);
fprintf('---------------------------------------------------------------------------------------\n');
% Volumes
fprintf('%-35s | %15.2f | %15.2f | m^3\n', 'Fuel Volume (V_fuel)', Data_Init.V_fuel_req, Data_Opt.V_fuel_req);
fprintf('%-35s | %15.2f | %15.2f | m^3\n', 'Tank Capacity (V_tank)', Data_Init.V_tank_cap, Data_Opt.V_tank_cap);
fprintf('---------------------------------------------------------------------------------------\n');

% --- TABLE 2: PHYSICS @ DESIGN POINT ---
fprintf('\n\n=======================================================================================\n');
fprintf('TABLE 2: PHYSICS & AERODYNAMICS (@ DESIGN POINT)\n');
fprintf('=======================================================================================\n');
fprintf('%-35s | %-15s | %-15s | %-8s\n', 'Metric', 'Initial', 'Optimized', 'Unit');
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.2f | %15.2f | m/s\n', 'Velocity (V)', Data_Init.V_cr, Data_Opt.V_cr);
fprintf('%-35s | %15.0f | %15.0f | m\n', 'Altitude (h)', Data_Init.h_cr, Data_Opt.h_cr);
fprintf('%-35s | %15.4f | %15.4f | -', 'Mach Number', Data_Init.M_cr, Data_Opt.M_cr);
fprintf('\n%-35s | %15.1f | %15.1f | Pa\n', 'Dynamic Pressure (q)', Data_Init.q_cr, Data_Opt.q_cr);
fprintf('%-35s | %15.2e | %15.2e | -', 'Reynolds Number (Re)', Data_Init.Re, Data_Opt.Re);
fprintf('\n---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.4e | %15.4e | 1/s\n', 'SFC (C_T)', Data_Init.Ct, Data_Opt.Ct);
fprintf('%-35s | %15.4f | %15.4f | -\n', 'Prop. Efficiency (eta)', Data_Init.eta, Data_Opt.eta);
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.4f | %15.4f | -\n', 'Lift Coeff (C_L)', Data_Init.CL, Data_Opt.CL);
fprintf('%-35s | %15.4f | %15.4f | deg\n', 'Angle of Attack (alpha)', Data_Init.alpha, Data_Opt.alpha);
fprintf('%-35s | %15.5f | %15.5f | -\n', 'Wing Drag Coeff (C_D,wing)', Data_Init.CD_wing, Data_Opt.CD_wing);
fprintf('%-35s | %15.5f | %15.5f | -\n', 'Induced Drag Coeff (C_Di,wing)', Data_Init.CD_ind, Data_Opt.CD_ind);
fprintf('%-35s | %15.2f | %15.2f | -\n', 'Efficiency (C_L/C_D)', Data_Init.L_D, Data_Opt.L_D);
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.1f | %15.1f | N\n', 'Drag A-W Force (D_AW)', Data_Init.D_AW, Data_Opt.D_AW);
fprintf('%-35s | %15.5f | %15.5f | -\n', 'Drag A-W Coeff (C_D,AW)', Data_Init.CD_AW, Data_Opt.CD_AW);
fprintf('---------------------------------------------------------------------------------------\n');

% --- TABLE 3: GEOMETRY ---
fprintf('\n\n=======================================================================================\n');
fprintf('TABLE 3: GEOMETRY DETAILS\n');
fprintf('=======================================================================================\n');
fprintf('%-35s | %-15s | %-15s | %-8s\n', 'Metric', 'Initial', 'Optimized', 'Unit');
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.2f | %15.2f | m^2\n', 'Wing Area (S)', Data_Init.S, Data_Opt.S);
fprintf('%-35s | %15.2f | %15.2f | m\n', 'MAC', Data_Init.MAC, Data_Opt.MAC);
fprintf('%-35s | %15.1f | %15.1f | kg/m^2\n', 'Wing Loading (W_TO_max/S)', Data_Init.WS, Data_Opt.WS);
fprintf('%-35s | %15.2f | %15.2f | -\n', 'Aspect Ratio (AR)', Data_Init.AR, Data_Opt.AR);
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.2f | %15.2f | deg\n', 'LE Sweep Angle', Data_Init.Sweep_LE, Data_Opt.Sweep_LE);
fprintf('%-35s | %15.2f | %15.2f | m\n', 'Inboard Span (b_in)', Data_Init.b_in, Data_Opt.b_in);
fprintf('%-35s | %15.2f | %15.2f | m\n', 'Outboard Span (b_out)', Data_Init.b_out, Data_Opt.b_out);
fprintf('---------------------------------------------------------------------------------------\n');

%% =========================================================================
%% 4. GENERATE PLOTS (Single Interface - 2x3 Grid)
%% =========================================================================
c_init = [0, 0.4470, 0.7410]; % Blue
c_opt  = [0.8500, 0.3250, 0.0980]; % Orange
lw = 1.5; ms = 4;

% Create one big figure
figure('Name', 'MDO Part 2 Results Overview', 'Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);

% --- 1. Spanwise Lift Distribution (Design Point) ---
subplot(2,3,1);
plot(Data_Init.Y_aero, Data_Init.cl_dist .* Data_Init.chord_dist, 'o-', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_init); hold on;
plot(Data_Opt.Y_aero, Data_Opt.cl_dist .* Data_Opt.chord_dist, 's-', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_opt);
title({'Spanwise Lift (c \cdot C_l)'; '@ Design Point'}, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Spanwise y [m]', 'FontSize', 10); ylabel('Lift Load [m]', 'FontSize', 10);
legend('Initial', 'Optimized', 'Location', 'southwest', 'FontSize', 8); grid on; set(gca, 'GridAlpha', 0.3);

% --- 2. Spanwise Lift Distribution (Critical) ---
subplot(2,3,2);
plot(Data_Init.Y_load, Data_Init.L_load, 'o-', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_init); hold on;
plot(Data_Opt.Y_load, Data_Opt.L_load, 's-', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_opt);
title({'Spanwise Lift Load'; '@ Critical Conditions'}, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Spanwise y [m]', 'FontSize', 10); ylabel('Lift Force [N/m]', 'FontSize', 10);
legend('Initial', 'Optimized', 'Location', 'southwest', 'FontSize', 8); grid on; set(gca, 'GridAlpha', 0.3);

% --- 3. Spanwise Drag Coefficients (Breakdown) ---
subplot(2,3,3);
% Calculate Section Drag LOAD (c * Cd)
Prof_Init_Dist = (Data_Init.cd_dist - Data_Init.cd_ind_dist) .* Data_Init.chord_dist;
Ind_Init_Dist  = Data_Init.cd_ind_dist .* Data_Init.chord_dist;
Prof_Opt_Dist = (Data_Opt.cd_dist - Data_Opt.cd_ind_dist) .* Data_Opt.chord_dist;
Ind_Opt_Dist  = Data_Opt.cd_ind_dist .* Data_Opt.chord_dist;

plot(Data_Init.Y_aero, Ind_Init_Dist, 'o--', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms); hold on;
plot(Data_Init.Y_aero, Prof_Init_Dist, 's-', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_init);
plot(Data_Opt.Y_aero, Ind_Opt_Dist, 'o--', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms);
plot(Data_Opt.Y_aero, Prof_Opt_Dist, 's-', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_opt);

title({'Spanwise Drag (c \cdot C_d)'; '@ Design Point'}, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Spanwise y [m]', 'FontSize', 10); ylabel('c \cdot C_d [m]', 'FontSize', 10);
legend('Init Induced', 'Init Prof+Wave', 'Opt Induced', 'Opt Prof+Wave', 'Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3);

% --- 4. Wing Airfoil Comparison (Simplified Legend) ---
subplot(2,3,4);
[xu_0, xl_0] = get_airfoil_coords(x_initial(8:end)); 
[xu_opt, xl_opt] = get_airfoil_coords(x_optimal(8:end));
% Dummy handles for simple legend
h1 = plot(nan, nan, '--', 'Color', c_init, 'LineWidth', lw); hold on;
h2 = plot(nan, nan, '-', 'Color', c_opt, 'LineWidth', lw+0.5);
% Actual Plots
plot(xu_0(1,:), xu_0(2,:), '--', 'Color', c_init, 'LineWidth', lw);
plot(xl_0(1,:), xl_0(2,:), '--', 'Color', c_init, 'LineWidth', lw);
plot(xu_opt(1,:), xu_opt(2,:), '-', 'Color', c_opt, 'LineWidth', lw+0.5);
plot(xl_opt(1,:), xl_opt(2,:), '-', 'Color', c_opt, 'LineWidth', lw+0.5);
axis equal; grid on; title('Normalized Airfoil', 'FontSize', 11, 'FontWeight', 'bold'); 
xlabel('x/c', 'FontSize', 10); ylabel('z/c', 'FontSize', 10); 
legend([h1, h2], 'Initial', 'Optimized', 'Location', 'best', 'FontSize', 8);

% --- 5. Wing Planform (Full, No Spars, Correct Legend) ---
subplot(2,3,5);
hold on; % CRITICAL: Must be ON before dummy plot
% Dummy handles for legend
h_p1 = plot(nan, nan, '--', 'Color', c_init, 'LineWidth', lw);
h_p2 = plot(nan, nan, '-', 'Color', c_opt, 'LineWidth', lw+0.5);

plot_planform_full(x_initial, '--', c_init, lw);
plot_planform_full(x_optimal, '-', c_opt, lw+0.5);

axis equal; grid on; title('Wing Planform', 'FontSize', 11, 'FontWeight', 'bold'); 
xlabel('Span y [m]', 'FontSize', 10); ylabel('Chord x [m]', 'FontSize', 10); 
legend([h_p1, h_p2], 'Initial', 'Optimized', 'Location', 'south', 'FontSize', 8);

% --- 6. Wing Isometric View (Fixed & Legend Fix) ---
subplot(2,3,6);
hold on;
% Dummy handles for legend
h1_3d = plot(nan, nan, 's', 'MarkerFaceColor', c_init, 'Color', c_init, 'MarkerSize', 10);
h2_3d = plot(nan, nan, 's', 'MarkerFaceColor', c_opt, 'Color', c_opt, 'MarkerSize', 10);
% Plots
plot_wing_3d(x_initial, c_init, 0.3); 
plot_wing_3d(x_optimal, c_opt, 0.45);
axis equal; grid on; title('3D Wing (Isometric)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 10); ylabel('y', 'FontSize', 10); zlabel('z', 'FontSize', 10);
view(135, 25); % Isometric View
legend([h1_3d, h2_3d], 'Initial', 'Optimized', 'Location', 'northeast', 'FontSize', 8);
b_max = max(Data_Init.b, Data_Opt.b);
xlim([-b_max/2-2, b_max/2+2]); ylim([0, 10]); zlim([-2, 4]);

%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================
function [D] = analyze_design_full(x, label)
    b = x(1); c_r = x(2); c_k = x(3); c_t = x(4);
    M_cr = x(5); h_cr = x(6); W_fuel = x(7); CST = x(8:19);
    
    b_k = 4.36 + 3.95/2;
    W_AminusW = 450000;
    spar_locs = [0.2, 0.6]; tank_limits = [0, 0.85];
    rho_fuel = 0.81715e3; 
    W_AW = W_AminusW; 
    
    % Geometry Calc
    S_half = (c_r + c_k) * b_k / 2 + (c_k + c_t) * (b/2 - b_k) / 2;
    S_wing = 2*S_half;
    AR = b^2 / S_wing;
    int_c2_dy = (b_k/3) * (c_r^2 + c_r*c_k + c_k^2) + ...
                ((b/2 - b_k)/3) * (c_k^2 + c_k*c_t + c_t^2);
    MAC = (1 / S_half) * int_c2_dy;
    
    % Sweep Calc
    x_le_r = 0;
    x_te_r = c_r;
    x_te_k = x_te_r + b_k * tan(deg2rad(0.01)); % Fixed TE sweep
    x_le_k = x_te_k - c_k;
    Sweep_LE = rad2deg(atan((x_le_k - x_le_r) / b_k));

    % Weight Convergence
    global W_wing
    W_wing = 6344*9.81; 
    V_MO_ref = 0.82 * sqrt(1.4*287*216.65); 
    n_max = 2.5;
    mat_props.upper = [7.1e10, 2795, 4.8e8, 4.6e8];
    mat_props.lower = mat_props.upper;
    mat_props.front = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
    mat_props.rear = mat_props.upper;
    engine_data.count = 1; engine_data.y_location = 4.7; engine_data.weight = 1969;
    
    cst2dat(CST, 'EMWET/optimized_airfoil');
    airfoils.root = 'EMWET/optimized_airfoil';
    airfoils.kink = 'EMWET/optimized_airfoil';
    airfoils.tip  = 'EMWET/optimized_airfoil';

    for iter = 1:5
        W_total = W_AminusW + 2*W_wing + W_fuel;
        MTOW = W_total / 9.81;
        ZFW = (W_AminusW + 2*W_wing) / 9.81;
        [Y_load, L_load, M_load] = loads(0.01, b_k, 5, 0, -1, -3, ...
                           n_max, V_MO_ref, W_AminusW, h_cr, ...
                           b, c_r, c_k, c_t, CST, W_wing, W_fuel);
        W_wing_new = 9.81 * structures(Y_load, L_load, M_load, ...
                            0.01, b_k, 5, ...
                            b, c_r, c_k, c_t, ...
                            MTOW, ZFW, n_max, S_wing/2, ...
                            mat_props, spar_locs, ...
                            tank_limits, engine_data, airfoils, ...
                            'temp_analysis', false);
        W_wing = W_wing_new;
    end
    
    % Aerodynamics
    [~, rho_cr, T_cr] = stdatm(h_cr);
    a_cr = sqrt(1.4 * 287 * T_cr);
    V_cr = M_cr * a_cr;
    Re = rho_cr * V_cr * MAC / (1.789e-5 * (288.15 + 110.4) / (T_cr + 110.4) * (T_cr / 288.15)^(1.5));
    q_cr = 0.5 * rho_cr * V_cr^2;
    
    Res_Aero = run_Q3D_internal(b, c_r, c_k, c_t, b_k, CST, h_cr, M_cr, W_total, S_wing);
    
    % Performance
    CD0_rest = 0.012; 
    D_AW = CD0_rest * q_cr * S_wing;
    CD_total = Res_Aero.CDwing + CD0_rest;
    L_D = Res_Aero.CLwing / CD_total;
    
    V_cr_ref = 230.1563; h_cr_ref = 11278; C_T_ref = 1.8639e-4;
    eta = 1 * exp(-((V_cr - V_cr_ref)^2)/(2*70^2) - ((h_cr - h_cr_ref)^2)/(2*2500^2));
    Ct = C_T_ref / eta;
    W_start = W_total;
    W_end = (1 - W_fuel/W_total) * W_start / 0.938;
    Range = (V_cr / Ct) * L_D * log(W_start / W_end);
    
    W_fuel_capacity = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);
    V_tank_cap = (W_fuel_capacity / 9.81) / (rho_fuel * 0.93); 
    V_fuel_req = (W_fuel / 9.81) / rho_fuel;
    
    % Constraints (Approx for output)
    b_orig = 34.0; c_r_orig = 7.0; c_k_orig = 3.7; c_t_orig = 1.6;
    S_orig = (c_r_orig+c_k_orig)*b_k + (c_k_orig+c_t_orig)*(b_orig/2-b_k);
    WL_current = W_total/S_wing; WL_orig = W_total/S_orig;
    Constr_WL = (WL_current - WL_orig)/WL_orig;
    W_fuel_max_orig = 16085 * 9.81;
    Constr_Vol = (W_fuel - W_fuel_capacity)/W_fuel_capacity;

    % Pack Data
    D.Range = Range; D.b = b; D.c_r = c_r; D.c_k = c_k; D.c_t = c_t;
    D.M_cr = M_cr; D.h_cr = h_cr; D.W_fuel = W_fuel; D.W_TO_max = W_total;
    D.W_wing = W_wing; D.W_CO2 = (W_fuel/9.81) * 3.16;
    D.V_fuel_req = V_fuel_req; D.V_tank_cap = V_tank_cap;
    D.V_cr = V_cr; D.q_cr = q_cr; D.Re = Re;
    D.CL = Res_Aero.CLwing; D.alpha = Res_Aero.Alpha;
    D.CD_wing = Res_Aero.CDwing;
    D.CD_AW = CD0_rest; D.D_AW = D_AW; D.W_AW = W_AW;
    D.L_D = L_D; D.Ct = Ct; D.eta = eta;
    D.S = S_wing; D.AR = AR; D.WS = WL_current; D.MAC = MAC;
    D.Sweep_LE = Sweep_LE; D.b_in = b_k; D.b_out = b/2 - b_k;
    D.Y_load = Y_load; D.L_load = L_load; D.n_max = n_max;
    D.Constr_WL = Constr_WL; D.Constr_Vol = Constr_Vol;

    % Extract Arrays (No Approximation)
    if isfield(Res_Aero, 'Wing')
        D.Y_aero = Res_Aero.Wing.Yst;
        D.cl_dist = Res_Aero.Wing.cl;
        D.chord_dist = Res_Aero.Wing.chord;
        
        if isfield(Res_Aero.Wing, 'cd'), D.cd_dist = Res_Aero.Wing.cd; else, D.cd_dist = zeros(size(D.Y_aero)); end
        if isfield(Res_Aero.Wing, 'cd_i'), D.cd_ind_dist = Res_Aero.Wing.cd_i; elseif isfield(Res_Aero.Wing, 'cdi'), D.cd_ind_dist = Res_Aero.Wing.cdi; else, D.cd_ind_dist = zeros(size(D.Y_aero)); end
        
        if isfield(Res_Aero, 'CDind'), D.CD_ind = Res_Aero.CDind; else, D.CD_ind = trapz(D.Y_aero, D.chord_dist .* D.cd_ind_dist) * 2 / S_wing; end
    else
        % Fallback for older Q3D versions
        D.Y_aero = Res_Aero.Yst; D.cl_dist = Res_Aero.cl; D.chord_dist = Res_Aero.chord;
        D.cd_dist = Res_Aero.cd; D.cd_ind_dist = Res_Aero.cd_i;
        D.CD_ind = Res_Aero.CDind;
    end
end

function [Res] = run_Q3D_internal(b, c_r, c_k, c_t, b_k, CST, h_cr, M_cr, W_total, S_wing)
    x_le_r = 0; x_te_r = c_r; 
    x_te_k = x_te_r + b_k * tan(deg2rad(0.01)); x_le_k = x_te_k - c_k;
    tan_sweep_le = (x_le_k - x_le_r) / b_k;
    x_le_t = x_le_r + (b/2) * tan_sweep_le;
    z_k = b_k * tan(deg2rad(5)); z_t = (b/2) * tan(deg2rad(5));

    AC.Wing.Geom = [x_le_r, 0, 0, c_r, 0;
                    x_le_k, b_k, z_k, c_k, -1.0;
                    x_le_t, b/2, z_t, c_t, -3];
    AC.Wing.inc = 0;   
    AC.Wing.Airfoils = [CST; CST; CST];
    AC.Wing.eta = [0; b_k/(b/2); 1];
    AC.Visc = 1; % FORCE VISCOUS
    AC.Aero.MaxIterIndex = 500; % INCREASED to wait for convergence

    [~, rho_cr, T_cr] = stdatm(h_cr);
    a_cr = sqrt(1.4 * 287 * T_cr);
    AC.Aero.V = M_cr * a_cr;
    AC.Aero.rho = rho_cr;
    AC.Aero.alt = h_cr;
    AC.Aero.M = M_cr;
    AC.Aero.CL = W_total / (0.5 * rho_cr * AC.Aero.V^2 * S_wing);

    q3dFolder = fullfile(pwd, 'Q3D');
    mockFile = fullfile(q3dFolder, 'atmosisa.m');
    fid = fopen(mockFile, 'w');
    if fid ~= -1
        fprintf(fid, 'function [T, a, P, rho] = atmosisa(h)\n    [P, rho, T] = stdatm(h);\n    gamma = 1.4; R = 287.058;\n    a = sqrt(gamma * R * T);\nend\n');
        fclose(fid);
    end
    
    originalD = pwd;
    cd(q3dFolder); 
    if ispc, [~, ~] = system('taskkill /F /IM xfoil.exe 2>NUL'); pause(0.2); end
    try
        Res = Q3D_solver(AC);
    catch ME
        cd(originalD); if exist(mockFile, 'file'), delete(mockFile); end; rethrow(ME);
    end
    cd(originalD); if exist(mockFile, 'file'), delete(mockFile); end
end

function [xu, xl] = get_airfoil_coords(CST)
    x = linspace(0,1,100);
    N1 = 0.5; N2 = 1.0;
    C_class = (x.^N1) .* ((1-x).^N2);
    Au = CST(1:6); Al = CST(7:12);
    S_u = 0; S_l = 0;
    for i = 1:6
       n_fac = factorial(6-1)/(factorial(i-1)*factorial(6-i));
       S_u = S_u + Au(i) * n_fac * (x.^(i-1)) .* ((1-x).^(6-i));
       S_l = S_l + Al(i) * n_fac * (x.^(i-1)) .* ((1-x).^(6-i));
    end
    yu = C_class .* S_u; yl = C_class .* S_l;
    xu = [x; yu]; xl = [x; yl];
end

function plot_planform_full(x, style, color, lw)
    b = x(1); c_r = x(2); c_k = x(3); c_t = x(4);
    b_k = 4.36 + 3.95/2;
    sweep_te = 0.01;
    x_r = 0; 
    x_k = c_r + b_k*tan(deg2rad(sweep_te)) - c_k;
    sweep_le = atan(x_k/b_k);
    x_t = (b/2)*tan(sweep_le);
    
    % Right Wing Coordinates
    X = [0, x_k, x_t, x_t+c_t, x_k+c_k, c_r, 0];
    Y = [0, b_k, b/2, b/2, b_k, 0, 0];
    
    % Left Wing Coordinates (Mirror Y)
    Y_left = -Y;
    
    % Plot Right
    plot(Y, X, style, 'Color', color, 'LineWidth', lw);
    % Plot Left
    plot(Y_left, X, style, 'Color', color, 'LineWidth', lw);
end

function plot_wing_3d(x, color, alpha)
    b = x(1); c_r = x(2); c_k = x(3); c_t = x(4);
    CST = x(8:19);
    b_k = 4.36 + 3.95/2; sweep_te = 0.01; dihedral = 5;
    x_te_r = c_r; x_te_k = x_te_r + b_k * tan(deg2rad(sweep_te));
    x_le_k = x_te_k - c_k; tan_sweep_le = x_le_k / b_k; x_le_t = (b/2) * tan_sweep_le;
    z_k = b_k * tan(deg2rad(dihedral)); z_t = (b/2) * tan(deg2rad(dihedral));
    [xu, xl] = get_airfoil_coords(CST);
    x_af = [xu(1,:), fliplr(xl(1,:))]; z_af = [xu(2,:), fliplr(xl(2,:))];
    Stations = [0, 0, 0, c_r, 0; b_k, x_le_k, z_k, c_k, -1.0; b/2, x_le_t, z_t, c_t, -3.0];
    X_surf = []; Y_surf = []; Z_surf = [];
    for i = 1:size(Stations, 1)
        y_loc = Stations(i,1); x_le = Stations(i,2); z_le = Stations(i,3); chord = Stations(i,4); twist = Stations(i,5);
        theta = deg2rad(-twist); R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        coords = R * [x_af; z_af];
        X_sect = x_le + coords(1,:) * chord; Z_sect = z_le + coords(2,:) * chord; Y_sect = ones(size(X_sect)) * y_loc;
        X_surf = [X_surf; X_sect]; Y_surf = [Y_surf; Y_sect]; Z_surf = [Z_surf; Z_sect];
    end
    
    % Draw Right Wing
    surf(Y_surf, X_surf, Z_surf, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    % Draw Left Wing (Mirrored)
    surf(-Y_surf, X_surf, Z_surf, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    
    % Brighter lighting setup
    camlight('headlight'); 
    material dull; 
    lighting gouraud;
end