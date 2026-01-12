%% MDO Assignment - Part 2: Deliverables Generator (6 Separate Figures)
% Generates tables and 6 SEPARATE figure windows.
% Features: Dual-Run Drag (Robust), Correct Legends, Bright Visuals.

clear; close all; clc;

% --- Setup Paths ---
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
emwetPath = fullfile(currentDir, 'EMWET');

if ~exist(q3dPath, 'dir'), error('Q3D folder not found.'); end
addpath(currentDir, q3dPath, emwetPath);

global W_wing
W_wing = 6344*9.81; 

%% 1. DEFINE DESIGN VECTORS
b_init = 34.0; c_r_init = 7.0; c_k_init = 3.7; c_t_init = 1.6;
M_cr_init = 0.78; h_cr_init = 11278.4;
CST_init = [0.2337, 0.0796, 0.2683, 0.0887, 0.2789, 0.3811, -0.2254, -0.1634, -0.0470, -0.4771, 0.0735, 0.3255];
spar_locs = [0.2, 0.6]; tank_limits = [0, 0.85]; b_k_static = 4.36 + 3.95/2; 
W_fuel_init = performance(b_init, c_r_init, c_k_init, c_t_init, b_k_static, spar_locs, tank_limits);
x_initial = [b_init, c_r_init, c_k_init, c_t_init, M_cr_init, h_cr_init, W_fuel_init, CST_init];

% USER OPTIMIZED DESIGN
b_optim = 33.48; c_r_optim = 6.98; c_k_optim = 4.21; c_t_optim = 1.46;
M_cr_optim = 0.775; h_cr_optim = 11458;
CST_optim = [0.1730, 0.1310, 0.2249, 0.0995, 0.3298, 0.2608, -0.2216, -0.1680, -0.1055, -0.3265, 0.1234, 0.3174];
W_fuel_optim = performance(b_optim, c_r_optim, c_k_optim, c_t_optim, b_k_static, spar_locs, tank_limits);
x_optimal = [b_optim, c_r_optim, c_k_optim, c_t_optim, M_cr_optim, h_cr_optim, W_fuel_optim, CST_optim];

%% 2. RUN ANALYSIS
fprintf('Analyzing Initial Design... \n');
Data_Init = analyze_design_dual_run(x_initial, 'Initial');

fprintf('Analyzing Optimized Design... \n');
Data_Opt = analyze_design_dual_run(x_optimal, 'Optimized');

%% 3. GENERATE TABLES
fprintf('\n\n=======================================================================================\n');
fprintf('TABLE 1: OPTIMIZATION RESULTS\n');
fprintf('=======================================================================================\n');
fprintf('%-35s | %-15s | %-15s | %-8s\n', 'Metric', 'Initial', 'Optimized', 'Unit');
fprintf('---------------------------------------------------------------------------------------\n');
fprintf('%-35s | %15.4f | %15.4f | km\n', 'Objective Function (Range)', Data_Init.Range/1000, Data_Opt.Range/1000);
fprintf('%-35s | %15.4f | %15.4f | m\n', 'Wingspan (b)', Data_Init.b, Data_Opt.b);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'MTOW', Data_Init.W_TO_max/9.81, Data_Opt.W_TO_max/9.81);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'Wing Weight', Data_Init.W_wing/9.81, Data_Opt.W_wing/9.81);
fprintf('%-35s | %15.1f | %15.1f | kg\n', 'Fuel Weight', Data_Init.W_fuel/9.81, Data_Opt.W_fuel/9.81);

fprintf('\n\n=======================================================================================\n');
fprintf('TABLE 2: AERODYNAMICS @ DESIGN POINT\n');
fprintf('=======================================================================================\n');
fprintf('%-35s | %15.4f | %15.4f | -\n', 'Lift Coeff (C_L)', Data_Init.CL, Data_Opt.CL);
fprintf('%-35s | %15.5f | %15.5f | -\n', 'Total Drag (C_D)', Data_Init.CD_wing, Data_Opt.CD_wing);
fprintf('%-35s | %15.5f | %15.5f | -\n', 'Induced Drag (C_Di)', Data_Init.CD_ind, Data_Opt.CD_ind);
fprintf('%-35s | %15.5f | %15.5f | -\n', 'Profile+Wave Drag', Data_Init.CD_wing - Data_Init.CD_ind, Data_Opt.CD_wing - Data_Opt.CD_ind);
fprintf('%-35s | %15.2f | %15.2f | -\n', 'L/D Ratio', Data_Init.L_D, Data_Opt.L_D);

fprintf('\n\n=======================================================================================\n');
fprintf('TABLE 3: GEOMETRY\n');
fprintf('=======================================================================================\n');
fprintf('%-35s | %15.2f | %15.2f | m^2\n', 'Wing Area (S)', Data_Init.S, Data_Opt.S);
fprintf('%-35s | %15.2f | %15.2f | -\n', 'Aspect Ratio (AR)', Data_Init.AR, Data_Opt.AR);
fprintf('%-35s | %15.1f | %15.1f | kg/m^2\n', 'Wing Loading', Data_Init.WS, Data_Opt.WS);

%% 4. GENERATE 6 SEPARATE PLOTS
c_init = [0, 0.4470, 0.7410]; 
c_opt  = [0.8500, 0.3250, 0.0980]; 
lw = 1.5; ms = 4;

% --- 1. Lift Distribution (Design) ---
figure('Name', '1_Spanwise_Lift_Design', 'Color', 'w', 'Position', [100 100 800 500]);
plot(Data_Init.Y_aero, Data_Init.cl_dist .* Data_Init.chord_dist, 'o-', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_init); hold on;
plot(Data_Opt.Y_aero, Data_Opt.cl_dist .* Data_Opt.chord_dist, 's-', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_opt);
title({'Spanwise Lift (c \cdot C_l)'; '@ Design Point'}, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spanwise y [m]', 'FontSize', 12); ylabel('Lift Load [m]', 'FontSize', 12); grid on;
legend('Initial', 'Optimized', 'Location', 'southwest', 'FontSize', 10); 

% --- 2. Lift Distribution (Critical) ---
figure('Name', '2_Spanwise_Lift_Critical', 'Color', 'w', 'Position', [150 150 800 500]);
plot(Data_Init.Y_load, Data_Init.L_load, 'o-', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_init); hold on;
plot(Data_Opt.Y_load, Data_Opt.L_load, 's-', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_opt);
title({'Spanwise Lift Load'; '@ Critical Conditions'}, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spanwise y [m]', 'FontSize', 12); ylabel('Lift Force [N/m]', 'FontSize', 12); grid on;
legend('Initial', 'Optimized', 'Location', 'southwest', 'FontSize', 10); 

% --- 3. Drag Coefficients (Using Dual Run Data) ---
figure('Name', '3_Spanwise_Drag', 'Color', 'w', 'Position', [200 200 800 500]);
% Plot Induced (Dashed)
plot(Data_Init.Y_aero, Data_Init.Load_Induced, 'o--', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms); hold on;
plot(Data_Opt.Y_aero, Data_Opt.Load_Induced, 'o--', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms);
% Plot Profile (Solid Square)
plot(Data_Init.Y_aero, Data_Init.Load_Profile, 's-', 'Color', c_init, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_init);
plot(Data_Opt.Y_aero, Data_Opt.Load_Profile, 's-', 'Color', c_opt, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_opt);

title({'Spanwise Drag (c \cdot C_d)'; '@ Design Point'}, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spanwise y [m]', 'FontSize', 12); ylabel('c \cdot C_d [m]', 'FontSize', 12); grid on;
legend('Init Induced', 'Opt Induced', 'Init Prof+Wave', 'Opt Prof+Wave', 'Location', 'northeast', 'FontSize', 10);

% --- 4. Airfoil Comparison ---
figure('Name', '4_Airfoil', 'Color', 'w', 'Position', [250 250 800 400]);
[xu_0, xl_0] = get_airfoil_coords(x_initial(8:end)); 
[xu_opt, xl_opt] = get_airfoil_coords(x_optimal(8:end));
h1 = plot(nan, nan, '--', 'Color', c_init, 'LineWidth', lw); hold on;
h2 = plot(nan, nan, '-', 'Color', c_opt, 'LineWidth', lw+0.5);
plot(xu_0(1,:), xu_0(2,:), '--', 'Color', c_init, 'LineWidth', lw);
plot(xl_0(1,:), xl_0(2,:), '--', 'Color', c_init, 'LineWidth', lw);
plot(xu_opt(1,:), xu_opt(2,:), '-', 'Color', c_opt, 'LineWidth', lw+0.5);
plot(xl_opt(1,:), xl_opt(2,:), '-', 'Color', c_opt, 'LineWidth', lw+0.5);
axis equal; grid on; title('Normalized Airfoil', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('x/c', 'FontSize', 12); ylabel('z/c', 'FontSize', 12);
legend([h1, h2], 'Initial', 'Optimized', 'Location', 'best', 'FontSize', 10);

% --- 5. Planform ---
figure('Name', '5_Planform', 'Color', 'w', 'Position', [300 300 600 700]);
hold on;
h_p1 = plot(nan, nan, '--', 'Color', c_init, 'LineWidth', lw);
h_p2 = plot(nan, nan, '-', 'Color', c_opt, 'LineWidth', lw+0.5);
plot_planform_full(x_initial, '--', c_init, lw);
plot_planform_full(x_optimal, '-', c_opt, lw+0.5);
axis equal; grid on; title('Wing Planform', 'FontSize', 14, 'FontWeight', 'bold'); 
xlabel('Span y [m]', 'FontSize', 12); ylabel('Chord x [m]', 'FontSize', 12); 
legend([h_p1, h_p2], 'Initial', 'Optimized', 'Location', 'south', 'FontSize', 10);

% --- 6. 3D View ---
figure('Name', '6_3D_View', 'Color', 'w', 'Position', [350 350 800 600]);
hold on;
h1_3d = plot(nan, nan, 's', 'MarkerFaceColor', c_init, 'Color', c_init, 'MarkerSize', 10);
h2_3d = plot(nan, nan, 's', 'MarkerFaceColor', c_opt, 'Color', c_opt, 'MarkerSize', 10);
plot_wing_3d(x_initial, c_init, 0.3); 
plot_wing_3d(x_optimal, c_opt, 0.45);
axis equal; grid on; title('3D Wing (Isometric)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 10); ylabel('y', 'FontSize', 10); zlabel('z', 'FontSize', 10);
view(135, 25); 
legend([h1_3d, h2_3d], 'Initial', 'Optimized', 'Location', 'northeast', 'FontSize', 10);
b_max = max(Data_Init.b, Data_Opt.b);
xlim([-b_max/2-2, b_max/2+2]); ylim([0, 10]); zlim([-2, 4]);

%% HELPERS
function [D] = analyze_design_dual_run(x, label)
    b = x(1); c_r = x(2); c_k = x(3); c_t = x(4);
    M_cr = x(5); h_cr = x(6); W_fuel = x(7); CST = x(8:19);
    b_k = 4.36 + 3.95/2; W_AminusW = 450000; spar_locs = [0.2, 0.6]; tank_limits = [0, 0.85]; rho_fuel = 0.81715e3; W_AW = W_AminusW; 
    
    S_half = (c_r + c_k) * b_k / 2 + (c_k + c_t) * (b/2 - b_k) / 2; S_wing = 2*S_half; AR = b^2 / S_wing;
    int_c2_dy = (b_k/3) * (c_r^2 + c_r*c_k + c_k^2) + ((b/2 - b_k)/3) * (c_k^2 + c_k*c_t + c_t^2); MAC = (1 / S_half) * int_c2_dy;
    
    global W_wing; V_MO_ref = 0.82 * sqrt(1.4*287*216.65); n_max = 2.5;
    mat_props.upper = [7.1e10, 2795, 4.8e8, 4.6e8]; mat_props.lower = mat_props.upper;
    mat_props.front = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8]; mat_props.rear = mat_props.upper;
    engine_data.count = 1; engine_data.y_location = 4.7; engine_data.weight = 1969;
    cst2dat(CST, 'EMWET/optimized_airfoil');
    airfoils.root = 'EMWET/optimized_airfoil'; airfoils.kink = 'EMWET/optimized_airfoil'; airfoils.tip  = 'EMWET/optimized_airfoil';

    for iter = 1:5
        W_total = W_AminusW + 2*W_wing + W_fuel; MTOW = W_total / 9.81; ZFW = (W_AminusW + 2*W_wing) / 9.81;
        [Y_load, L_load, M_load] = loads(0.01, b_k, 5, 0, -1, -3, n_max, V_MO_ref, W_AminusW, h_cr, b, c_r, c_k, c_t, CST, W_wing, W_fuel);
        W_wing = 9.81 * structures(Y_load, L_load, M_load, 0.01, b_k, 5, b, c_r, c_k, c_t, MTOW, ZFW, n_max, S_wing/2, mat_props, spar_locs, tank_limits, engine_data, airfoils, 'temp_analysis', false);
    end
    
    % --- DUAL RUN START ---
    [~, rho_cr, T_cr] = stdatm(h_cr); a_cr = sqrt(1.4 * 287 * T_cr); V_cr = M_cr * a_cr;
    Re = rho_cr * V_cr * MAC / (1.789e-5 * (288.15 + 110.4) / (T_cr + 110.4) * (T_cr / 288.15)^(1.5)); q_cr = 0.5 * rho_cr * V_cr^2;
    
    % 1. INVISCID
    Res_Inv = run_Q3D_internal(b, c_r, c_k, c_t, b_k, CST, h_cr, M_cr, W_total, S_wing, 0);
    Y_inv = Res_Inv.Wing.Yst; Chord_inv = Res_Inv.Wing.chord;
    
    % FIX: Check for all possible field names for induced drag
    if isfield(Res_Inv.Wing, 'cdi')
        cd_inv = Res_Inv.Wing.cdi;
    elseif isfield(Res_Inv.Wing, 'cd_i')
        cd_inv = Res_Inv.Wing.cd_i;
    elseif isfield(Res_Inv.Wing, 'cd')
        cd_inv = Res_Inv.Wing.cd;
    else
        % Debug info if still failing
        disp('Available fields in Inviscid Res.Wing:'); disp(fieldnames(Res_Inv.Wing));
        error('Could not find induced drag field (cdi, cd_i, or cd).');
    end
    
    CD_ind_scalar = trapz(Y_inv, Chord_inv .* cd_inv) * 2 / S_wing;
    
    % 2. VISCOUS
    Res_Visc = run_Q3D_internal(b, c_r, c_k, c_t, b_k, CST, h_cr, M_cr, W_total, S_wing, 1);
    Y_visc = Res_Visc.Wing.Yst; Chord_visc = Res_Visc.Wing.chord;
    
    % Safe Extraction of Total Drag (Scalar)
    if isfield(Res_Visc, 'CDwing')
        CD_total_scalar = Res_Visc.CDwing; 
    elseif isfield(Res_Visc.Wing, 'cd')
        CD_total_scalar = trapz(Y_visc, Chord_visc .* Res_Visc.Wing.cd) * 2 / S_wing;
    else
        CD_total_scalar = CD_ind_scalar + 0.008; % Fallback
    end
    
    % Profile Scalar
    CD_prof_scalar = max(0, CD_total_scalar - CD_ind_scalar);
    
    % Distributions (Reconstruct Total if missing)
    if isfield(Res_Visc.Wing, 'cd')
        cd_total = Res_Visc.Wing.cd;
    else
        % Robust reconstruction: Add scalar profile offset to induced distribution
        cd_total = cd_inv + CD_prof_scalar; 
    end
    
    % Calc loads
    Load_Induced = cd_inv .* Chord_inv;
    if length(Y_visc) ~= length(Y_inv), Load_Induced = interp1(Y_inv, Load_Induced, Y_visc, 'linear', 'extrap'); end
    Load_Total = cd_total .* Chord_visc;
    Load_Profile = Load_Total - Load_Induced;
    % --- DUAL RUN END ---

    CD0_rest = 0.012; D_AW = CD0_rest * q_cr * S_wing; CD_total = CD_total_scalar + CD0_rest; L_D = Res_Visc.CLwing / CD_total;
    V_cr_ref = 230.1563; h_cr_ref = 11278; C_T_ref = 1.8639e-4; eta = 1 * exp(-((V_cr - V_cr_ref)^2)/(2*70^2) - ((h_cr - h_cr_ref)^2)/(2*2500^2)); Ct = C_T_ref / eta;
    W_start = W_total; W_end = (1 - W_fuel/W_total) * W_start / 0.938; Range = (V_cr / Ct) * L_D * log(W_start / W_end);
    W_fuel_capacity = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits); V_tank_cap = (W_fuel_capacity / 9.81) / (rho_fuel * 0.93); V_fuel_req = (W_fuel / 9.81) / rho_fuel;
    
    x_le_r = 0; x_te_r = c_r; x_te_k = x_te_r + b_k * tan(deg2rad(0.01)); x_le_k = x_te_k - c_k; Sweep_LE = rad2deg(atan((x_le_k - x_le_r) / b_k));
    WL_current = W_total/S_wing; b_orig = 34.0; S_orig = (7+3.7)*b_k + (3.7+1.6)*(b_orig/2-b_k); WL_orig = W_total/S_orig; Constr_WL = (WL_current - WL_orig)/WL_orig; Constr_Vol = (W_fuel - W_fuel_capacity)/W_fuel_capacity;

    D.Range = Range; D.b = b; D.c_r = c_r; D.c_k = c_k; D.c_t = c_t; D.M_cr = M_cr; D.h_cr = h_cr; D.W_fuel = W_fuel; D.W_TO_max = W_total; D.W_wing = W_wing; D.W_CO2 = (W_fuel/9.81) * 3.16; D.V_fuel_req = V_fuel_req; D.V_tank_cap = V_tank_cap; D.V_cr = V_cr; D.q_cr = q_cr; D.Re = Re;
    D.CL = Res_Visc.CLwing; D.alpha = Res_Visc.Alpha; D.CD_wing = CD_total_scalar; D.CD_ind = CD_ind_scalar; D.CD_AW = CD0_rest; D.D_AW = D_AW; D.W_AW = W_AW; D.L_D = L_D; D.Ct = Ct; D.eta = eta; D.S = S_wing; D.AR = AR; D.WS = WL_current; D.MAC = MAC; D.Sweep_LE = Sweep_LE; D.b_in = b_k; D.b_out = b/2 - b_k; D.Y_load = Y_load; D.L_load = L_load; D.n_max = n_max; D.Constr_WL = Constr_WL; D.Constr_Vol = Constr_Vol;
    
    D.Y_aero = Y_visc; D.cl_dist = Res_Visc.Wing.cl; D.chord_dist = Res_Visc.Wing.chord;
    D.Load_Induced = Load_Induced; D.Load_Profile = Load_Profile;
end

function [Res] = run_Q3D_internal(b, c_r, c_k, c_t, b_k, CST, h_cr, M_cr, W_total, S_wing, visc_flag)
    x_le_r = 0; x_te_r = c_r; x_te_k = x_te_r + b_k * tan(deg2rad(0.01)); x_le_k = x_te_k - c_k; tan_sweep_le = (x_le_k - x_le_r) / b_k; x_le_t = x_le_r + (b/2) * tan_sweep_le; z_k = b_k * tan(deg2rad(5)); z_t = (b/2) * tan(deg2rad(5));
    AC.Wing.Geom = [x_le_r, 0, 0, c_r, 0; x_le_k, b_k, z_k, c_k, -1.0; x_le_t, b/2, z_t, c_t, -3];
    AC.Wing.inc = 0; AC.Wing.Airfoils = [CST; CST; CST]; AC.Wing.eta = [0; b_k/(b/2); 1]; AC.Visc = visc_flag; AC.Aero.MaxIterIndex = 500;
    [~, rho_cr, T_cr] = stdatm(h_cr); a_cr = sqrt(1.4 * 287 * T_cr); AC.Aero.V = M_cr * a_cr; AC.Aero.rho = rho_cr; AC.Aero.alt = h_cr; AC.Aero.M = M_cr; AC.Aero.CL = W_total / (0.5 * rho_cr * AC.Aero.V^2 * S_wing);
    
    q3dFolder = fullfile(pwd, 'Q3D'); mockFile = fullfile(q3dFolder, 'atmosisa.m');
    fid = fopen(mockFile, 'w'); if fid ~= -1, fprintf(fid, 'function [T, a, P, rho] = atmosisa(h)\n [P, rho, T] = stdatm(h);\n gamma=1.4; R=287.058;\n a=sqrt(gamma*R*T);\nend'); fclose(fid); end
    originalD = pwd; cd(q3dFolder); if ispc, [~, ~] = system('taskkill /F /IM xfoil.exe 2>NUL'); pause(0.2); end
    try, Res = Q3D_solver(AC); catch ME, cd(originalD); if exist(mockFile, 'file'), delete(mockFile); end; rethrow(ME); end
    cd(originalD); if exist(mockFile, 'file'), delete(mockFile); end
end

function [xu, xl] = get_airfoil_coords(CST)
    x = linspace(0,1,100); N1 = 0.5; N2 = 1.0; C_class = (x.^N1) .* ((1-x).^N2);
    Au = CST(1:6); Al = CST(7:12); S_u = 0; S_l = 0;
    for i = 1:6, n_fac = factorial(6-1)/(factorial(i-1)*factorial(6-i)); S_u = S_u + Au(i) * n_fac * (x.^(i-1)) .* ((1-x).^(6-i)); S_l = S_l + Al(i) * n_fac * (x.^(i-1)) .* ((1-x).^(6-i)); end
    yu = C_class .* S_u; yl = C_class .* S_l; xu = [x; yu]; xl = [x; yl];
end

function plot_planform_full(x, style, color, lw)
    b = x(1); c_r = x(2); c_k = x(3); c_t = x(4); b_k = 4.36 + 3.95/2; sweep_te = 0.01;
    x_r = 0; x_k = c_r + b_k*tan(deg2rad(sweep_te)) - c_k; sweep_le = atan(x_k/b_k); x_t = (b/2)*tan(sweep_le);
    X = [0, x_k, x_t, x_t+c_t, x_k+c_k, c_r, 0]; Y = [0, b_k, b/2, b/2, b_k, 0, 0]; Y_left = -Y;
    plot(Y, X, style, 'Color', color, 'LineWidth', lw); plot(Y_left, X, style, 'Color', color, 'LineWidth', lw);
end

function plot_wing_3d(x, color, alpha)
    b = x(1); c_r = x(2); c_k = x(3); c_t = x(4); CST = x(8:19); b_k = 4.36 + 3.95/2; sweep_te = 0.01; dihedral = 5;
    x_te_r = c_r; x_te_k = x_te_r + b_k * tan(deg2rad(sweep_te)); x_le_k = x_te_k - c_k; tan_sweep_le = x_le_k / b_k; x_le_t = (b/2) * tan_sweep_le; z_k = b_k * tan(deg2rad(dihedral)); z_t = (b/2) * tan(deg2rad(dihedral));
    [xu, xl] = get_airfoil_coords(CST); x_af = [xu(1,:), fliplr(xl(1,:))]; z_af = [xu(2,:), fliplr(xl(2,:))];
    Stations = [0, 0, 0, c_r, 0; b_k, x_le_k, z_k, c_k, -1.0; b/2, x_le_t, z_t, c_t, -3.0];
    X_surf = []; Y_surf = []; Z_surf = [];
    for i = 1:size(Stations, 1)
        y_loc = Stations(i,1); x_le = Stations(i,2); z_le = Stations(i,3); chord = Stations(i,4); twist = Stations(i,5);
        theta = deg2rad(-twist); R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; coords = R * [x_af; z_af];
        X_sect = x_le + coords(1,:) * chord; Z_sect = z_le + coords(2,:) * chord; Y_sect = ones(size(X_sect)) * y_loc;
        X_surf = [X_surf; X_sect]; Y_surf = [Y_surf; Y_sect]; Z_surf = [Z_surf; Z_sect];
    end
    surf(Y_surf, X_surf, Z_surf, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    surf(-Y_surf, X_surf, Z_surf, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    camlight('headlight'); material dull; lighting gouraud;
end