function [Range] = optimize(x)
% OPTIMIZE Wrapper function to run MDO analysis for given design variables
% Includes caching to avoid redundant calculations during optimization

% --- CACHING MECHANISM ---
persistent cache_x cache_Range cache_count;
if isempty(cache_count)
    cache_count = 0;
    cache_x = [];
    cache_Range = [];
end

% Check if this input has been evaluated before
cache_tol = 1e-5;  % Tolerance for matching cached inputs
for i = 1:size(cache_x, 1)
    if max(abs(cache_x(i,:) - x(:)')) < cache_tol
        Range = cache_Range(i);
        fprintf('[CACHE HIT] Returning cached result for evaluation #%d\n', i);
        return;
    end
end

% Extract design variables from input vector x
% x = [b, c_r, c_k, c_t, M_cr, h_cr, W_fuel, CST(1:12)]
b = x(1);           % Total wingspan (m)
c_r = x(2);         % Root chord (m)
c_k = x(3);         % Kink chord (m)
c_t = x(4);         % Tip chord (m)
M_cr = x(5);        % Cruise Mach number [-]
h_cr = x(6);        % Cruise altitude (m)
W_fuel = x(7);      % Fuel weight (N)
CST = x(8:19);      % CST airfoil parameters (12 values: 6 upper + 6 lower)
CST = reshape(CST, 1, []); % Defensive: always row vector

%%  SECTION 1: CONSTANT VALUES
b_k = 4.36 + 3.95/2;          % Spanwise location of kink (m) [Estimated, drawing]

% Wing angles (in degrees)
sweep_te_k = 0.01;    % Trailing edge sweep angle at kink (deg) [Estimated for aerodynamic efficiency]
dihedral = 5;       % Dihedral angle (deg) [Typical A320 value]

% Twist (Washout) - CRITICAL FOR CONVERGENCE
% Tips must be twisted down to prevent shock-induced tip stall at Mach 0.78
twist_r = 0;        % Twist at root (deg) [Standard]
twist_k = -1.0;     % Twist at kink (deg) [Typical washout]
twist_t = -3;       % Twist at tip (deg) [Typical washout]

M_MO_ref = 0.82;   % Maximum operating Mach number [-] [Specification]
n_max = 2.5;        % Maximum load factor (Structural) [CS-25 requirement]

gamma = 1.4;        % Ratio of specific heats [-] [ISA]
R = 287.058;        % Specific gas constant (J/kg-K) [ISA]

W_AminusW = 450000; % Aircraft weight minus wing (N) [Estimated for A320 class]

% Spar and tank limits
spar_locs = [0.2, 0.6];  % Front and rear spar locations [% chord, typical wing box]
tank_limits = [0, 0.85];  % Fuel tank spanwise limits [% half-span, typical]

% Engine data
engine_data.count = 1;      % Engines per half-wing [Standard twin-engine config]
engine_data.y_location = 4.7;  % Engine spanwise position (m) [~35% span, typical]
engine_data.weight = 1969;  % Engine weight (kg) [CFM56-5B class turbofan]

% Materials
m_up = [7.1e10, 2795, 4.8e8, 4.6e8];
m_fr = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.upper = m_up;
mat_props.lower = m_up;
mat_props.front = m_fr;
mat_props.rear  = m_up;

C_T_ref = 1.8639e-4;    % Reference specific fuel consumption (1/s)

%% SECTION 3: INITIAL GUESS VALUES FOR WEIGHTS

% Weight breakdown (in Newtons)
global W_wing
W_wing = 6344*9.81;     % Wing weight (N) - initial guess (~5100 kg) [Typical for A320]
% W_fuel is now a design variable extracted from x(7)

%%  SECTION 4: START OF OPTIMIZER LOOP

% Structures parameters
MTOW = (W_AminusW + 2*W_wing + W_fuel) / 9.81; % kg
ZFW = (W_AminusW + 2*W_wing) / 9.81;           % kg
S_ref = (c_r + c_k) * b_k + (c_k + c_t) * (b/2 - b_k); % Area

% Airfoils
cst2dat(CST, 'EMWET/optimized_airfoil');  % Generate .dat file for airfoil in emwet folder
airfoils.root = 'EMWET/optimized_airfoil';
airfoils.kink = 'EMWET/optimized_airfoil';
airfoils.tip  = 'EMWET/optimized_airfoil';

% Calculate cruise speed and atmospheric conditions using ISA
[~, ~, T_cr] = stdatm(h_cr);  % Get temperature at cruise altitude (K)
a_cr = sqrt(gamma * R * T_cr);  % Speed of sound at cruise altitude (m/s)
V_cr = M_cr * a_cr;    % Cruise speed (m/s)
V_MO_ref = M_MO_ref * a_cr;  % Maximum operating speed (m/s)

% Engine reference conditions
V_cr_ref = 230.1563;        % Reference cruise speed (m/s) - matches actual cruise
h_cr_ref = 1.1278e+04;        % Reference cruise altitude (m)

% % Plot wing geometry for verification
% plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, dihedral);
% plot_airfoil(CST);

%% SECTION 5: MDA LOOP FOR WING WEIGHT CONVERGENCE

% Convergence parameters
max_iter = 20;
tol = 0.01;  % 1% convergence tolerance
converged = false;
iter = 0;

while ~converged && iter < max_iter
    iter = iter + 1;
    
    W_wing_old = W_wing;
    
    % Loads analysis - Use actual loads() function
    [Y, L, M] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                      n_max, V_MO_ref, W_AminusW, h_cr,...
                      b, c_r, c_k, c_t, CST, W_wing, W_fuel);

    % Structures analysis
    W_wing_new = 9.81 * structures(Y, L, M, ...
                            sweep_te_k, b_k, dihedral, ...
                            b, c_r, c_k, c_t, ...
                            MTOW, ZFW, 1, S_ref, ...
                            mat_props, spar_locs, ...
                            tank_limits, engine_data, airfoils, ...
                            'a320_main');
    
    err_wing = abs(W_wing_new - W_wing_old) / W_wing_old;
        
    if err_wing < tol
        converged = true;
    end
    
    % Update for next iteration
    W_wing = W_wing_new;
    MTOW = (W_AminusW + 2*W_wing + W_fuel) / 9.81;
    ZFW = (W_AminusW + 2*W_wing) / 9.81;
end

%% SECTION 6: FINAL PERFORMANCE CALCULATION

% Clear any file handles and pause briefly to allow Windows to release file locks
% This is needed because Q3D_solver creates temporary files that may not be
% immediately released on Windows systems
pause(0.1);

% Aerodynamic analysis
[Cl, Cd] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                M_cr, W_AminusW, h_cr,...
                b, c_r, c_k, c_t, CST, W_wing, W_fuel);

% Add drag contribution from aircraft minus wing (fuselage, tail, nacelles, etc.)
% Typical A320 zero-lift drag breakdown: wing ~40%, rest ~60%
% Estimate Cd0 contribution from non-wing components
Cd0_rest = 0.012;  % Parasite drag coefficient for fuselage, tail, nacelles, etc.
Cd = Cd + Cd0_rest;

% Calculate L/D ratio from aerodynamic coefficients
if isnan(Cd) || Cd <= 0
    L_D_ratio = 0;
else
    L_D_ratio = Cl / Cd;
end

% Define weights for cruise segment
W_TO_max = W_AminusW + 2 * W_wing + W_fuel;
W_start_cr = W_TO_max;  % Weight at start of cruise (after taxi, takeoff, climb)
W_end_cr = (1 - W_fuel / W_TO_max) * W_start_cr / (0.938);

% Calculate performance using Breguet equations
[Range, ~, ~] = calculatePerformance(V_cr, h_cr, L_D_ratio, ...
    W_start_cr, W_end_cr, W_TO_max, ...
    V_cr_ref, h_cr_ref, C_T_ref);

fprintf('b=%.2f c_r=%.2f c_k=%.2f c_t=%.2f M=%.3f h=%.0f W_f=%.0f\n', b, c_r, c_k, c_t, M_cr, h_cr, W_fuel/9.81);
fprintf('R=%.0fkm W_w=%.0fkg L/D=%.2f MTOW=%.0fkg\n', Range/1000, W_wing/9.81, L_D_ratio, MTOW);

end