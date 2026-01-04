function [Range] = optimize(x)
% OPTIMIZE Wrapper function to run MDO analysis for given design variables
% Extract design variables from input vector x
b = x(1);           % Total wingspan (m)
c_r = x(2);         % Root chord (m)
c_k = x(3);         % Kink chord (m)
c_t = x(4);         % Tip chord (m)
M_cr = x(5);        % Cruise Mach number [-]
h_cr = x(6);        % Cruise altitude (m)
CST = x(7:16);      % CST airfoil parameters (10 values)

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

W_AminusW = 400000; % Aircraft weight minus wing (N) [Estimated for A320 class]

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
W_wing = 60000;     % Wing weight (N) - initial guess (~5100 kg) [Typical for A320]
W_fuel = 150000;    % Fuel weight (N) [Estimated for medium range]

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

% Engine reference conditions NEEDS CHECK
V_cr_ref = V_cr;        % Reference cruise speed (m/s) - matches actual cruise
h_cr_ref = h_cr;        % Reference cruise altitude (m)

% % Plot wing geometry for verification
% plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, dihedral);
% plot_airfoil(CST);

% Evaluate fuel tank volume (performance)
W_fuel = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);
fprintf('Initial fuel capacity: %.2f kg\n', W_fuel / 9.81);

%% SECTION 5: MDA LOOP FOR WING WEIGHT CONVERGENCE

% Convergence parameters
max_iter = 20;
tol = 0.001;  % 0.1% convergence tolerance
converged = false;
iter = 0;

fprintf('\n========== STARTING MDO CONVERGENCE LOOP ==========\n');

while ~converged && iter < max_iter
    iter = iter + 1;
    fprintf('\n--- Iteration %d ---\n', iter);
    
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
    
    fprintf('  W_wing: %.2f kg (change: %.2f%%)\n', W_wing_new/9.81, err_wing*100);
    
    if err_wing < tol
        converged = true;
        fprintf('>>> CONVERGED after %d iterations <<<\n', iter);
    end
    
    % Update for next iteration
    W_wing = W_wing_new;
    MTOW = (W_AminusW + 2*W_wing + W_fuel) / 9.81;
    ZFW = (W_AminusW + 2*W_wing) / 9.81;
end

if ~converged
    fprintf('>>> WARNING: Did not converge after %d iterations <<<\n', max_iter);
end

fprintf('===================================================\n\n');

%% SECTION 6: FINAL PERFORMANCE CALCULATION

% Aerodynamic analysis
[Cl, Cd] = aero(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                M_cr, W_AminusW, h_cr,...
                b, c_r, c_k, c_t, CST, W_wing_new, W_fuel);

% Calculate L/D ratio from aerodynamic coefficients
L_D_ratio = Cl / Cd;

% Define weights for cruise segment
W_TO_max = W_AminusW + 2 * W_wing_new + W_fuel;
W_start_cr = W_TO_max;  % Weight at start of cruise (after taxi, takeoff, climb)
W_end_cr = (1 - W_fuel / W_TO_max) * W_start_cr / (0.938);

% Calculate performance using Breguet equations
[Range, ~, ~] = calculatePerformance(V_cr, h_cr, L_D_ratio, ...
    W_start_cr, W_end_cr, W_TO_max, ...
    V_cr_ref, h_cr_ref, C_T_ref);

end