function [Y, L, M] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,... 
                           n_max, V_MO_ref, W_AminusW, h_cr,...
                           b, c_r, c_k, c_t, CST, W_wing, W_fuel)
% This function calculates the spanwise load distribution on a wing.
% It uses the Q3D_solver based on the provided geometric and flight parameters.

% Define the AC (Aircraft) structure for the Q3D_solver

% --- Wing Planform Geometry Calculation ---
% The geometry is defined by creating a matrix with the properties of each
% station (root, kink, tip).
% Columns are: [x_LE, y, z, chord, twist_angle]

% Root Station (y=0) is at the origin
x_le_r = 0;
x_te_r = c_r; 

% Kink Station (y=b_k)
x_te_k = x_te_r + b_k * tan(deg2rad(sweep_te_k));
x_le_k = x_te_k - c_k;

% Tip Station (y=b/2)
% The leading edge sweep is assumed to be constant from root to tip.
tan_sweep_le = (x_le_k - x_le_r) / b_k;
x_le_t = x_le_r + (b/2) * tan_sweep_le;

% Z-Coordinates based on dihedral
z_k = b_k * tan(deg2rad(dihedral));
z_t = (b/2) * tan(deg2rad(dihedral));

% Assemble the Geometry Matrix
%                x_LE (m)    y (m)     z (m)    chord(m)    twist angle (deg) 
AC.Wing.Geom = [x_le_r,    0,        0,       c_r,        twist_r;
                x_le_k,    b_k,      z_k,     c_k,        twist_k;
                x_le_t,    b/2,      z_t,     c_t,        twist_t];

% --- Atmospheric Properties Calculation (ISA Model for Troposphere, h < 11km) ---
% Constants for International Standard Atmosphere
% Properties at cruise altitude h_cr using standard atmosphere model
[~, rho_cr, T_cr] = stdatm(h_cr);

% Additional atmospheric properties
R = 287.058;    % Specific gas constant for dry air [J/(kg*K)]
gamma = 1.4;    % Ratio of specific heats for air
a_cr = sqrt(gamma * R * T_cr); % Speed of sound at altitude [m/s]

% Reference conditions for viscosity calculation
T0 = 288.15;    % Sea level temperature [K]

% Dynamic viscosity using Sutherland's Law
mu_0 = 1.789e-5; % Reference viscosity at T0 [Pa.s]
S_suth = 110.4;  % Sutherland's constant [K]
mu_cr = mu_0 * (T0 + S_suth) / (T_cr + S_suth) * (T_cr / T0)^(1.5);

% --- Aerodynamic and Flight Condition Setup ---

% Wing incidence angle (degree)
AC.Wing.inc  = 0;   
            
% Airfoil coefficients from CST parameter
AC.Wing.Airfoils   = [CST;
                      CST;
                      CST];
AC.Wing.eta = [0;b_k/(b/2);1];  % Spanwise location of the airfoil sections

% Analysis settings
AC.Visc  = 0;              % 0 for inviscid (potential flow) analysis
AC.Aero.MaxIterIndex = 150;

% Flight Condition using calculated atmospheric properties
AC.Aero.V     = V_MO_ref;  % Flight speed (m/s)
AC.Aero.rho   = rho_cr;    % Air density at cruise altitude
AC.Aero.alt   = h_cr;      % Flight altitude (m)
AC.Aero.M     = V_MO_ref / a_cr; % Mach number at altitude

% Calculate Reynolds number based on Mean Aerodynamic Chord (MAC)
S_half = (c_r + c_k) * b_k / 2 + (c_k + c_t) * (b/2 - b_k) / 2; % Semi-span area
int_c2_dy = (b_k/3) * (c_r^2 + c_r*c_k + c_k^2) + ...
            ((b/2 - b_k)/3) * (c_k^2 + c_k*c_t + c_t^2);
MAC = (1 / S_half) * int_c2_dy;
AC.Aero.Re = (rho_cr * V_MO_ref * MAC) / mu_cr;

% Calculate the required lift coefficient (CL) for the flight condition
W_total = W_AminusW + W_wing + W_fuel; % Total weight (N)
S_wing = 2 * S_half; % Total wing planform area
q = 0.5 * rho_cr * AC.Aero.V^2; % Dynamic pressure
Required_Lift = n_max * W_total;   % Required lift force (Weight * load factor)
AC.Aero.CL = Required_Lift / (q * S_wing);

% Handle directory navigation for Q3D
originalDir = pwd;
loadsPath = fileparts(mfilename('fullpath'));
q3dPath = fullfile(loadsPath, 'Q3D');

% Ensure we're in the Q3D directory before calling the solver
if exist(q3dPath, 'dir')
    cd(q3dPath);
else
    error('Q3D path does not exist: %s', q3dPath);
end

% Verify Storage folder exists
storagePath = fullfile(q3dPath, 'Storage');
if ~exist(storagePath, 'dir')
    mkdir(storagePath);
end

%% Run the aerodynamic solver
try
    % Use evalc to suppress command-line output from the solver
    [~, Res] = evalc('Q3D_solver(AC)');
catch ME
    cd(originalDir);
    rethrow(ME);
end

% Return to original directory
cd(originalDir);

%% Extract and convert results to forces and moments
% Get the distribution of chord length at each station from the results
c_dist = Res.Wing.chord; 

% Y: Spanwise stations (m)
Y = Res.Wing.Yst;

% L: Sectional lift force per unit span (N/m)
% L' = q * c * cl
L = q .* c_dist .* Res.Wing.cl;

% M: Sectional pitching moment per unit span about quarter-chord (Nm/m)
% M' = q * c^2 * cm_c/4
M = q .* (c_dist.^2) .* Res.Wing.cm_c4;

Y = [Y; b/2];
L = [L; 0];
M = [M; 0];

end