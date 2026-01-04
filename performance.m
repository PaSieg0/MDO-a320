function W_fuel = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits)
% This function calculates the fuel weight based on available tank volume.
% Currently uses calculation from calculate_wing_tank_volume.m
%
% Inputs:
%   b           - Total wingspan (m)
%   c_r         - Root chord (m)
%   c_k         - Kink chord (m)
%   c_t         - Tip chord (m)
%   b_k         - Spanwise location of kink (m)
%   spar_locs   - [front_spar, rear_spar] as fraction of chord
%   tank_limits - [start_eta, end_eta] as fraction of half-span
%
% Outputs:
%   W_fuel      - Fuel weight (N)

% Fuel density (kg/m^3)
rho_fuel = 804; % Jet A-1 fuel density

% Airfoil data file
filename = 'b737a.dat';

% Apply tank spanwise limits
half_span = b / 2;
y_tank_start = tank_limits(1) * half_span;  % Start of fuel tank (m)
y_tank_end = tank_limits(2) * half_span;    % End of fuel tank (m)

% Calculate volume using new function
[Vol_full_one_wing, ~, ~] = calculate_wing_tank_volume(c_r, c_k, c_t, b_k, b, filename, spar_locs(1), spar_locs(2));

% Scale by spanwise fraction
span_fraction = (y_tank_end - y_tank_start) / half_span;
V_total = Vol_full_one_wing * 2 * span_fraction;

% Calculate fuel mass and weight
m_fuel = rho_fuel * V_total * 0.85; % 0.85 = volume utilization factor
W_fuel = m_fuel * 9.81; % Convert to Newtons

end