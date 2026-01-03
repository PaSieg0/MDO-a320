function W_fuel = performance(b, c_r, c_k, c_t, b_k, tank_limits)
% This function calculates the fuel weight based on available tank volume.
% Currently uses a simplified placeholder calculation.
%
% Inputs:
%   b           - Total wingspan (m)
%   c_r         - Root chord (m)
%   c_k         - Kink chord (m)
%   c_t         - Tip chord (m)
%   b_k         - Spanwise location of kink (m)
%   tank_limits - [front_limit, rear_limit] as fraction of chord
%
% Outputs:
%   W_fuel      - Fuel weight (N)

% Fuel density (kg/m^3)
rho_fuel = 804; % Jet A-1 fuel density

% Calculate approximate tank volume
% Simplified calculation: assume tank occupies a fraction of the wing box
tank_depth_fraction = tank_limits(2) - tank_limits(1); % Chordwise fraction

% Volume calculation for inner wing (root to kink)
c_avg_inner = (c_r + c_k) / 2;
V_inner = b_k * c_avg_inner * tank_depth_fraction * 0.15; % 0.15 = assumed height fraction

% Volume calculation for outer wing (kink to tip)
c_avg_outer = (c_k + c_t) / 2;
span_outer = b/2 - b_k;
V_outer = span_outer * c_avg_outer * tank_depth_fraction * 0.12; % 0.12 = thinner section

% Total volume (both half-wings)
V_total = 2 * (V_inner + V_outer);

% Calculate fuel mass and weight
m_fuel = rho_fuel * V_total * 0.85; % 0.85 = volume utilization factor
W_fuel = m_fuel * 9.81; % Convert to Newtons

end