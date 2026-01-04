function W_fuel = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits)
% This function calculates the fuel weight based on available tank volume.
% Currently uses a simplified placeholder calculation.
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

% Apply tank spanwise limits
half_span = b / 2;
y_tank_start = tank_limits(1) * half_span;  % Start of fuel tank (m)
y_tank_end = tank_limits(2) * half_span;    % End of fuel tank (m)

% Tank chordwise depth fraction
tank_depth_fraction = spar_locs(2) - spar_locs(1);

% Height-to-chord ratio (assumed constant, typical for wing box)
h_to_c = 0.10;

% Integrate volume: V = integral of (c^2 * h_to_c * tank_depth_fraction) dy
% Chord varies linearly in each segment

% Initialize volume
V_total = 0;

% Inner segment: root to kink
if y_tank_end > 0 && y_tank_start < b_k
    y1 = max(y_tank_start, 0);
    y2 = min(y_tank_end, b_k);
    if y2 > y1
        % Linear chord variation: c(y) = c_r + (c_k - c_r) * y / b_k
        % Integral of c(y)^2 from y1 to y2:
        % = integral[(c_r + (c_k - c_r) * y / b_k)^2] dy
        % = [c_r^2 * y + 2*c_r*(c_k-c_r)/(2*b_k) * y^2 + (c_k-c_r)^2/(3*b_k^2) * y^3]
        V_inner = (c_r^2 * (y2 - y1) + ...
                   c_r * (c_k - c_r) / b_k * (y2^2 - y1^2) + ...
                   (c_k - c_r)^2 / (3 * b_k^2) * (y2^3 - y1^3)) * ...
                  h_to_c * tank_depth_fraction;
        V_total = V_total + V_inner;
    end
end

% Outer segment: kink to tip
if y_tank_end > b_k && y_tank_start < half_span
    y1 = max(y_tank_start, b_k);
    y2 = min(y_tank_end, half_span);
    if y2 > y1
        % Linear chord variation: c(y) = c_k + (c_t - c_k) * (y - b_k) / (half_span - b_k)
        % Shift coordinate: z = y - b_k, dz = dy, span_outer = half_span - b_k
        z1 = y1 - b_k;
        z2 = y2 - b_k;
        span_outer = half_span - b_k;
        
        V_outer = (c_k^2 * (z2 - z1) + ...
                   c_k * (c_t - c_k) / span_outer * (z2^2 - z1^2) + ...
                   (c_t - c_k)^2 / (3 * span_outer^2) * (z2^3 - z1^3)) * ...
                  h_to_c * tank_depth_fraction;
        V_total = V_total + V_outer;
    end
end

% Account for both half-wings
V_total = 2 * V_total;

% Calculate fuel mass and weight
m_fuel = rho_fuel * V_total * 0.85; % 0.85 = volume utilization factor
W_fuel = m_fuel * 9.81; % Convert to Newtons

end