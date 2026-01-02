function [R, W_fuel, C_T] = calculatePerformance(V, h, L_D_ratio, W_start_cr, W_end_cr, W_TO_max, V_cr_ref, h_cr_ref, C_T_ref)
% calculatePerformance - Computes aircraft range and fuel weight based on Breguet equations.
%
% This function implements the performance calculations as specified in the
% "AE4-205 MDO for Aerospace Applications 2025/2026" homework assignment.
%
% SYNTAX:
%   [R, W_fuel, C_T] = calculatePerformance(V, h, L_D_ratio, W_start_cr, W_end_cr, W_TO_max, V_cr_ref, h_cr_ref, C_T_ref)
%
% INPUTS:
%   V           - Actual cruise speed (m/s)
%   h           - Actual cruise altitude (m)
%   L_D_ratio   - Lift-to-drag ratio at cruise conditions (-)
%   W_start_cr  - Aircraft weight at the start of the cruise phase (N)
%   W_end_cr    - Aircraft weight at the end of the cruise phase (N)
%   W_TO_max    - Maximum Take-Off Weight (N)
%   V_cr_ref    - Reference cruise speed for engine sizing (m/s)
%   h_cr_ref    - Reference cruise altitude for engine sizing (m)
%   C_T_ref     - Reference specific fuel consumption (N/Ns or 1/s)
%
% OUTPUTS:
%   R           - Mission range (m)
%   W_fuel      - Total fuel weight required for the mission (N)
%   C_T         - Actual specific fuel consumption at cruise conditions (N/Ns or 1/s)
%
% See also: AE4-205 MDO Homework Assignment, Page 5-6.

%% 1. Calculate Propulsive Efficiency Factor (eta)
% This factor accounts for performance degradation of the engine when
% operating outside the reference flight condition.
% Source: Page 5 of the assignment document.
eta = exp(-((V - V_cr_ref)^2 / 2.70^2) - ((h - h_cr_ref)^2 / 22500^2));

%% 2. Calculate Actual Specific Fuel Consumption (C_T)
% Source: Page 5 of the assignment document.
C_T = C_T_ref / eta;

%% 3. Calculate Mission Range (R) using Breguet Equation
% Source: Page 5 of the assignment document.
% The primary objective function to be maximized.
% Note: The range R will be in meters.
R = (V / C_T) * L_D_ratio * log(W_start_cr / W_end_cr);

%% 4. Calculate Total Fuel Weight (W_fuel)
% This formula estimates the total mission fuel from the cruise weight
% fraction and maximum take-off weight.
% The factor 0.938 accounts for fuel fractions for mission segments 
% outside the cruise stage (e.g., taxi, take-off, climb, descent).
% Source: Page 5 of the assignment document.
W_fuel = (1 - 0.938 * (W_end_cr / W_start_cr)) * W_TO_max;

end