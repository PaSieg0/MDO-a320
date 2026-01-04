function [P, rho, T] = stdatm(h)
% STDATM Standard Atmosphere model (ISA)
%
% Calculates atmospheric properties at a given altitude according to the
% International Standard Atmosphere (ISA) model.
%
% INPUTS:
%   h   - Geometric altitude (m), can be scalar or array
%
% OUTPUTS:
%   P   - Static pressure (Pa)
%   rho - Air density (kg/m^3)
%   T   - Static temperature (K)
%
% Valid for altitudes from sea level to 47,000 m
%
% Reference: ISO 2533:1975, ICAO Standard Atmosphere

    % Constants
    g0 = 9.80665;       % Standard gravity (m/s^2)
    R = 287.058;        % Specific gas constant for air (J/kg-K)
    
    % Sea level conditions
    T0 = 288.15;        % Sea level temperature (K)
    P0 = 101325;        % Sea level pressure (Pa)
    
    % Initialize outputs (same size as input)
    P = zeros(size(h));
    rho = zeros(size(h));
    T = zeros(size(h));
    
    % Process each altitude value
    for i = 1:numel(h)
        alt = h(i);
        
        if alt < 0
            error('Altitude must be non-negative');
        elseif alt <= 11000
            % Troposphere (0 to 11 km)
            L = -0.0065;            % Temperature lapse rate (K/m)
            T(i) = T0 + L * alt;
            P(i) = P0 * (T(i) / T0)^(-g0 / (L * R));
            rho(i) = P(i) / (R * T(i));
            
        elseif alt <= 20000
            % Lower Stratosphere (11 to 20 km) - Isothermal layer
            T_11 = 216.65;          % Temperature at 11 km (K)
            P_11 = 22632.1;         % Pressure at 11 km (Pa)
            
            T(i) = T_11;
            P(i) = P_11 * exp(-g0 * (alt - 11000) / (R * T_11));
            rho(i) = P(i) / (R * T(i));
            
        elseif alt <= 32000
            % Upper Stratosphere (20 to 32 km)
            L = 0.001;              % Temperature lapse rate (K/m)
            T_20 = 216.65;          % Temperature at 20 km (K)
            P_20 = 5474.89;         % Pressure at 20 km (Pa)
            
            T(i) = T_20 + L * (alt - 20000);
            P(i) = P_20 * (T(i) / T_20)^(-g0 / (L * R));
            rho(i) = P(i) / (R * T(i));
            
        elseif alt <= 47000
            % Middle Stratosphere (32 to 47 km)
            L = 0.0028;             % Temperature lapse rate (K/m)
            T_32 = 228.65;          % Temperature at 32 km (K)
            P_32 = 868.019;         % Pressure at 32 km (Pa)
            
            T(i) = T_32 + L * (alt - 32000);
            P(i) = P_32 * (T(i) / T_32)^(-g0 / (L * R));
            rho(i) = P(i) / (R * T(i));
            
        else
            error('Altitude exceeds valid range (0 to 47,000 m)');
        end
    end
    
    % Reshape outputs to match input shape
    P = reshape(P, size(h));
    rho = reshape(rho, size(h));
    T = reshape(T, size(h));
end
