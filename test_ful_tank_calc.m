%% Main Script to Calculate Wing Tank Volume
clc; clear; close all;

% --- 1. SETUP PARAMETERS ---
airfoil_file = 'b737a.dat';  % Ensure this file is in your folder

% Tank Limits (% of chord)
limit_start = 0.15;
limit_end   = 0.85;

% Wing Geometry (Example values - REPLACE with your actual data)
c_r  = 6.0;   % Root chord (m)
c_k  = 3.8;   % Kink chord (m)
c_t  = 1.5;   % Tip chord (m)
b_k   = 5.5;   % Span to kink (distance from root to kink) (m)
b  = 34.1;  % Total wingspan (tip-to-tip) (m)


% --- 2. CALCULATE WING TANK VOLUME ---
fprintf('\nCalculating Wing Volume based on geometry...\n');

% Call the function
try
    [Vol_Total, Vol_Inner, Vol_Outer] = calculate_wing_tank_volume( ...
        c_r, c_k, c_t, ...
        b_k, b, ...
        airfoil_file, ...
        limit_start, limit_end);
catch ME
    fprintf('Error calculating wing tank volume: %s\n', ME.message);
    return;
end


% --- 4. DISPLAY FINAL RESULTS ---
fprintf('--------------------------------------\n');
fprintf('       WING TANK VOLUME REPORT        \n');
fprintf('--------------------------------------\n');
fprintf('Tank Limits:         %.0f%% to %.0f%% chord\n', limit_start*100, limit_end*100);
fprintf('--------------------------------------\n');
fprintf('Inner Tank Volume:   %.4f m^3\n', Vol_Inner);
fprintf('Outer Tank Volume:   %.4f m^3\n', Vol_Outer);
fprintf('Total Volume (1 Wing): %.4f m^3\n', Vol_Total);
fprintf('--------------------------------------\n');
fprintf('TOTAL FUEL VOLUME:   %.4f m^3 (Both Wings)\n', Vol_Total * 2);
fprintf('--------------------------------------\n');