function [c,ceq] = constraints(x)
    % Extract design variables
    % x = [b, c_r, c_k, c_t, M_cr, h_cr, W_fuel, CST(1:12)]
        b = x(1);           % Total wingspan (m)
        c_r = x(2);         % Root chord (m)
        c_k = x(3);         % Kink chord (m)
        c_t = x(4);         % Tip chord (m)
        W_fuel = x(7);      % Fuel weight (N) - design variable
        CST = x(8:19);      % CST airfoil parameters
        
        global W_wing
        % Fallback if W_wing is not set
        if isempty(W_wing)
            W_wing = 6344 * 9.81;  % Default value
        end
        
        % Generate airfoil file for current CST (needed for fuel tank volume calc)
        cst2dat(CST, 'EMWET/optimized_airfoil');
        
        % Constant values
        b_k = 4.36 + 3.95/2;  % Spanwise location of kink (m)
        W_AminusW = 450000;   % Aircraft weight minus wing (N)
        
        % Calculate wing area for current design
        S_current = (c_r + c_k) * b_k + (c_k + c_t) * (b/2 - b_k);
        
        % Original design parameters
        b_orig = 34.0;
        c_r_orig = 7.0;
        c_k_orig = 3.7;
        c_t_orig = 1.6;
        
        % Calculate original wing area
        S_orig = (c_r_orig + c_k_orig) * b_k + (c_k_orig + c_t_orig) * (b_orig/2 - b_k);
        
        % Approximate total weight (using rough estimates)
        W_total = W_AminusW + 2*W_wing + W_fuel;  % N (aircraft - wing + 2*wing + fuel)
        
        % Wing loading constraint: WL_current <= WL_original
        % c1 > 0 means constraint violated
        WL_current = W_total / S_current;
        WL_orig = W_total / S_orig;
        c1 = WL_current - WL_orig;
        
        % Fuel weight constraint
        W_fuel_max = 16084.17 * 9.81; % Maximum fuel weight for original design (N)
        % Use performance function to calculate actual fuel capacity
        spar_locs = [0.2, 0.6];  % Front and rear spar locations
        tank_limits = [0, 0.85]; % Fuel tank spanwise limits
        
        % Calculate fuel capacity for current design
        fprintf('b=%.2f c_r=%.2f c_k=%.2f c_t=%.2f\n', b, c_r, c_k, c_t);
        W_fuel_max_current = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);
        
        % Fuel capacity constraint: W_fuel <= min(W_fuel_max_current, W_fuel_max)
        % c2 > 0 means constraint violated (fuel exceeds tank capacity)
        fprintf('Max Fuel Capacity Current Design: %.2f kg\n', W_fuel_max_current/9.81);
        c2 = W_fuel - min(W_fuel_max_current, W_fuel_max);

        % Inequality constraints
        c = [c1; c2];
        
        % No equality constraints
        ceq = [];
end
