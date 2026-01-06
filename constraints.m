function [c,ceq] = constraints(x)
    % Extract design variables
    % x = [b, c_r, c_k, c_t, M_cr, h_cr, W_fuel, CST(1:10)]
        b = x(1);           % Total wingspan (m)
        c_r = x(2);         % Root chord (m)
        c_k = x(3);         % Kink chord (m)
        c_t = x(4);         % Tip chord (m)
        W_fuel = x(7);      % Fuel weight (N) - design variable
        
        % Constant values
        b_k = 4.36 + 3.95/2;  % Spanwise location of kink (m)
        W_AminusW = 400000;   % Aircraft weight minus wing (N)
        
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
        W_total = W_AminusW + 2*6081.08*9.81 + 16275.44*9.81;  % N (aircraft - wing + 2*wing + fuel)
        
        % Wing loading constraint: WL_current <= WL_original
        % c1 > 0 means constraint violated
        WL_current = W_total / S_current;
        WL_orig = W_total / S_orig;
        c1 = WL_current - WL_orig;
        
        % Fuel weight constraint
        W_fuel_max = 16275.44 * 9.81; % Maximum fuel weight for original design (N)
        % Use performance function to calculate actual fuel capacity
        spar_locs = [0.2, 0.6];  % Front and rear spar locations
        tank_limits = [0, 0.85]; % Fuel tank spanwise limits
        
        % Calculate fuel capacity for current design
        W_fuel_max_current = performance(b, c_r, c_k, c_t, b_k, spar_locs, tank_limits);
        
        % Fuel capacity constraint: W_fuel == W_fuel_max_current (equality)
        % ceq = 0 means constraint satisfied
        ceq1 = W_fuel - min(W_fuel_max_current, W_fuel_max);

        % Inequality constraints
        c = [c1];
        
        % Equality constraints
        ceq = [ceq1];
end
