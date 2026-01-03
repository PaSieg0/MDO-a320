% Test sensitivity of loads to wing weight
clear; clc;

fprintf('========== TESTING LOAD SENSITIVITY TO W_WING ==========\n\n');

%% Parameters from main.m
b = 34.1;
b_k = 5.5;
c_r = 7.2;
c_k = 4.5;
c_t = 2.0;
sweep_te_k = 3.67;
dihedral = 0.48;
twist_r = 0;
twist_k = -1.0;
twist_t = -3;
n_max = 2.5;
V_MO_ref = 250;
W_AminusW = 400000;
h_cr = 10000;
W_fuel = 150000;
CST = [0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797];

%% Test different wing weights
W_wing_values = [10000, 30000, 50000, 70000, 90000];  % N
results = zeros(length(W_wing_values), 4);  % Store: mean L, max L, mean M, total W

fprintf('Testing loads() with different wing weights...\n\n');

for i = 1:length(W_wing_values)
    W_wing = W_wing_values(i);
    
    [Y, L, M] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,...
                      n_max, V_MO_ref, W_AminusW, h_cr,...
                      b, c_r, c_k, c_t, CST, W_wing, W_fuel);
    
    W_total = W_AminusW + W_wing + W_fuel;
    results(i, :) = [mean(L), max(L), mean(M), W_total];
    
    fprintf('Test %d: W_wing = %.2f kg (%.0f N)\n', i, W_wing/9.81, W_wing);
    fprintf('  W_total = %.2f kg\n', W_total/9.81);
    fprintf('  Mean L  = %.2e N/m\n', mean(L));
    fprintf('  Max L   = %.2e N/m\n', max(L));
    fprintf('  Mean M  = %.2e Nm/m\n\n', mean(M));
end

%% Analysis
fprintf('========== SENSITIVITY ANALYSIS ==========\n');
fprintf('W_wing (kg) | W_total (kg) | Mean L (N/m) | Change in L\n');
fprintf('------------------------------------------------------------\n');
for i = 1:length(W_wing_values)
    if i == 1
        change_str = '(baseline)';
    else
        pct_change = 100 * (results(i,1) - results(1,1)) / results(1,1);
        change_str = sprintf('%+.2f%%', pct_change);
    end
    fprintf(' %9.2f  |  %10.2f  |  %10.2e  | %s\n', ...
            W_wing_values(i)/9.81, results(i,4)/9.81, results(i,1), change_str);
end

%% Check proportionality
fprintf('\n========== PROPORTIONALITY CHECK ==========\n');
W_total_ratio = results(:,4) / results(1,4);
L_ratio = results(:,1) / results(1,1);

fprintf('If loads scale with total weight, these ratios should match:\n');
fprintf('W_wing (kg) | W_total ratio | Mean L ratio | Match?\n');
fprintf('--------------------------------------------------------\n');
for i = 1:length(W_wing_values)
    match = abs(W_total_ratio(i) - L_ratio(i)) < 0.01;
    match_str = 'YES';
    if ~match
        match_str = 'NO';
    end
    fprintf(' %9.2f  |     %.4f     |    %.4f     | %s\n', ...
            W_wing_values(i)/9.81, W_total_ratio(i), L_ratio(i), match_str);
end

fprintf('\n========== CONCLUSION ==========\n');
if all(abs(W_total_ratio - L_ratio) < 0.01)
    fprintf('✓ Loads scale correctly with total aircraft weight.\n');
    fprintf('  The loads() function properly couples wing weight changes.\n');
    fprintf('  Mean L changes by %.2f%% when W_wing changes by %.2f%%\n', ...
            100*(L_ratio(end)-1), 100*(W_wing_values(end)/W_wing_values(1)-1));
else
    fprintf('✗ Loads do NOT scale proportionally with total weight.\n');
    fprintf('  This indicates an issue with the loads() coupling.\n');
end
