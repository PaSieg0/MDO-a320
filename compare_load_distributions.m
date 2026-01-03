% Compare load distributions: EMWET-working vs Q3D loads()
clear; clc;

fprintf('========== COMPARING LOAD DISTRIBUTIONS ==========\n\n');

%% Common parameters
b = 34.1;
b_k = 5.5;
c_r = 7.2;
c_k = 4.5;
c_t = 2.0;
sweep_te_k = 10;
dihedral = 5;
twist_r = 0;
twist_k = -1.0;
twist_t = -3;
n_max = 2.5;
V_MO_ref = 250;
W_AminusW = 400000;
h_cr = 10000;
W_wing = 50000;
W_fuel = 150000;
CST = [0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797];

%% 2. Generate Q3D loads() distribution
[Y_q3d, L_q3d, M_q3d] = loads(sweep_te_k, b_k, dihedral, twist_r, twist_k, twist_t,...
                              n_max, V_MO_ref, W_AminusW, h_cr,...
                              b, c_r, c_k, c_t, CST, W_wing, W_fuel);

% Add wingtip point at semispan (Y = b/2) with zero loads
Y_q3d = [Y_q3d; b/2];
L_q3d = [L_q3d; 0];
M_q3d = [M_q3d; 0];

fprintf('Q3D distribution:\n');
fprintf('  Points: %d\n', length(Y_q3d));
fprintf('  L range: %.2e to %.2e N/m\n', min(L_q3d), max(L_q3d));
fprintf('  M range: %.2e to %.2e Nm/m\n', min(M_q3d), max(M_q3d));
if any(M_q3d > 0)
    fprintf('  M sign changes: YES\n\n');
else
    fprintf('  M sign changes: NO\n\n');
end

%% 1. Generate EMWET-working distribution (simple elliptical)
W_total = W_AminusW + W_wing + W_fuel;
Y_simple = linspace(0, b/2, 15)';
eta_simple = 2*Y_simple/b;
L_simple = (W_total / (0.5 * pi * b/2)) * sqrt(1 - eta_simple.^2) * 2;  % Elliptical
M_simple = (-4.5e4 + 4.0e4 * eta_simple) * 2;  % Linear transition

fprintf('Simple elliptical distribution:\n');
fprintf('  Points: %d\n', length(Y_simple));
fprintf('  L range: %.2e to %.2e N/m\n', min(L_simple), max(L_simple));
fprintf('  M range: %.2e to %.2e Nm/m\n', min(M_simple), max(M_simple));
if any(M_simple > 0)
    fprintf('  M sign changes: YES\n\n');
else
    fprintf('  M sign changes: NO\n\n');
end

%% 3. Plot comparison
figure('Position', [100, 100, 1200, 800]);

% Lift distribution
subplot(2,2,1)
plot(Y_simple, L_simple, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Y_q3d, L_q3d, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('Spanwise Position Y (m)', 'FontSize', 12);
ylabel('Lift per unit span L (N/m)', 'FontSize', 12);
title('Lift Distribution Comparison', 'FontSize', 14, 'FontWeight', 'bold');
legend('Simple (works with EMWET)', 'Q3D loads()', 'Location', 'best');
xlim([0, b/2]);

% Moment distribution
subplot(2,2,2)
plot(Y_simple, M_simple, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Y_q3d, M_q3d, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
yline(0, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Spanwise Position Y (m)', 'FontSize', 12);
ylabel('Moment per unit span M (Nm/m)', 'FontSize', 12);
title('Moment Distribution Comparison', 'FontSize', 14, 'FontWeight', 'bold');
legend('Simple (works with EMWET)', 'Q3D loads()', 'Zero', 'Location', 'best');
xlim([0, b/2]);

% Normalized eta plots
subplot(2,2,3)
eta_q3d = 2*Y_q3d/b;
plot(eta_simple, L_simple, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(eta_q3d, L_q3d, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('Normalized Span \eta = 2y/b', 'FontSize', 12);
ylabel('Lift per unit span L (N/m)', 'FontSize', 12);
title('Lift vs Normalized Span', 'FontSize', 14, 'FontWeight', 'bold');
legend('Simple', 'Q3D', 'Location', 'best');
xlim([0, 1]);

subplot(2,2,4)
plot(eta_simple, M_simple, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(eta_q3d, M_q3d, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
yline(0, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('Normalized Span \eta = 2y/b', 'FontSize', 12);
ylabel('Moment per unit span M (Nm/m)', 'FontSize', 12);
title('Moment vs Normalized Span', 'FontSize', 14, 'FontWeight', 'bold');
legend('Simple', 'Q3D', 'Zero', 'Location', 'best');
xlim([0, 1]);

sgtitle('Load Distribution Comparison: Simple vs Q3D', 'FontSize', 16, 'FontWeight', 'bold');

%% 4. Statistical comparison
fprintf('========== STATISTICAL COMPARISON ==========\n');
fprintf('Lift (L):\n');
fprintf('  Simple mean: %.2e N/m\n', mean(L_simple));
fprintf('  Q3D mean:    %.2e N/m (%.1f%% difference)\n', ...
        mean(L_q3d), 100*abs(mean(L_q3d)-mean(L_simple))/mean(L_simple));

fprintf('\nMoment (M):\n');
fprintf('  Simple mean: %.2e Nm/m\n', mean(M_simple));
fprintf('  Q3D mean:    %.2e Nm/m\n', mean(M_q3d));
fprintf('  Simple min/max: %.2e to %.2e Nm/m\n', min(M_simple), max(M_simple));
fprintf('  Q3D min/max:    %.2e to %.2e Nm/m\n', min(M_q3d), max(M_q3d));

fprintf('\n========== KEY DIFFERENCE ==========\n');
if all(M_simple < 0) && all(M_q3d < 0)
    fprintf('Both distributions have all-negative moments\n');
elseif any(M_simple > 0) && all(M_q3d < 0)
    fprintf('>>> CRITICAL: Simple has sign transition, Q3D does not!\n');
    fprintf('    Simple M crosses zero at eta = %.3f\n', eta_simple(find(M_simple>0, 1)));
    fprintf('    This is why EMWET works with Simple but not Q3D\n');
elseif any(M_simple > 0) && any(M_q3d > 0)
    fprintf('Both distributions have sign transition\n');
end

%% 5. Test both distributions with EMWET
fprintf('\n========== TESTING WITH EMWET ==========\n');

% Common EMWET parameters
MTOW = (W_AminusW + W_wing + W_fuel) / 9.81;
ZFW = (W_AminusW + W_wing) / 9.81;
S_ref = 122.6;
mat_props.upper = [7.1e10, 2795, 4.8e8, 4.6e8];
mat_props.lower = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.front = [7.37765e10, 2795.68, 3.24065e8, 2.68905e8];
mat_props.rear = [7.1e10, 2795, 4.8e8, 4.6e8];
spar_locs = [0.15, 0.65];
tank_limits = [0.1, 0.9];
engine_data.count = 1;
engine_data.y_location = 4.7;
engine_data.weight = 1969;
airfoils.root = 'b737a';
airfoils.kink = 'b737a';
airfoils.tip = 'b737a';

fprintf('\nTest 1: Simple distribution (works with EMWET)\n');
W_wing_simple = 9.81 * structures(Y_simple, L_simple, M_simple, ...
                                   sweep_te_k, b_k, dihedral, ...
                                   b, c_r, c_k, c_t, ...
                                   MTOW, ZFW, n_max, S_ref, ...
                                   mat_props, spar_locs, ...
                                   tank_limits, engine_data, airfoils, ...
                                   'test_simple');
fprintf('  Result: W_wing = %.2f kg\n', W_wing_simple/9.81);

fprintf('\nTest 2: Q3D distribution\n');
W_wing_q3d = 9.81 * structures(Y_q3d, L_q3d, M_q3d, ...
                                sweep_te_k, b_k, dihedral, ...
                                b, c_r, c_k, c_t, ...
                                MTOW, ZFW, n_max, S_ref, ...
                                mat_props, spar_locs, ...
                                tank_limits, engine_data, airfoils, ...
                                'test_q3d');
fprintf('  Result: W_wing = %.2f kg\n', W_wing_q3d/9.81);

fprintf('\n========== EMWET RESULTS COMPARISON ==========\n');
fprintf('Simple distribution: %.2f kg\n', W_wing_simple/9.81);
fprintf('Q3D distribution:    %.2f kg\n', W_wing_q3d/9.81);
fprintf('Difference:          %.2f kg (%.1f%%)\n', ...
        abs(W_wing_simple-W_wing_q3d)/9.81, ...
        100*abs(W_wing_simple-W_wing_q3d)/W_wing_simple);

if abs(W_wing_simple - W_wing_q3d) < 1
    fprintf('\n>>> PROBLEM: Both give same weight (minimum gauge)\n');
    fprintf('    EMWET is not responding to load differences\n');
else
    fprintf('\n>>> SUCCESS: Different weights indicate EMWET is sizing properly\n');
end
