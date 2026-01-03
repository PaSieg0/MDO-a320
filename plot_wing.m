function plot_wing(b, c_r, c_k, c_t, b_k, sweep_te_k, ~)
% PLOT_WING Plots a top view of the right wing planform
%
% Inputs:
%   b           - Total wingspan (m)
%   c_r         - Root chord (m)
%   c_k         - Kink chord (m)
%   c_t         - Tip chord (m)
%   b_k         - Spanwise location of kink (m)
%   sweep_te_k  - Trailing edge sweep angle at kink (deg)
%   dihedral    - Dihedral angle (deg) - not used in top view
%
% Output:
%   Figure showing top view of right wing planform

% Root Station (y=0)
x_le_r = 0;
x_te_r = c_r;

% Kink Station (y=b_k)
x_te_k = x_te_r + b_k * tand(sweep_te_k);
x_le_k = x_te_k - c_k;

% Tip Station (y=b/2)
% Assume constant leading edge sweep from root to tip
tan_sweep_le = (x_le_k - x_le_r) / b_k;
x_le_t = x_le_r + (b/2) * tan_sweep_le;
x_te_t = x_le_t + c_t;

% Define wing outline coordinates (clockwise from root LE)
x_wing = [x_le_r, x_le_k, x_le_t, x_te_t, x_te_k, x_te_r, x_le_r];
y_wing = [0, b_k, b/2, b/2, b_k, 0, 0];

% Create figure
figure('Name', 'Wing Planform - Top View', 'NumberTitle', 'off');
hold on;
grid on;
axis equal;

% Plot wing outline
fill(x_wing, y_wing, [0.7 0.85 1], 'EdgeColor', [0 0.4 0.7], 'LineWidth', 2);

% Plot reference lines
plot([x_le_r, x_le_k, x_le_t], [0, b_k, b/2], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Leading Edge');
plot([x_te_r, x_te_k, x_te_t], [0, b_k, b/2], 'b--', 'LineWidth', 1.5, 'DisplayName', 'Trailing Edge');

% Add markers at key stations
plot([x_le_r, x_te_r], [0, 0], 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot([x_le_k, x_te_k], [b_k, b_k], 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot([x_le_t, x_te_t], [b/2, b/2], 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Labels
xlabel('Chordwise Position (m)', 'FontSize', 12);
ylabel('Spanwise Position (m)', 'FontSize', 12);
title('Right Wing Planform - Top View', 'FontSize', 14, 'FontWeight', 'bold');

% Add annotations
text(c_r/2, -0.5, sprintf('Root: c_r=%.2fm', c_r), 'HorizontalAlignment', 'center');
text(x_le_k + c_k/2, b_k + 0.5, sprintf('Kink: c_k=%.2fm', c_k), 'HorizontalAlignment', 'center');
text(x_le_t + c_t/2, b/2 + 0.5, sprintf('Tip: c_t=%.2fm', c_t), 'HorizontalAlignment', 'center');

legend('Location', 'best');

end