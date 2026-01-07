function plot_airfoil(CST)
% PLOT_AIRFOIL Plots an airfoil from CST (Class Shape Transformation) coefficients
%
% Inputs:
%   CST        - Vector of CST coefficients [Upper_A0, ..., Upper_An, Lower_A0, ..., Lower_An]
%                First half: upper surface coefficients
%                Second half: lower surface coefficients
%
% Output:
%   Figure showing the airfoil shape

% Number of points for smooth curve
n_points = 200;

% Generate x coordinates (0 to 1)
x = linspace(0, 1, n_points);

% Split CST into upper and lower surface coefficients
n_total = length(CST);
n_half = n_total / 2;
CST_upper = CST(1:n_half);
CST_lower = CST(n_half+1:end);

% Use D_airfoil2 to get airfoil coordinates
Au = CST_upper;  % Upper surface Bernstein parameters
Al = CST_lower;  % Lower surface Bernstein parameters
X = x';          % X-coordinates as column vector

[Xtu, Xtl, ~, ~, ~, ~] = D_airfoil2(Au, Al, X);

% Extract y coordinates from D_airfoil2 output
y_upper = Xtu(:, 2)';
y_lower = Xtl(:, 2)';

% Create figure
figure('Name', 'Airfoil Shape', 'NumberTitle', 'off');
hold on;
grid on;
axis equal;

% Plot airfoil
plot(x, y_upper, 'b-', 'LineWidth', 2, 'DisplayName', 'Upper Surface');
plot(x, y_lower, 'r-', 'LineWidth', 2, 'DisplayName', 'Lower Surface');

% Add chord line
plot([0, 1], [0, 0], 'k--', 'LineWidth', 1);

% Formatting
xlabel('x/c', 'FontSize', 12);
ylabel('y/c', 'FontSize', 12);
xlim([-0.05, 1.05]);
legend('Location', 'best');

% Add leading and trailing edge markers
plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(1, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

hold off;

end