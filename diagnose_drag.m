%% DIAGNOSTIC SCRIPT: Drag Calculation Debugger (Fixed Integration)
% Fixes: Calculates CD_ind manually from spanwise distribution.

clear; close all; clc;

% --- 1. Setup Paths ---
currentDir = fileparts(mfilename('fullpath'));
addpath(currentDir);
addpath(fullfile(currentDir, 'Q3D'));
addpath(fullfile(currentDir, 'EMWET'));

fprintf('=== STARTING DRAG DIAGNOSTIC ===\n');

% --- 2. Define Your Optimized Inputs ---
b = 33.48;
c_r = 6.98;
c_k = 4.21;
c_t = 1.46;
M_cr = 0.775;
h_cr = 11458;
% Your Optimized CST
CST = [0.1730, 0.1310, 0.2249, 0.0995, 0.3298, 0.2608, -0.2216, -0.1680, -0.1055, -0.3265, 0.1234, 0.3174];

% Estimated Weight
W_fuel = 16085 * 9.81; 
W_AminusW = 450000;
W_wing_est = 6344 * 9.81; 
W_total = W_AminusW + 2*W_wing_est + W_fuel;

% --- 3. Check Airfoil Shape ---
fprintf('Checking Airfoil Geometry...\n');
[xu, xl] = get_airfoil_coords(CST);
figure('Name', 'Diagnostic: Airfoil Shape');
plot(xu(1,:), xu(2,:), 'b-'); hold on;
plot(xl(1,:), xl(2,:), 'r-');
axis equal; grid on; title('Optimized Airfoil Shape');
legend('Upper', 'Lower');

% Warning if crossing (common cause of bad drag)
if any(xl(2,:) > xu(2,:))
    warning('CRITICAL: Lower surface crosses Upper surface! This often causes bad drag data.');
end

% --- 4. Prepare Q3D Inputs ---
b_k = 4.36 + 3.95/2;
x_le_r = 0; x_te_r = c_r; 
x_te_k = x_te_r + b_k * tan(deg2rad(0.01)); x_le_k = x_te_k - c_k;
tan_sweep_le = (x_le_k - x_le_r) / b_k;
x_le_t = x_le_r + (b/2) * tan_sweep_le;
z_k = b_k * tan(deg2rad(5)); z_t = (b/2) * tan(deg2rad(5));

AC.Wing.Geom = [x_le_r, 0, 0, c_r, 0;
                x_le_k, b_k, z_k, c_k, -1.0;
                x_le_t, b/2, z_t, c_t, -3];
AC.Wing.inc = 0;   
AC.Wing.Airfoils = [CST; CST; CST];
AC.Wing.eta = [0; b_k/(b/2); 1];

[~, rho_cr, T_cr] = stdatm(h_cr);
AC.Aero.V = M_cr * sqrt(1.4 * 287 * T_cr);
AC.Aero.rho = rho_cr;
AC.Aero.alt = h_cr;
AC.Aero.M = M_cr;
S_wing = 2 * ((c_r+c_k)*b_k/2 + (c_k+c_t)*(b/2-b_k)/2);
AC.Aero.CL = W_total / (0.5 * rho_cr * AC.Aero.V^2 * S_wing);

% --- 5. Run Q3D (VISCOUS MODE) ---
fprintf('\nRunning Q3D (Viscous=1, Iter=500)...\n');
AC.Visc = 1; 
AC.Aero.MaxIterIndex = 500; 

% Mock atmosisa
q3dFolder = fullfile(currentDir, 'Q3D');
mockFile = fullfile(q3dFolder, 'atmosisa.m');
fid = fopen(mockFile, 'w');
if fid ~= -1
    fprintf(fid, 'function [T, a, P, rho] = atmosisa(h)\n [P, rho, T] = stdatm(h);\n gamma=1.4; R=287.058;\n a=sqrt(gamma*R*T);\nend');
    fclose(fid);
end

originalD = pwd;
cd(q3dFolder);
if ispc, [~,~]=system('taskkill /F /IM xfoil.exe 2>NUL'); pause(0.2); end

try
    Res = Q3D_solver(AC);
    
    fprintf('\n--- RESULTS ---\n');
    fprintf('CL_wing:  %.5f\n', Res.CLwing);
    fprintf('CD_wing:  %.5f  (Total)\n', Res.CDwing);
    
    % --- FIX: Calculate Induced Drag Manually ---
    if isfield(Res, 'Wing')
        Y = Res.Wing.Yst;
        chord = Res.Wing.chord;
        % Check variable name variations
        if isfield(Res.Wing, 'cd_i')
            cdi = Res.Wing.cd_i;
        elseif isfield(Res.Wing, 'cdi')
            cdi = Res.Wing.cdi;
        else
            error('Could not find induced drag distribution (cd_i or cdi) in Res.Wing');
        end
        
        % Integrate: CD_ind = (2/S) * Integral( c(y) * cdi(y) dy )
        CD_ind = trapz(Y, chord .* cdi) * 2 / S_wing;
        fprintf('CD_ind:   %.5f  (Calculated from distribution)\n', CD_ind);
        
        % Calculate Profile
        Profile = Res.CDwing - CD_ind;
        fprintf('Profile:  %.5f  (Total - Induced)\n', Profile);
        
    else
        error('Res structure missing "Wing" field.');
    end
    
catch ME
    fprintf('\n[CRASH] Q3D_solver failed: %s\n', ME.message);
end

cd(originalD);
if exist(mockFile, 'file'), delete(mockFile); end

function [xu, xl] = get_airfoil_coords(CST)
    x = linspace(0,1,100);
    N1 = 0.5; N2 = 1.0;
    C_class = (x.^N1) .* ((1-x).^N2);
    Au = CST(1:6); Al = CST(7:12);
    S_u = 0; S_l = 0;
    for i = 1:6
       n_fac = factorial(6-1)/(factorial(i-1)*factorial(6-i));
       S_u = S_u + Au(i) * n_fac * (x.^(i-1)) .* ((1-x).^(6-i));
       S_l = S_l + Al(i) * n_fac * (x.^(i-1)) .* ((1-x).^(6-i));
    end
    yu = C_class .* S_u; yl = C_class .* S_l;
    xu = [x; yu]; xl = [x; yl];
end