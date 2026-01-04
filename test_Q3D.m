clear all;
close all;
clc;

% Setup paths
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
addpath(q3dPath);
cd(q3dPath);

%% Aerodynamic solver setting
% Wing planform geometry with sweep for transonic convergence
% x y z chord(m) twist angle (deg)
% Adding LE sweep (~25 deg) helps reduce effective Mach on airfoil sections
AC.Wing.Geom = [0      0    0    3.5   0; 
                6.76   14.5 0    1.4  -2];  % x = 14.5*tan(25deg) = 6.76m sweep
% Wing incidence angle (degree)
AC.Wing.inc = 0;

% Airfoil coefficients - Using thinner airfoil for transonic
% SC(2)-0610 supercritical-like CST coefficients (t/c ~ 10%)
CST_transonic = [0.1800 0.2800 0.2400 0.2200 0.2400 -0.1000 -0.1800 -0.1200 -0.0300 0.0600];
AC.Wing.Airfoils = [CST_transonic;
                    CST_transonic];
AC.Wing.eta = [0;1]; % Spanwise location of the airfoil sections

% Viscous vs inviscid
AC.Visc = 1; % 0 for inviscid and 1 for viscous analysis
AC.Aero.MaxIterIndex = 150;

% Flight Condition
AC.Aero.V = 230; % flight speed (m/s) - UNCHANGED
AC.Aero.rho = 1.225; % air density (kg/m3)
AC.Aero.alt = 0; % flight altitude (m)
AC.Aero.Re = 1.14e7; % reynolds number (based on mean aerodynamic chord)
AC.Aero.M = 0.78; % flight Mach number - UNCHANGED
AC.Aero.CL = 0.35; % Reduced CL target for better transonic convergence
% AC.Aero.Alpha = 2; % angle of attack - comment this line to run the code for given cl
Res = Q3D_solver(AC);

Cl = Res.CLwing; % lift coefficient 
Cd = Res.CDwing; % drag coefficient
fprintf('Lift Coefficient (Cl): %.4f\n', Cl);
fprintf('Drag Coefficient (Cd): %.4f\n', Cd);
fprintf('Lift-to-Drag Ratio (L/D): %.2f\n', Cl / Cd);