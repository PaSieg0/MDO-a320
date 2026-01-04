clear all;
close all;
clc;

% Setup paths
currentDir = fileparts(mfilename('fullpath'));
q3dPath = fullfile(currentDir, 'Q3D');
addpath(q3dPath);
cd(q3dPath);

%% Aerodynamic solver setting
% Wing planform geometry
% x y z chord(m) twist angle (deg)
AC.Wing.Geom = [0   0    0  3.5  0; 
                0.9 14.5 0  1.4  0];
% Wing incidence angle (degree)
AC.Wing.inc = 0;
% Airfoil coefficients input matrix
% | -> upper curve coeff. <-| | -> lower curve coeff. <-|
AC.Wing.Airfoils =[0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797;
0.2171 0.3450 0.2975 0.2685 0.2893 -0.1299 -0.2388 -0.1635 -0.0476 0.0797];
AC.Wing.eta = [0;1]; % Spanwise location of the airfoil sections
% Viscous vs inviscid
AC.Visc = 1; % 0 for inviscid and 1 for viscous analysis
AC.Aero.MaxIterIndex = 150;
% Flight Condition
AC.Aero.V = 230; % flight speed (m/s)
AC.Aero.rho = 1.225; % air density (kg/m3)
AC.Aero.alt = 0; % flight altitude (m)
AC.Aero.Re = 1.14e7; % reynolds number (based on mean aerodynamic chord)
AC.Aero.M = 0.78; % flight Mach number
AC.Aero.CL = 0.4; % lift coefficient - comment this line to run the code for given alpha%
% AC.Aero.Alpha = 2; % angle of attack - comment this line to run the code for given cl
Res = Q3D_solver(AC);

Cl = Res.CLwing; % lift coefficient 
Cd = Res.CDwing; % drag coefficient
fprintf('Lift Coefficient (Cl): %.4f\n', Cl);
fprintf('Drag Coefficient (Cd): %.4f\n', Cd);
fprintf('Lift-to-Drag Ratio (L/D): %.2f\n', Cl / Cd);