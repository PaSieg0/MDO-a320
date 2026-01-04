function cst2dat(A, filename, n_points)
% CST2DAT Converts CST parameters to a .dat airfoil file
%   cst2dat(A, filename) generates an airfoil .dat file from CST parameters
%   cst2dat(A, filename, n_points) specifies the number of points per surface
%
%   Inputs:
%       A        - CST coefficients [Au, Al] concatenated (Au first, then Al)
%       filename - Output filename (without .dat extension)
%       n_points - Number of points per surface (default: 50)
%
%   The CST (Class Shape Transformation) method uses Bernstein polynomials
%   to define the airfoil shape.

    if nargin < 3
        n_points = 50;
    end
    
    % Split A into upper and lower surface coefficients
    n_coeff = length(A) / 2;
    Au = A(1:n_coeff);
    Al = A(n_coeff+1:end);
    
    % Generate x coordinates (cosine spacing for better leading edge resolution)
    beta = linspace(0, pi, n_points);
    x = (1 - cos(beta)) / 2;
    
    % Calculate upper and lower surface coordinates
    y_upper = cst_curve(x, Au);
    y_lower = cst_curve(x, Al);
    
    % Combine coordinates (start from trailing edge, go around leading edge)
    x_coords = [flip(x), x(2:end)]';
    y_coords = [flip(y_upper), y_lower(2:end)]';
    
    % Write to .dat file
    filepath = fullfile(fileparts(mfilename('fullpath')), [filename, '.dat']);
    fid = fopen(filepath, 'w');
    
    if fid == -1
        error('Could not open file for writing: %s', filepath);
    end
    
    % Write coordinates in scientific notation format (similar to b737a.dat)
    for i = 1:length(x_coords)
        fprintf(fid, '  %14.7e  %14.7e\n', x_coords(i), y_coords(i));
    end
    
    fclose(fid);
    fprintf('Airfoil data written to: %s\n', filepath);
end

function y = cst_curve(x, A)
% CST_CURVE Calculates y coordinates using CST method
%   y = cst_curve(x, A) computes the y-coordinates for given x and CST coefficients
%
%   Inputs:
%       x - x-coordinates (0 to 1)
%       A - CST coefficients (Bernstein polynomial weights)
%
%   Output:
%       y - y-coordinates

    % Class function for airfoil: C(x) = x^N1 * (1-x)^N2
    % Standard values for airfoil: N1 = 0.5, N2 = 1.0
    N1 = 0.5;
    N2 = 1.0;
    
    % Class function
    C = (x.^N1) .* ((1 - x).^N2);
    
    % Shape function using Bernstein polynomials
    n = length(A) - 1;  % Order of Bernstein polynomial
    S = zeros(size(x));
    
    for i = 0:n
        % Bernstein polynomial basis
        K = nchoosek(n, i);
        B = K .* (x.^i) .* ((1 - x).^(n - i));
        S = S + A(i + 1) * B;
    end
    
    % CST representation: y = C(x) * S(x)
    y = C .* S;
end