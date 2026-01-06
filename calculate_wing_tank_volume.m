
function [Vol_Total, Vol_Inner, Vol_Outer] = calculate_wing_tank_volume( ...
    c_r, c_k, c_t, ...          % Chords: Root, Kink, Tip
    b_k, b, ...        % Spans: Span of kink (y distance), Total Span (b)
    filename, ...               % Airfoil data file
    limit_start, limit_end)    % Tank limits (e.g., 0.15, 0.85)

% CALCULATE_WING_TANK_VOLUME Calculates exact volume using integrated chords

    % 1. Calculate Average Height Ratio (Au) from airfoil data
    % Read the Data
    try
        data = load(filename);
    catch
        error('Could not load file. Ensure filename is correct and has no text headers.');
    end
    
    x = data(:,1);
    
    % Split into Upper and Lower Surfaces
    % Find the Leading Edge (minimum x value)
    [~, le_idx] = min(x);
    
    % Separate the arrays (Top: start to LE, Bottom: LE to end)
    upper_data = data(1:le_idx, :);
    lower_data = data(le_idx:end, :);
    
    % Sort so X is strictly increasing for interpolation
    upper_data = sortrows(upper_data, 1);
    lower_data = sortrows(lower_data, 1);
    
    % Remove duplicate points to ensure unique x values for interpolation
    upper_data = unique(upper_data, 'rows', 'stable');
    lower_data = unique(lower_data, 'rows', 'stable');
    
    % Define the query points
    % Create 1000 points between the start and end limits
    query_x = linspace(limit_start, limit_end, 1000);
    
    % Interpolate both surfaces to the query points
    upper_y = interp1(upper_data(:,1), upper_data(:,2), query_x, 'linear');
    lower_y = interp1(lower_data(:,1), lower_data(:,2), query_x, 'linear');
    
    % Calculate Thickness (Height)
    thickness = upper_y - lower_y;
    
    % Calculate Area and Average Height
    % trapz(x, y) integrates y with respect to x
    Area = trapz(query_x, thickness);
    Width = limit_end - limit_start;
    
    Au = Area / Width;

    % 2. Calculate the non-dimensional tank length
    delta_xc = limit_end - limit_start; 
    
    % 3. Define Section Dimensions
    % Inner Section (Root to Kink)
    l_inner = b_k; 
    
    % Outer Section (Kink to Tip)
    % assuming b is tip-to-tip span
    l_outer = 0.85*(b / 2) - b_k;
    
    % 4. Calculate Volume of Inner Section (Root -> Kink)
    % Formula: V = K * L/3 * (Ca^2 + Ca*Cb + Cb^2)
    Vol_Inner = (Au * delta_xc * l_inner / 3) * (c_r^2 + c_r*c_k + c_k^2);
    
    % 5. Calculate Volume of Outer Section (Kink -> Tip)
    Vol_Outer = (Au * delta_xc * l_outer / 3) * (c_k^2 + c_k*c_t + c_t^2);
    
    % Factor to account for structural components (ribs, spars, etc.)
    n_struct = 0.93;

    % 6. Apply structural factor to both sections
    Vol_Inner = n_struct * Vol_Inner;
    Vol_Outer = n_struct * Vol_Outer;

    % 7. Total Volume (One Wing)
    Vol_Total = Vol_Inner + Vol_Outer;

end

