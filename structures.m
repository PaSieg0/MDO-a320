function [W_wing] = structures(Y, L, M, ...
                               sweep_te_k, b_k, dihedral, ...
                               b, c_r, c_k, c_t, ...
                               MTOW, ZFW, n_max, S_ref, ...
                               mat_props, spar_locs, ...
                               tank_limits, engine_data, airfoils, ...
                               caseName, verbose)
    % STRUCTURES Calculates wing weight using EMWET.
    %
    % Inputs:
    %   Y, L, M     : Arrays for Span (m), Lift (N/m), Moment (Nm/m)
    %   Geom Vars   : sweep_te_k, b_k, dihedral, b, c_r, c_k, c_t
    %   Weights     : MTOW, ZFW, n_max, S_ref
    %   Props       : Structs/Vectors for materials, spars, tanks, engines
    %   caseName    : (Optional) String name for files. Default: 'a320'
    %   verbose     : (Optional) Boolean. Set true to print output. Default: false
    %
    % Output:
    %   W_wing      : Computed wing weight (kg)

    % Handle optional arguments
    if nargin < 20 || isempty(caseName)
        caseName = 'a320';
    end
    if nargin < 21 || isempty(verbose)
        verbose = false;
    end
    
    % Handle directory navigation for EMWET
    originalDir = pwd;
    structuresPath = fileparts(mfilename('fullpath'));
    emwetPath = fullfile(structuresPath, 'EMWET 1.5');
    if exist(emwetPath, 'dir')
        cd(emwetPath);
    end

    %% 1. GEOMETRY CALCULATION
    % ---------------------------------------------------------
    % Replicate logic from loads.m/aero.m to define X, Y, Z for EMWET
    % EMWET defines geometry by the Leading Edge (LE) coordinates.
    
    % Root Station (y=0)
    x_le_r = 0;
    x_te_r = c_r; 
    
    % Kink Station (y=b_k)
    % Calculate x_le_k based on Trailing Edge Sweep
    x_te_k = x_te_r + b_k * tan(deg2rad(sweep_te_k));
    x_le_k = x_te_k - c_k;
    
    % Tip Station (y=b/2)
    % LE sweep is assumed constant from root to tip (loads.m logic)
    tan_sweep_le = (x_le_k - x_le_r) / b_k;
    x_le_t = x_le_r + (b/2) * tan_sweep_le;
    
    % Z-Coordinates based on dihedral
    z_r = 0;
    z_k = b_k * tan(deg2rad(dihedral));
    z_t = (b/2) * tan(deg2rad(dihedral));
    
    % Build Geometry Matrix: [Chord, X_le, Y, Z_le, FrontSpar, RearSpar]
    % Note: Spar locations are passed in spar_locs = [front_pct, rear_pct]
    geo_matrix = [
        c_r, x_le_r, 0,   z_r, spar_locs(1), spar_locs(2);
        c_k, x_le_k, b_k, z_k, spar_locs(1), spar_locs(2);
        c_t, x_le_t, b/2, z_t, spar_locs(1), spar_locs(2)
    ];

    %% 2. GENERATE .INIT FILE
    % ---------------------------------------------------------
    fid = fopen([caseName '.init'], 'w');
    
    % 1. Weights
    fprintf(fid, '%.2f %.2f\n', MTOW, ZFW);
    % 2. Load Factor
    fprintf(fid, '%.2f\n', n_max);
    % 3. Area, Span, Planform Sections, Airfoil Sections
    fprintf(fid, '%.2f %.2f 3 3\n', S_ref, b);
    
    % 4. Airfoil Definitions (Normalized span eta, Filename)
    fprintf(fid, '0 %s\n', airfoils.root);
    fprintf(fid, '%.4f %s\n', (b_k/(b/2)), airfoils.kink);
    fprintf(fid, '1 %s\n', airfoils.tip);
    
    % 5. Geometry Sections
    for i = 1:3
        fprintf(fid, '%.4f %.4f %.4f %.4f %.2f %.2f\n', ...
            geo_matrix(i,1), geo_matrix(i,2), geo_matrix(i,3), ...
            geo_matrix(i,4), geo_matrix(i,5), geo_matrix(i,6));
    end
    
    % 6. Fuel Tank (start eta, end eta)
    fprintf(fid, '%.2f %.2f\n', tank_limits(1), tank_limits(2));
    
    % 7. Engines
    fprintf(fid, '%d\n', engine_data.count);
    eta_eng = engine_data.y_location / (b/2);
    fprintf(fid, '%.4f %.2f\n', eta_eng, engine_data.weight);
    
    % 8. Materials (Upper, Lower, Front, Rear)
    % Format: E, density, yield_tens, yield_comp
    fprintf(fid, '%.4e %.2f %.4e %.4e\n', mat_props.upper);
    fprintf(fid, '%.4e %.2f %.4e %.4e\n', mat_props.lower);
    fprintf(fid, '%.4e %.2f %.4e %.4e\n', mat_props.front);
    fprintf(fid, '%.4e %.2f %.4e %.4e\n', mat_props.rear);
    
    % 9. Panel Efficiency and Rib Pitch (Standard defaults)
    fprintf(fid, '0.96 0.5\n');
    
    % 10. Display Option
    % If verbose is true, EMWET displays results (1), else (0)
    disp_opt = 0;
    if verbose
        disp_opt = 1; 
    end
    fprintf(fid, '%d\n', disp_opt); 
    
    fclose(fid);
    
    %% 3. GENERATE .LOAD FILE
    % ---------------------------------------------------------
    % Convert inputs Y (m), L (N/m), M (Nm/m) to EMWET format
    % EMWET expects: eta, Force (N), Moment (Nm)
    
    eta_vec = Y ./ (b/2);
    
    fid = fopen([caseName '.load'], 'w');
    for i = 1:length(eta_vec)
        fprintf(fid, '%.4f %.4e %.4e\n', eta_vec(i), L(i), M(i));
    end
    fclose(fid);
    
    %% 4. EXECUTE EMWET
    % ---------------------------------------------------------
    if verbose
        fprintf('Running EMWET for: %s...\n', caseName);
    end
    
    EMWET(caseName); 
    
    %% 5. READ OUTPUT
    % ---------------------------------------------------------
    outputFile = [caseName '.weight'];
    
    if exist(outputFile, 'file')
        fid = fopen(outputFile, 'r');
        line = fgetl(fid); 
        fclose(fid);
        
        % Extract number: "Wing total weight(kg) 4436.79"
        W_wing = sscanf(line, 'Wing total weight(kg) %f');
        
        if verbose
            fprintf('--------------------------------------\n');
            fprintf('Wing Weight (W_wing): %.2f kg\n', W_wing);
            fprintf('--------------------------------------\n');
        end
    else
        error('Output file was not created.');
    end
    
    % Return to original directory
    cd(originalDir);
end