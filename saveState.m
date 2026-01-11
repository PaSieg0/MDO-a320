function stop = saveState(x_norm, optimValues, state, denormalize, constraints_fun)
%saveState An output function to save the history of an optimization to both
%a .mat file and a .csv file for portability.

% Use a persistent variable to store the history
persistent history

% Defensive check: If 'history' is not a struct (e.g., it was cleared),
% re-initialize it to prevent dot indexing errors.
if ~isstruct(history)
    history.x = [];
    history.fval = [];
    history.c = [];
    history.iteration = [];
end

stop = false;

switch state
    case 'init'
        fprintf('\n>> saveState called with state: init\n');
        % Initialize/reset the history structure
        history.x = [];         % Stores normalized design vectors
        history.fval = [];      % Stores objective function values
        history.c = [];         % Stores inequality constraint values
        history.iteration = []; % Stores iteration numbers
        
        % --- Setup for CSV logging ---
        csv_filename = 'optimization_history.csv';
        try
            fid = fopen(csv_filename, 'w');
            if fid == -1
                error('Cannot open CSV file for writing: %s', csv_filename);
            end
            % Construct and write the header row
            num_x_vars = 18; % There are 18 design variables
            x_headers = cell(1, num_x_vars);
            for i = 1:num_x_vars
                x_headers{i} = ['x' num2str(i)];
            end
            header_str = strjoin(['Iteration', 'Objective', 'C1', 'C2', 'C3', x_headers], ',');
            fprintf(fid, '%s\n', header_str);
            fclose(fid);
        catch ME
            warning('Could not create CSV log file: %s', ME.message);
        end

    case 'iter'
        fprintf('\n>> saveState called with state: iter (Iteration %d)\n', optimValues.iteration);
        % At each iteration, append the current data to the history.
        % This now includes iteration 0.

        % Calculate constraint values for the current point
        x_phys = denormalize(x_norm);
        [c, ~] = constraints_fun(x_phys);
        
        % --- Display the design vector for the current iteration ---
        fprintf('--- DEBUG: Normalized Vector (x_norm) ---\n');
        for i = 1:length(x_norm)
            fprintf('x_norm(%2d): %+.6f\n', i, x_norm(i));
        end
        
        fprintf('\n--- DEBUG: Physical Vector (x_phys) ---\n');
        for i = 1:length(x_phys)
            fprintf('x_phys(%2d): %g\n', i, x_phys(i));
        end
        fprintf('--- End of Vector Display ---\n');
        
        % Append all data to the persistent history struct
        history.x = [history.x, x_norm];
        history.fval = [history.fval, optimValues.fval];
        history.c = [history.c, c];
        history.iteration = [history.iteration, optimValues.iteration];

        % --- Append latest iteration to CSV file ---
        csv_filename = 'optimization_history.csv';
        try
            fid = fopen(csv_filename, 'a');
            if fid == -1
                error('Cannot open CSV file for appending: %s', csv_filename);
            end
            % Prepare data row, ensuring all parts are row vectors before concatenation
            data_row = [optimValues.iteration, optimValues.fval, c(:)', x_norm(:)'];
            % Write data using a generic format specifier
            fprintf(fid, ['%d,' repmat('%.8g,', 1, length(data_row)-1) '%.8g\n'], data_row);
            fclose(fid);
        catch ME
            warning('Could not write to CSV log file: %s', ME.message);
        end

        % --- Save .mat file on every iteration for robustness ---
        save('optimization_history.mat', 'history');
        
    case 'done'
        fprintf('\n>> saveState called with state: done\n');
        % When optimization is finished, assign the final history to the base 
        % workspace for immediate access. The .mat and .csv files are already up-to-date.
        if ~isempty(history) && ~isempty(history.iteration)
            fprintf('Optimization finished. Final history saved in "optimization_history.mat" and "optimization_history.csv".\n');
            assignin('base', 'optimization_history', history);
            fprintf('History also assigned to workspace variable "optimization_history".\n');
        else
            fprintf('No iteration history was recorded.\n');
        end
end

end
