function stop = plotter(x_norm, optimValues, state, denormalize, value_to_plot)
% PLOTTER plots a specific value during fmincon optimization.
% Each value gets its own figure to prevent conflicts.

persistent history_map figure_map
stop = false;

% Use the name of the value as the key for the history and figure
key = value_to_plot;

switch state
    case 'init'
        if isempty(history_map)
            history_map = containers.Map();
            figure_map = containers.Map();
        end
        
        % Initialize history for this specific plot
        history.iter = [];
        history.val = [];
        history_map(key) = history;
        
        % Create and store a handle to a new figure
        fig = figure('Name', ['Convergence: ' get_title(key)]);
        figure_map(key) = fig;
        
        % Setup the plot
        hold on;
        grid on;
        title(get_title(key));
        xlabel('Iteration');
        ylabel(get_ylabel(key));

    case 'iter'
        if optimValues.iteration == 0
            return;
        end
        
        % Retrieve the history and figure for this plot
        history = history_map(key);
        fig = figure_map(key);
        
        % Check if the figure was closed by the user
        if ~ishandle(fig)
            stop = true;
            fprintf('Plotter: Figure for "%s" was closed. Stopping.\n', key);
            return;
        end
        
        % --- Calculate the value to plot ---
        current_val = 0;
        switch key
            case 'fval'
                current_val = optimValues.fval;
            case {'c1', 'c2', 'c3'}
                % Only calculate constraints once per iteration if needed
                x_phys = denormalize(x_norm);
                [c, ~] = constraints(x_phys);
                if strcmp(key, 'c1')
                    current_val = c(1);
                elseif strcmp(key, 'c2')
                    current_val = c(2);
                else % c3
                    current_val = c(3);
                end
        end

        % Append data to history
        history.iter = [history.iter, optimValues.iteration];
        history.val = [history.val, current_val];
        history_map(key) = history;

        % --- Update the plot ---
        figure(fig); % Make sure we are drawing on the correct figure
        plot(history.iter, history.val, '-o', 'MarkerSize', 4);
        
        % For constraints, draw a zero line
        if startsWith(key, 'c')
            line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--');
        end
        
        drawnow;

    case 'done'
        % Clean up maps if you want, or leave them for next run
        % clear history_map figure_map
end

end

function str = get_title(key)
    switch key
        case 'fval', str = 'Objective Function';
        case 'c1',   str = 'Constraint c1: Wing Loading';
        case 'c2',   str = 'Constraint c2: Fuel Capacity';
        case 'c3',   str = 'Constraint c3: Max Fuel';
        otherwise,   str = 'Unknown';
    end
end

function str = get_ylabel(key)
    switch key
        case 'fval', str = 'Objective Value';
        case 'c1',   str = 'c1 Value';
        case 'c2',   str = 'c2 Value';
        case 'c3',   str = 'c3 Value';
        otherwise,   str = 'Value';
    end
end
