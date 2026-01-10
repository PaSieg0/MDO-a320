function stop = plotConvergence(x_norm, optimValues, state, denormalize)
% PLOTCONVERGENCE plots the objective function and all constraints during optimization.
%
%   STOP = PLOTCONVERGENCE(X,OPTIMVALUES,STATE,DENORMALIZE) where X is the
%   current point, OPTIMVALUES is a structure containing data from the
%   current iteration, STATE is the current state of the algorithm, and
%   DENORMALIZE is a function handle to denormalize the design variables.

persistent history

stop = false;

switch state
    case 'init'
        % Create a figure, arrange subplots
        figure('Name', 'Optimization Convergence');
        
        % Subplot for objective function
        subplot(2, 2, 1);
        title('Objective Function');
        xlabel('Iteration');
        ylabel('Objective Value');
        grid on;
        hold on;

        % Subplot for Constraint 1
        subplot(2, 2, 2);
        title('Constraint c1: Wing Loading');
        xlabel('Iteration');
        ylabel('c1 Value');
        grid on;
        hold on;

        % Subplot for Constraint 2
        subplot(2, 2, 3);
        title('Constraint c2: Fuel Capacity');
        xlabel('Iteration');
        ylabel('c2 Value');
        grid on;
        hold on;

        % Subplot for Constraint 3
        subplot(2, 2, 4);
        title('Constraint c3: Max Fuel');
        xlabel('Iteration');
        ylabel('c3 Value');
        grid on;
        hold on;
        
        % Initialize history
        history.fval = [];
        history.c = [];
        history.iter = [];

    case 'iter'
        % Don't plot anything if it's the first evaluation (iter 0)
        if optimValues.iteration == 0
            return;
        end

        % Denormalize x to get physical values
        x_phys = denormalize(x_norm);
        
        % Calculate constraints at the current point
        [c_current, ~] = constraints(x_phys);
        
        % Append data to history
        history.fval = [history.fval, optimValues.fval];
        history.c = [history.c, c_current]; % Appends a column
        history.iter = [history.iter, optimValues.iteration];
        
        % --- Update Plots ---
        
        % Plot objective function
        subplot(2, 2, 1);
        plot(history.iter, history.fval, 'b-o', 'MarkerSize', 4);
        
        % Plot constraints
        subplot(2, 2, 2);
        plot(history.iter, history.c(1,:), 'r-o', 'MarkerSize', 4);
        line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--'); % Zero line

        subplot(2, 2, 3);
        plot(history.iter, history.c(2,:), 'g-o', 'MarkerSize', 4);
        line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--'); % Zero line

        subplot(2, 2, 4);
        plot(history.iter, history.c(3,:), 'm-o', 'MarkerSize', 4);
        line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--'); % Zero line
        
        % Force plot to update
        drawnow;

    case 'done'
        % Release figure hold on all subplots
        for i = 1:4
            subplot(2, 2, i);
            hold off;
        end
end
end
