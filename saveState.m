function stop = saveState(x_norm, optimValues, state, denormalize, constraints_fun)
%saveState An output function to save the history of an optimization.
%
%   STOP = saveState(X,OPTIMVALUES,STATE,DENORMALIZE,CONSTRAINTS_FUN) saves 
%   the current state of the optimizer to be post-processed.

% Use a persistent variable to store the history
persistent history

stop = false;

switch state
    case 'init'
        % Initialize the history structure
        history.x = [];         % Stores normalized design vectors
        history.fval = [];      % Stores objective function values
        history.c = [];         % Stores inequality constraint values
        history.iteration = []; % Stores iteration numbers
        
    case 'iter'
        % At each iteration, append the current data to the history
        if optimValues.iteration > 0 % Don't save the initial (0th) iteration if it's just a test eval
            % Calculate constraint values for the current point
            x_phys = denormalize(x_norm);
            [c, ~] = constraints_fun(x_phys);
            
            % Append all data
            history.x = [history.x, x_norm];
            history.fval = [history.fval, optimValues.fval];
            history.c = [history.c, c];
            history.iteration = [history.iteration, optimValues.iteration];
        end
        
    case 'done'
        % When optimization is finished, save the history to a .mat file
        % and also assign it to the base workspace for immediate access.
        if ~isempty(history) && ~isempty(history.iteration)
            fprintf('Saving optimization history to optimization_history.mat...\n');
            save('optimization_history.mat', 'history');
            assignin('base', 'optimization_history', history);
            fprintf('History saved. Access it in the workspace via the "optimization_history" variable.\n');
        else
            fprintf('No iteration history was recorded.\n');
        end
end

end
