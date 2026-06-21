function test_nested_outputfcn()
cost_history = [];

f = @(x) (x(1)-1)^2 + (x(2)-2)^2;
lb = [-10;-10]; ub = [10;10];

options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'MaxIterations', 100, ...
    'OutputFcn', @record_iterate);

try
    [x, ~] = fmincon(f, [0;0], [], [], [], [], lb, ub, [], options);
    fprintf('fmincon with nested OutputFcn OK: x=[%.3f, %.3f], %d iterates\n', ...
        x(1), x(2), numel(cost_history));
catch e
    fprintf('fmincon with nested OutputFcn FAILED: %s\n', e.message);
end

    function stop = record_iterate(~, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            cost_history(end + 1, 1) = optimValues.fval;
        end
    end
end
