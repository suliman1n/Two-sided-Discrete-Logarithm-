function PlotSuccessRatesVsDimensionsDLP(N,num_trials)
    % N: Maximum dimension to test
    % num_trials: Number of trials for each dimension

    % Define the range of dimensions to test
    n_values = 2:N;

    % Initialize a vector to store success rates
    success_rates = zeros(size(n_values));

    for i = 1:length(n_values)
        % Calculate the success rate for the current dimension
        n = n_values(i);
        success_rate = successrateDLP(n,-100,100,n^3,num_trials);

        % Store the success rate in the vector
        success_rates(i) = success_rate;
    end

    % Create a plot of n vs. success rates
    plot(n_values, success_rates, 'o-');
    xlabel('Dimension (n)');
    ylabel('Success Rate');
    title('Success Rate vs. Dimension');
end