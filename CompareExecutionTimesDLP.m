function CompareExecutionTimesDLP(N, num_trials)
    

    
    n_values = 2:N;

    
    avg_csr_times = zeros(length(n_values), 1);
    avg_trialerror_times = zeros(length(n_values), 1);

    for i = 1:length(n_values)
        n = n_values(i);
        csr_times = zeros(num_trials, 1);
        trialerror_times = zeros(num_trials, 1);

        for trial = 1:num_trials
            
            D_1 = randi([-1000, 1000], n);
            D_2 = randi([-1000, 1000], n);
            [CritX,lX] = CriticalCycle(D_1);
            [CritY,lY] = CriticalCycle(D_2);
            alpha = randi([1, 1000]);
            t_1 = randi([(n-1)*lX+1, 2*((n-1)*lX+1)]);
            t_2 = randi([(n-1)*lY+1, 2*((n-1)*lY+1)]);
            
            U = alpha+TropMulti(fastpowermaxplus(D_1, t_1), fastpowermaxplus(D_2, t_2));

           
             
            
            t_start = tic;
            [t_1_,t_2_,beta_1] = ModTwoSidedCSR(U,D_1,TropId(n),D_2);
            csr_times(trial) = toc(t_start);

            t_start = tic;
            [s_1_,s_2_,beta_2] = ModTwoSidedTrialandErrorpairs(U, D_1, D_2, 2*((n-2)*lY+n+1));
            trialerror_times(trial) = toc(t_start);

            % Check if AttackCSR succeeded
        
        end

        % Calculate the average execution times over trials for this dimension
        avg_csr_times(i) = customMean(csr_times);
        avg_trialerror_times(i) = customMean(trialerror_times);
    end

    % Create a plot of n vs. average execution times
    plot(n_values, avg_trialerror_times, 's-', 'DisplayName', 'Algorithm 1');
    hold on;
    plot(n_values, avg_csr_times, '^-', 'DisplayName', 'Algorithm 2');
    hold off;

    xlabel('Dimension (n)');
    ylabel('Average Execution Time (seconds)');
    title('Average Execution Times vs. Dimension');
    legend('Location', 'Northwest');
end

function avg = customMean(values)
    % Custom mean calculation excluding NaN values
    non_nan_values = values(~isnan(values));
    avg = mean(non_nan_values);
end