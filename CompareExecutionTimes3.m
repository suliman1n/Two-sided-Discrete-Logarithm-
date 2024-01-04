function CompareExecutionTimes3(N, num_trials)
    % N: Maximum dimension to test
    % num_trials: Number of trials for each dimension

    % Define the range of dimensions to test
    n_values = 2:N;

    % Initialize matrices to store average execution times
    avg_attack_csr_times2 = zeros(length(n_values), 1);
    avg_attack_csr_times3 = zeros(length(n_values), 1);

    for i = 1:length(n_values)
        n = n_values(i);
        attack_csr_times2 = zeros(num_trials, 1);
        attack_csr_times3 = zeros(num_trials, 1);

        for trial = 1:num_trials
            % Generate random pair matrices
            X = GenerateRandomPairMatrix(n, -100, 100);
            Y = GenerateRandomPairMatrix(n, -100, 100);

            Xt = ReduceMatrix(X);
            Yt = ReduceMatrix(Y);
            [CritX,lZ] = CriticalCycle(Xt);
            [CritY,lW] = CriticalCycle(Yt);

        


            k = randi([(2*n-1)*lZ,2*((2*n-1)*lZ)]);
            l = randi([(2*n-1)*lW,2*((2*n-1)*lW)]);
            r = randi([(2*n-1)*lZ,2*((2*n-1)*lZ)]);
            s = randi([(2*n-1)*lW,2*((2*n-1)*lW)]);
            

            

            c = randi([-100, 100]);
            d = randi([-100, 100]);
            p = randi([-100, 100]);
            q = randi([-100, 100]);
            
            % Measure the time for GenerateKey()
            [key, A_p, B_q] = GenerateKey(X, Y, k, l, r, s, p, c, q, d);

            % Measure the time for Attack()
            

            % Measure the time for AttackCSR()
            t_start = tic;
            key_attack_csr2 = AttackCSR(A_p, B_q, X, Y);
            attack_csr_times2(trial) = toc(t_start);

            t_start = tic;
            key_attack_csr3 = Attack3(A_p, B_q, X, Y);
            attack_csr_times3(trial) = toc(t_start);

            % Check if AttackCSR succeeded
            if isequal(key, key_attack_csr2) && isequal(key, key_attack_csr3)
                % Continue with the rest of the code as before
            else
                % Exclude this trial from the average calculation
                attack_csr_times2(trial) = NaN;
                attack_csr_times3(trial) = NaN;
          
            end
        end

        % Calculate the average execution times over trials for this dimension
        avg_attack_csr_times2(i) = customMean(attack_csr_times2);
        avg_attack_csr_times3(i) = customMean(attack_csr_times3);
    end

    % Create a plot of n vs. average execution times
    plot(n_values, avg_attack_csr_times2, 's-', 'DisplayName', 'Attack with 2n-sized matrices');
    hold on;
    plot(n_values, avg_attack_csr_times3, '^-', 'DisplayName', 'Attack with n-sized matrices');
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