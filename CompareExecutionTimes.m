function CompareExecutionTimes(N, num_trials)
    % N: Maximum dimension to test
    % num_trials: Number of trials for each dimension

    % Define the range of dimensions to test
    n_values = 2:N;

    % Initialize matrices to store average execution times
    avg_generate_key_times = zeros(length(n_values), 1);
    avg_attack_times = zeros(length(n_values), 1);
    avg_attack_csr_times = zeros(length(n_values), 1);

    for i = 1:length(n_values)
        n = n_values(i);
        generate_key_times = zeros(num_trials, 1);
        attack_times = zeros(num_trials, 1);
        attack_csr_times = zeros(num_trials, 1);

        for trial = 1:num_trials
            % Generate random pair matrices
            t_start = tic;

            %k = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);
            %l = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);
            %r = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);
            %s = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);

            

            c = randi([-100, 100]);
            d = randi([-100, 100]);
            p = randi([-100, 100]);
            q = randi([-100, 100]);
            X = GenerateRandomPairMatrix(n, -100, 100);
            Y = GenerateRandomPairMatrix(n, -100, 100);
            Xt = TransformToTypicalMatrix(X);
            Yt = TransformToTypicalMatrix(Y);
            [CritX,lZ] = CriticalCycle(Xt);
            [CritY,lW] = CriticalCycle(Yt);
            k = randi([(2*n-1)*lZ,2*((2*n-1)*lZ)]);
            l = randi([(2*n-1)*lW,2*((2*n-1)*lW)]);
            r = randi([(2*n-1)*lZ,2*((2*n-1)*lZ)]);
            s = randi([(2*n-1)*lW,2*((2*n-1)*lW)]);

            % Measure the time for GenerateKey()
            [key, A_p, B_q] = GenerateKey(X, Y, k, l, r, s, p, c, q, d);
            generate_key_times(trial) = toc(t_start);

            % Measure the time for Attack()
            t_start = tic;
            key_attack = Attack(A_p, B_q, X, Y, 500);
            attack_times(trial) = toc(t_start);

            % Measure the time for AttackCSR()
            t_start = tic;
            key_attack_csr = AttackCSR(A_p, B_q, X, Y);
            attack_csr_times(trial) = toc(t_start);

            % Check if AttackCSR succeeded
            if isequal(key, key_attack_csr)
                % Continue with the rest of the code as before
            else
                % Exclude this trial from the average calculation
                attack_csr_times(trial) = NaN;
                attack_times(trial)=NaN;
                generate_key_times(trial)=NaN;
            end
        end

        % Calculate the average execution times over trials for this dimension
        avg_generate_key_times(i) = customMean(generate_key_times);
        avg_attack_times(i) = customMean(attack_times);
        avg_attack_csr_times(i) = customMean(attack_csr_times);
    end

    % Create a plot of n vs. average execution times
    plot(n_values, avg_generate_key_times, 'o-', 'DisplayName', 'Key generation');
    hold on;
    plot(n_values, avg_attack_times, 's-', 'DisplayName', 'Attack using the non-unique solution');
    plot(n_values, avg_attack_csr_times, '^-', 'DisplayName', 'Attack using CSR');
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