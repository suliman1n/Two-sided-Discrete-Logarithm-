function success_rate = SuccessRate(n, max_t, num_instances)
    count = 0;

    for instance = 1:num_instances
        minValue = -100;
        maxValue = 100;

        X = GenerateRandomPairMatrix(n, minValue, maxValue);
        Y = GenerateRandomPairMatrix(n, minValue, maxValue);

        %k = randi([1, 100]);
        %l = randi([1, 100]);
        %r = randi([1, 100]);
        %s = randi([1, 100]);

        %k = 2;
        %l = 3;
        %r = 2;
        %s = 3;

        k = randi([1, n^5]);
        l = randi([1, n^5]);
        r = randi([1, n^5]);
        s = randi([1, n^5]);

        %k = randi([(n-1)^2+1, 2*((n-1)^2+1)]);
        %l = randi([(n-1)^2+1, 2*((n-1)^2+1)]);
        %r = randi([(n-1)^2+1, 2*((n-1)^2+1)]);
        %s = randi([(n-1)^2+1, 2*((n-1)^2+1)]);

        %k = randi([(2*n-1)^2+1, 2*((2*n-1)^2+1)]);
        %l = randi([(2*n-1)^2+1, 2*((2*n-1)^2+1)]);
        %r = randi([(2*n-1)^2+1, 2*((2*n-1)^2+1)]);
        %s = randi([(2*n-1)^2+1, 2*((2*n-1)^2+1)]);

        c = randi([-100, 100]);
        d = randi([-100, 100]);
        p = randi([-100, 100]);
        q = randi([-100, 100]);

        [key, A_p, B_q] = GenerateKey(X, Y, k, l, r, s, p, c, q, d);
        key_attack = Attack(A_p, B_q, X, Y, max_t);

        if isequal(key, key_attack)
            count = count + 1;
        end
    end

    success_rate = count / num_instances;
end