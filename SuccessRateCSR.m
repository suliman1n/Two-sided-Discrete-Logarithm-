function success_rate = SuccessRateCSR(n, num_instances)
    count = 0;

    for instance = 1:num_instances
        minValue = -100;
        maxValue = 100;

        X = GenerateRandomPairMatrix(n, minValue, maxValue);
        Y = GenerateRandomPairMatrix(n, minValue, maxValue);
        Xt = TransformToTypicalMatrix(X);
        Yt = TransformToTypicalMatrix(Y);
        [CritX,lZ] = CriticalCycle(Xt);
        [CritY,lW] = CriticalCycle(Yt);

        %k = 3;
        %l = 3;
        %r = 3;
        %s = 3;

        %k = randi([1, (n-1)^2+1]);
        %l = randi([1, (n-1)^2+1]);
        %r = randi([1, (n-1)^2+1]);
        %s = randi([1, (n-1)^2+1]);

        %k = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);
        %l = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);
        %r = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);
        %s = randi([(n-1)^2+1+1, 2*((n-1)^2+1)]);

        %k = randi([(2*n-1)^2+1+1, 2*((2*n-1)^2+1)]);
        %l = randi([(2*n-1)^2+1+1, 2*((2*n-1)^2+1)]);
        %r = randi([(2*n-1)^2+1+1, 2*((2*n-1)^2+1)]);
        %s = randi([(2*n-1)^2+1+1, 2*((2*n-1)^2+1)]);
        


        k = randi([1, (n-1)*lZ]);
        l = randi([1, (n-1)*lW]);
        r = randi([1, (n-1)*lZ]);
        s = randi([1, (n-1)*lW]);

        %k = randi([(n-1)*lZ, (2*n-1)*lZ]);
        %l = randi([(n-1)*lW, (2*n-1)*lW]);
        %r = randi([(n-1)*lZ, (2*n-1)*lZ]);
        %s = randi([(n-1)*lW, (2*n-1)*lW]);

        %k = randi([(2*n-1)*lZ,2*((2*n-1)*lZ)]);
        %l = randi([(2*n-1)*lW,2*((2*n-1)*lW)]);
        %r = randi([(2*n-1)*lZ,2*((2*n-1)*lZ)]);
        %s = randi([(2*n-1)*lW,2*((2*n-1)*lW)]);

        c = randi([-100, 100]);
        d = randi([-100, 100]);
        p = randi([-100, 100]);
        q = randi([-100, 100]);

        [key, A_p, B_q] = GenerateKey(X, Y, k, l, r, s, p, c, q, d);
        key_attack = AttackCSR(A_p, B_q, X, Y);

        if isequal(key, key_attack)
            count = count + 1;
        end
    end

    success_rate = count / num_instances;
end