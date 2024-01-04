function success_rate = SuccessRate3(n, num_instances)
    count = 0;

    for instance = 1:num_instances
        minValue = -100;
        maxValue = 100;

        X = GenerateRandomPairMatrix(n, minValue, maxValue);
        Y = GenerateRandomPairMatrix(n, minValue, maxValue);
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

        [key, A_p, B_q] = GenerateKey(X, Y, k, l, r, s, p, c, q, d);
        key_attack = Attack3(A_p, B_q, X, Y);

        if isequal(key, key_attack)
            count = count + 1;
        end
    end

    success_rate = count / num_instances;
end