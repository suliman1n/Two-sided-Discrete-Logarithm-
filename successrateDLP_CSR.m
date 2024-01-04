function successRate = successrateDLP_CSR(n, mm, mM,numInstances)
    successCount = 0;

    for instance = 1:numInstances

        D_1 = randi([mm, mM], n);
        D_2 = randi([mm, mM], n);
        [CritX,lX] = CriticalCycle(D_1);
        [CritY,lY] = CriticalCycle(D_2);

        alpha = randi([1, 1000]);
        t_1 = randi([1, (n-1)*lX]);
        t_2 = randi([1, (n-1)*lY]);

        %t_1 = randi([(n-1)*lX+1, n^3]);
        %t_2 = randi([(n-1)*lY+1, n^3]);

        U = alpha+TropMulti(fastpowermaxplus(D_1, t_1), fastpowermaxplus(D_2, t_2));

        [t_1_, t_2_, beta] = ModTwoSidedCSR(U, D_1,TropId(n), D_2);

        check = U - TropMulti(fastpowermaxplus(D_1, t_1_), fastpowermaxplus(D_2, t_2_));
        if all(check == check(1))
            successCount = successCount + 1;
        end
    end

    successRate = successCount / numInstances;
end
