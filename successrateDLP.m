function successRate = successrateDLP(n, mm, mM, max_t, numInstances)
    successCount = 0;

    for instance = 1:numInstances
        alpha = randi([1, 1000]);
        t_1 = randi([1, max_t]);
        t_2 = randi([1, max_t]);

        D_1 = randi([mm, mM], n);
        D_2 = randi([mm, mM], n);

        U = alpha+TropMulti(fastpowermaxplus(D_1, t_1), fastpowermaxplus(D_2, t_2));

        [t_1_, t_2_, beta] = ModTwoSidedTrialandErrorpairs(U, D_1, D_2, n^2);

        check = U - TropMulti(fastpowermaxplus(D_1, t_1_), fastpowermaxplus(D_2, t_2_));
        if all(check == check(1))
            successCount = successCount + 1;
        end
    end

    successRate = successCount / numInstances;
end
