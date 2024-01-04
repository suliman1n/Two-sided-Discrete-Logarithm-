function reducedMatrix = ReduceMatrix(A)
    [m, n] = size(A);
    reducedMatrix = zeros(m, n);

    for i = 1:m
        for j = 1:n
            tuple = A{i, j};
            max_value = max(tuple);
            reducedMatrix(i, j) = max_value;
        end
    end
end