function resultMatrix = PairMatrixMultiply(X, Y)
    [mX, nX] = size(X);
    [mY, nY] = size(Y);

    if nX ~= mY
        error('Matrix dimensions are not compatible for multiplication.');
    end

    resultMatrix = cell(mX, nY);

    for i = 1:mX
        for j = 1:nY
            tuple = [-Inf, -Inf];

            for k = 1:nX
                product = TupleMultiply(X{i, k}, Y{k, j});
                tuple = TupleSum(tuple, product);
            end

            resultMatrix{i, j} = tuple;
        end
    end
end