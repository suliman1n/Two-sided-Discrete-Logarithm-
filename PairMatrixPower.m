function resultMatrix = PairMatrixPower(A, p)
    if p < 1
        error('Exponent p should be a positive integer.');
    end

    % Initialize the result matrix as A
    resultMatrix = A;

    % Compute the p-th power using repeated PairMatrixMultiply
    for i = 2:p
        resultMatrix = PairMatrixMultiply(resultMatrix, A);
    end
end