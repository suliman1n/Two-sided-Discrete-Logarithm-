function [resultMatrix, A_p, B_q] = GenerateKey(X, Y, k, l, r, s, p, c, q, d)
    % Compute PairMatrixPower for X with exponents k and r
    Xk = fastpowermaxpluspairt(X, k);
    Xr = fastpowermaxpluspairt(X, r);

    % Compute PairMatrixPower for Y with exponents s and l
    Ys = fastpowermaxpluspairt(Y, s);
    Yl = fastpowermaxpluspairt(Y, l);

    % Use PairMatrixMultiply to combine the results
    tempMatrix1 = PairMatrixMultiply(Xk, Xr);
    tempMatrix2 = PairMatrixMultiply(Ys, Yl);
    tempMatrix = PairMatrixMultiply(tempMatrix1, tempMatrix2);

    % Compute the sum q + p + (r-1)d + (k-1)c + (l-1)c + (s-1)d
    scalarValue = q + p + (r-1)*d + (k-1)*c + (l-1)*c + (s-1)*d;

    A_p = TropicalMultiScalarMatrix(PairMatrixMultiply(fastpowermaxpluspairt(X, k), fastpowermaxpluspairt(Y, l)), p + (k + l - 2) * c);
    B_q = TropicalMultiScalarMatrix(PairMatrixMultiply(fastpowermaxpluspairt(X, r), fastpowermaxpluspairt(Y, s)), q + (r + s - 2) * d);

    % Use TropicalMultiScalarMatrix to add the scalar to the result
    resultMatrix = TropicalMultiScalarMatrix(tempMatrix, scalarValue);
end