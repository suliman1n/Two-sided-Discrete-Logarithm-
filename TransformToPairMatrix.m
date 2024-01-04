function pairMatrix = TransformToPairMatrix(matrix)
    [m, n] = size(matrix);
    
    if mod(m, 2) ~= 0 || mod(n, 2) ~= 0
        error('Input matrix dimensions are not suitable for conversion to a pair matrix.');
    end
    
    pairMatrix = cell(m/2, n/2);
    
    for i = 1:m/2
        for j = 1:n/2
            tuple = [matrix(2*i-1, 2*j-1), matrix(2*i-1, 2*j)];
            pairMatrix{i, j} = tuple;
        end
    end
end