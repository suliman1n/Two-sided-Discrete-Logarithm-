function typicalMatrix = TransformToTypicalMatrix(pairMatrix)
    % Get the size of the input pair matrix
    [m, n] = size(pairMatrix);

    % Initialize the typical matrix with zeros and double the number of rows and columns
    typicalMatrix = zeros(2 * m, 2 * n);

    for i = 1:m
        for j = 1:n
            % Extract the tuple from the pair matrix
            tuple = pairMatrix{i, j};

            % Fill in the entries in the typical matrix
            % with the tuple and its reverse
            typicalMatrix(2*i-1, 2*j-1) = tuple(1);
            typicalMatrix(2*i-1, 2*j) = tuple(2);
            typicalMatrix(2*i, 2*j-1) = tuple(2);
            typicalMatrix(2*i, 2*j) = tuple(1);
        end
    end
end