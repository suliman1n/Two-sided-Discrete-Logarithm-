function resultMatrix = TropicalMultiScalarMatrix(A, c)
    % Check if A is a cell array
    if ~iscell(A)
        error('Input matrix A should be a cell array.');
    end

    % Get the dimensions of the input matrix
    [m, n] = size(A);

    % Initialize the result matrix with the same dimensions as A
    resultMatrix = cell(m, n);

    % Add the scalar c to all elements
    for i = 1:m
        for j = 1:n
            % Check if the element is a tuple (numeric array)
            if isnumeric(A{i, j}) && numel(A{i, j}) == 2
                resultMatrix{i, j} = [A{i, j}(1) + c, A{i, j}(2) + c];
            else
                error('Matrix element is not a tuple (numeric array).');
            end
        end
    end
end