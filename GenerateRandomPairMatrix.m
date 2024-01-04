function randomMatrix = GenerateRandomPairMatrix(n, minValue, maxValue)
    % Initialize the random matrix as a cell array
    randomMatrix = cell(n, n);
    
    for i = 1:n
        for j = 1:n
            % Generate random integers within the specified range
            a = randi([minValue, maxValue]);
            b = randi([minValue, maxValue]);
            
            % Store the tuple as a numeric array in a cell
            randomMatrix{i, j} = [a, b];
        end
    end
end