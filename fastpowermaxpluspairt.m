function A = fastpowermaxpluspairt(B, t)
    [n, m] = size(B);
    
    if (n ~= m)
        error("Dimension Error! Not a square matrix")
    end
    
    if t == 0
        A = TropId(2*n);
        A=TransformToPairMatrix(A);
        return
    end 
    
    exp = dec2bin(t);
    value = TropId(2*n);
    value=TransformToPairMatrix(value);
    
    for i = 1:length(exp)
        value = PairMatrixMultiply(value, value);
        
        if exp(i) == '1'
            value = PairMatrixMultiply(value, B);
        end
    end
    
    A = value;
end
