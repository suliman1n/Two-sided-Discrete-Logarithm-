function A = fastpowermaxplus(B, t)
    [n, m] = size(B);
    
    if (n ~= m)
        error("Dimension Error! Not a square matrix")
    end
    
    if t == 0
        A = TropId(n);
        return
    end 
    
    exp = dec2bin(t);
    value = TropId(n); 
    
    for i = 1:length(exp)
        value = TropMulti(value, value);
        
        if exp(i) == '1'
            value = TropMulti(value, B);
        end
    end
    
    A = value;
end
