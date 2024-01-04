%function result = Attack(A_p, B_q, X, Y, max_t)
    
    %A_p = TransformToTypicalMatrix(A_p);
    %B_q = TransformToTypicalMatrix(B_q);
    %X = TransformToTypicalMatrix(X);
    %Y = TransformToTypicalMatrix(Y);
    
    
    %[k, l, tau1] = ModTwoSidedTrialandErrorpairs(A_p, X, Y, max_t);
    
    
    %[r, s, tau2] = ModTwoSidedTrialandErrorpairs(B_q, X, Y, max_t);
    
    
    %Xkr = TropMulti(fastpowermaxplus(X, k), fastpowermaxplus(X, r));
    
    
    %Ysl = TropMulti(fastpowermaxplus(Y, s), fastpowermaxplus(Y, l));
    
    
    %combined = TropMulti(Xkr, Ysl);
    
    
    %transformed = TransformToPairMatrix(combined);
    
    
    %result = TropicalMultiScalarMatrix(transformed, tau1 + tau2);
%end


function result = Attack(A_p, B_q, X, Y, max_t)
    
    A_pt = TransformToTypicalMatrix(A_p);
    
    Xt = TransformToTypicalMatrix(X);
    Yt = TransformToTypicalMatrix(Y);
    
    
    [k, l, tau1] = ModTwoSidedTrialandErrorpairs(A_pt, Xt, Yt, max_t);
    
    
    
    
    temp=PairMatrixMultiply(PairMatrixMultiply(fastpowermaxpluspairt(X,k),B_q),fastpowermaxpluspairt(Y,l));
    
    
    
    result = TropicalMultiScalarMatrix(temp, tau1);
end