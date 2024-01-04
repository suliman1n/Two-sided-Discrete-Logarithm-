%function result = AttackCSR(A_p, B_q, X, Y)
    
    %A_p = TransformToTypicalMatrix(A_p);
    %B_q = TransformToTypicalMatrix(B_q);
    %X = TransformToTypicalMatrix(X);
    %Y = TransformToTypicalMatrix(Y);
    %M=TropId(size(X,1));
    
    
    %[k, l, tau1] = ModTwoSidedCSR(A_p, X,M, Y);
    
    
    %[r, s, tau2] = ModTwoSidedCSR(B_q, X,M, Y);
    
    %[CX,SX,RX] = CSR(X);
    %[CY,SY,RY] = CSR(Y);
    %lambda1 = MaxCycleMean(X);
    %lambda2 = MaxCycleMean(Y);
    %[CritX,lX] = CriticalCycle(X);
    %[CritY,lY] = CriticalCycle(Y);

    %TEMX1=TropMulti(CX,TropMulti(TropMatPower(SX,mod(k, lX)),RX));
    %TEMX2=TropMulti(CX,TropMulti(TropMatPower(SX,mod(r, lX)),RX));
    %TEMX=TropMulti(TEMX1,TEMX2);

    
    %Xkr = TropMulti(fastpowermaxplus(X, k), fastpowermaxplus(X, r));
    

    %TEMY1=TropMulti(CY,TropMulti(TropMatPower(SY,mod(s, lY)),RY));
    %TEMY2=TropMulti(CY,TropMulti(TropMatPower(SY,mod(l, lY)),RY));
    %TEMY=TropMulti(TEMY1,TEMY2);

    
    %Ysl = TropMulti(fastpowermaxplus(Y, s), fastpowermaxplus(Y, l));
    
    
    %combined = TropMulti(TEMX, TEMY);
    %combined = TropMulti(Xkr, Ysl);
    
    
    %transformed = TransformToPairMatrix(combined);
    
    
    %result = TropicalMultiScalarMatrix(transformed, tau1 + tau2+lambda1*k+lambda1*r+lambda2*s+lambda2*l);
    %result = TropicalMultiScalarMatrix(transformed, tau1 + tau2);
%end




function result = AttackCSR(A_p, B_q, X, Y)
    
    A_pt = TransformToTypicalMatrix(A_p);
    Xt = TransformToTypicalMatrix(X);
    Yt = TransformToTypicalMatrix(Y);
    M=TropId(size(Xt,1));
    
    [k, l, tau1] = ModTwoSidedCSR(A_pt, Xt,M, Yt);
    
    
    temp=PairMatrixMultiply(PairMatrixMultiply(fastpowermaxpluspairt(X,k),B_q),fastpowermaxpluspairt(Y,l));
    

    result = TropicalMultiScalarMatrix(temp, tau1);
    
end