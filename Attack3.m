%function result = Attack3(A_p, B_q, X, Y)
    
    %A_pr = ReduceMatrix(A_p);
    %B_qr = ReduceMatrix(B_q);
    %Xr = ReduceMatrix(X);
    %Yr = ReduceMatrix(Y);
    %M=TropId(size(Xr,1));
    
    
    %[k, l, tau1] = ModTwoSidedCSR(A_pr, Xr,M, Yr);
    
    
    %[r, s, tau2] = ModTwoSidedCSR(B_qr, Xr,M, Yr);

    %X = TransformToTypicalMatrix(X);
    %Y = TransformToTypicalMatrix(Y);

    %[CX,SX,RX] = CSR(X);
    %[CY,SY,RY] = CSR(Y);
    %lambda1 = MaxCycleMean(X);
    %lambda2 = MaxCycleMean(Y);
    %[CritX,lX] = CriticalCycle(X);
    %[CritY,lY] = CriticalCycle(Y);

    %TEMX1=TropMulti(CX,TropMulti(fastpowermaxplus(SX,mod(k, lX)),RX));
    %TEMX2=TropMulti(CX,TropMulti(fastpowermaxplus(SX,mod(r, lX)),RX));
    %TEMX=TropMulti(TEMX1,TEMX2);
    

    
    %Xkr = PairMatrixMultiply(fastpowermaxpluspairt(X, k), fastpowermaxpluspairt(X, r));
    


    
    %Ysl = PairMatrixMultiply(fastpowermaxpluspairt(Y, s), fastpowermaxpluspairt(Y, l));

    %TEMY1=TropMulti(CY,TropMulti(fastpowermaxplus(SY,mod(s, lY)),RY));
    %TEMY2=TropMulti(CY,TropMulti(fastpowermaxplus(SY,mod(l, lY)),RY));
    %TEMY=TropMulti(TEMY1,TEMY2);
    
    
    %combined = PairMatrixMultiply(Xkr, Ysl);
    %combined = TropMulti(TEMX, TEMY);
    
    
    %transformed = combined;
    %transformed = TransformToPairMatrix(combined);
    
    
    %result = TropicalMultiScalarMatrix(transformed, tau1 + tau2);
    %result = TropicalMultiScalarMatrix(transformed, tau1 + tau2+lambda1*k+lambda1*r+lambda2*s+lambda2*l);
%end



function result = Attack3(A_p, B_q, X, Y)
  
    A_pr = ReduceMatrix(A_p);
    Xr = ReduceMatrix(X);
    Yr = ReduceMatrix(Y);
    M=TropId(size(Xr,1));
    
    [k, l, tau1] = ModTwoSidedCSR(A_pr, Xr,M, Yr);
    
    temp=PairMatrixMultiply(PairMatrixMultiply(fastpowermaxpluspairt(X,k),B_q),fastpowermaxpluspairt(Y,l));
    
    result = TropicalMultiScalarMatrix(temp, tau1);
end