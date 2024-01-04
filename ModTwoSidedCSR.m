function [t_1,t_2,alpha] = ModTwoSidedCSR(U,A,M,B)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[n,m] = size(A);
[~,SA,RA] = CSR(A);
[CB,SB,~] = CSR(B)
lambda1 = MaxCycleMean(A);
lambda2 = MaxCycleMean(B);
[CritA,lA] = CriticalCycle(A)
[CritB,lB] = CriticalCycle(B)
t_1 = -1 
t_2 = -1 
alpha = 0 


for i = 1:lA
    for j =1:lB
        SARA = TropMulti(fastpowermaxplus(SA,i),RA)
        CBSB = TropMulti(CB,TropMatPower(SB,j))
        SRMCS = TropMulti(TropMulti(SARA,M),CBSB)
        temp = U - SRMCS
       checker = []

     for k=1:n 
            for d =1:n
                if max(CritA(k,:))~=-inf && max(CritB(:,d))~=-inf
                checker = [checker,temp(k,d)]
                end
            end
     end
        if all(abs(checker-checker(1))<0.000000005)
            r = checker(1)
            i 
            j
            i_ = rem(i,lA)
            j_ = rem(j,lB)
        

            
            T =  r - i_*lambda1 - j_*lambda2;
            
          

            A_ = lA*lambda1;
            B_ = lB*lambda2;
            x_bound = ((n-1)*lA-i_)/lA;
            y_bound = ((n-1)*lB-j_)/lB;

LB = [10*x_bound;10*y_bound;-inf]
UB = [inf;inf;inf]

            f = [1;1;0];
            intcon = [ 1 2];
            Aa = [ -1 0 0  ; 0 -1 0  ; -1 0 0  ; 0 -1 0 ] ;
            bb = [-2*x_bound;-2*y_bound;0;0];
            Aeq = [A_,B_,1];
            beq = T;
            [x1] = intlinprog(f,intcon,Aa,bb,Aeq,beq,LB,UB);
            x_ = x1(1)
            y_ = x1(2)
            alpha = x1(3)
       
            t_1 = x_*lA + i_
            t_2 = y_*lB + j_ 
            t_1*lambda1 + t_2*lambda2 + alpha 
            return

        end 


    end 
end
end 