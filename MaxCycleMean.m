function lambda = MaxCycleMean(B)
%MaxCycleMean - Copied from Sergey's Library, calculates the eigenvalue/
%max cycle mean of a tropical matrix

[sci, sizes]=TropSCC(B);
noofcomp=max(sci);
eigen=zeros(noofcomp,1);
for h=1:noofcomp
    T=find(sci==h);
    len=length(T);
    for i=1:len
        for j=1:len
            A(i,j)=B(T(i),T(j));
        end
    end
    
 [n,n] = size(A);  
 if n~=1
    for i = 1:n
       K(i,1) = A(i,1);
    end
C = A;
      for k=2:n+1
        A = TropMulti(C,A);  
          for i=1:n
            K(i,k) = A(i,1);
          end
      end
       for i=1:n
            for k=1:n
                 if k==1
                    curmin = (K(i,n+1) - K(i,k))/n;                    
                 else
                    curmin = min(curmin,((K(i,n+1) - K(i,k))/(n+1-k)));
                 end
            end
         if i==1 
             lambda = curmin;
         else
              lambda = max(lambda,curmin);
         end
       end
       eigen(h)=lambda;
 else
     eigen(h)=A;
 end
 A=[];
end
lambda=max(eigen);
end
