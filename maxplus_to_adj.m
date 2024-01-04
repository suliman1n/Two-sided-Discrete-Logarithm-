function C=maxplus_to_adj(A)
[m,n]=size(A);
for i=1:m
    for j=1:n
        if A(i,j)==-inf 
            C(i,j)=0;
       else
            C(i,j)=1;
        end
    end
end