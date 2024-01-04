function C= TropAdd(A,B)
%Tropical Addition - Adds two matrices using tropical addition
%   Detailed explanation goes here


[n,m] = size(A);
[n1,m1] = size(B);
if (n~=n1) || (m~=m1) 
    error("Dimension Error!")
end
C = zeros(n,m);
for i=1:n
    for j=1:m
        C(i,j) = max(A(i,j),B(i,j));
    end
end
end

        





