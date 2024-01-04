function A = TropId(n)
%TropIdent - Creates a Tropical Identity matrix of size n 
%   Detailed explanation goes here
A = zeros(n);
A(1:n,1:n) = -inf;
for i=1:n
    A(i,i)=0;
end