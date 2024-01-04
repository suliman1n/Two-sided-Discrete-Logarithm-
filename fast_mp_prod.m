function C = fast_mp_prod(A,B)
n = size(A,1);
C = -inf*ones(n);
A = transpose(A);
for i = 1:n
    C(i,:) = max(A(:,i) + B);
end