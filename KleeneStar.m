function Astar = KleeneStar(A)
%KleeneStar - Calculates the Kleene star of A  
lambda = MaxCycleMean(A);
if  abs(lambda) >0.0005 
    error("This Matrix will not converge as it has a non zero eigenvalue")
end

[n,n]=size(A);
for i=1:n
    for j=1:n
B(i,j)=A(i,j)-lambda;
    end
end

for i=1:n
    for j=1:n
        if i==j
             I(i,j)=0;
        else
             I(i,j)=-inf;
        end
    end
end

G=B;
g=B;
for i=1:n
    G=TropMulti(B,G);
    g=TropAdd(G,g);
    i=i+1;
end
gamma=g;
Astar=TropAdd(I,gamma);





end