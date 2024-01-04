function [C,S,R] = CSR(A)
%CSR -  Calculates the CSR terms of the tropical matrix A 
%   Detailed explanation goes here
% might include a second input as the critical graph, as can obviously gain
% different CSR expansions from different critical cycles, for now ill keep
% it like this until we implement the critical graph part
[n,n] = size(A);
lambda = MaxCycleMean(A);
[Crit,L] = CriticalCycle(A);
Abar= TropMatPower(A - lambda,L);
Ahat = KleeneStar(Abar);


C_identifier = TropId(n);
R_identifier = TropId(n);
C  =zeros(n);
R= zeros(n);
S = Crit - lambda;
%First we need to find a critical cycle, cant remember how we do this so
%just assume we have it, in the form of a matrix of the size of A, where a
%1 indicates the edge is in this critical cycle and a 0 indicates that it
%is not
% we could then also illustrate the critical graph as an actual graph of a
% matrix with the element, and then -inf if its not in the critical graph.


for i =1:n
    if max(Crit(:,i))==-inf  
        C_identifier(i,i) = -inf;
    end
    if max(Crit(i,:))==-inf
        R_identifier(i,i) = -inf;
    end
end
C = TropMulti(Ahat,C_identifier);
R = TropMulti(R_identifier,Ahat);


% Just a step to make sure there are not any -inf 
 p = find(isnan(R))';
 R(p) = -inf;

  q = find(isnan(C))';
 C(q) = -inf;



end