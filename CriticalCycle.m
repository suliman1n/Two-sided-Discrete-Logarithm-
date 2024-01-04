function [C,length] = CriticalCycle(A,index)
%CriticalCycle - Finds a critical cycle of graph , starting a node 1
    if ~exist('index','var')
       index = 1 
    end 
[n,n] = size(A);
C = zeros(n)+-inf;
Critnodes = [];
[a,b] = howardmaxplus(A);
history=index;

     
[m,i] = max(A(index,:)+b');
history=[history,i];

while size(unique(history))==size(history)
    [m,i]=max(A(i,:)+b');
    history=[history,i];
end
start = find(history==i);
start = start(1);
Critnodes = history(start:end-1);
[d,l]=size(Critnodes); 

for j=1:(l-1)
    C(Critnodes(j),Critnodes(j+1))= A(Critnodes(j),Critnodes(j+1));
end
     C(Critnodes(end),Critnodes(1))= A(Critnodes(end),Critnodes(1));
length = l;
end

    



