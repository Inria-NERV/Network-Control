function [co , cc]  = heuristiControlCentrality(A)
n = size(A,1);
cc = zeros(n,1);
G = digraph(A');
dout = outdegree(G);
din = indegree(G);
for k=1:n
    B = zeros(n);
    B(k,k)=1;
    C = ctrb(A,B);
    cc(k)= rank(C);
end
[B,I] = sort(cc + dout - din  , 'descend');
%[B,I] = sort(cc + n/2*(dout +1)/(din +1)  , 'descend');
cno = 1:n;
co = cno(I);
end