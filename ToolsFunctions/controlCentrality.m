function [co , cc]  = controlCentrality(A)
n = size(A,1);
cc = zeros(n,1);
for k=1:n
    B = zeros(n);
    B(k,k)=1;
    C = ctrb(A,B);
    cc(k)= rank(C);
end
[B,I] = sort(cc , 'descend');
cno = 1:n;
co = cno(I);
end