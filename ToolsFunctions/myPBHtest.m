function pbhTest = myPBHtest(A,B)
n = size(A,1);
pbhTest = 1;
lambda = eig(A);
for l=1:n
    if rank([(lambda(l)*eye(n) - A) B]) < n
        pbhTest = 0;
    end
end
end