function L = trajLength(x)
[nSamp, n] = size(x);
x0 = x(1,:);
L=0;
for kSamp=2:nSamp
    L = L + norm(x(kSamp , :) - x0);
    x0 = x(kSamp,:);
end
end