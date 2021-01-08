%% Function that gives the minimum nuber of driver nodes as the maximum geometric multiplicity

function nd = minNumDrivers(A)
N = size(A,1);
geoMul = zeros(N,1);
lamAlreadySeen = NaN(N,1);
lambda = eig(A);
for i=1:N
    lam = lambda(i);
	if sum(abs(lamAlreadySeen - lam*ones(N,1))< 10^-8) == 0
        lamAlreadySeen(i) = lam;
		geoMul(i) = N - rank(lam*eye(N)-A);
	end
end

[nd,I] = max(geoMul);
end