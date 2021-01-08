function s = localControl(x,x0,xf, localEps)

if localEps =='auto'
    localEps = norm(xf - x0)/100;
end
nSamps = size(x,1);
norms = zeros(nSamps,1);
for k = 1:nSamps
    norms(k) = norm(x(k,:)' - x0);
end
s = (max(norms) < 2*norm(xf - x0) + localEps);
end