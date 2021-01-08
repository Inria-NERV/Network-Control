function [xTarget,yTarget] = genTargetCircle(xf,eps)
nPoints = 100;
seg = - eps : 2*eps/nPoints : eps;
nPoints = length(seg);
xTarget = [(xf(1)*ones(1,nPoints) + seg) (xf(1)*ones(1,nPoints) + seg)];
yTarget = [(xf(2)*ones(1,nPoints) + sqrt(eps*eps*ones(1,nPoints) - seg.^2)) ...
    (xf(2)*ones(1,nPoints) - sqrt(eps*eps*ones(1,nPoints) - seg.^2)) ];
end