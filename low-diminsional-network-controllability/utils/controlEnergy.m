function e = controlEnergy(u , dt)
[nSamps , n] = size(u);
e = 0;
for k= 1 :nSamps
    e = e + (norm(u(k,:))^2) * dt;
end
% e = e/nSamps;
end