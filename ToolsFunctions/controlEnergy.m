function e = controlEnergy(u)
e = sum(u.^2 , 'all');
end