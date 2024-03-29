function dist = cosineDist(x,y)
dist = ((x')*y)/(norm(x) * norm(y));
end