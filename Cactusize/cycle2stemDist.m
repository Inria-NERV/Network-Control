function [xDist , xStem , xCycle] = cycle2stemDist(stem , cycle , distances)
nStem = length(stem);
nCycle = length(cycle);

xDist = Inf;
xStem = 0 ; 
xCycle = 0;
for kStem = 1:nStem
    for kCycle = 1:nCycle
        dist = distances(stem(kStem) , cycle(kCycle));
        if dist < xDist
            xDist = dist;
            xStem = stem(kStem);
            xCycle = cycle(kCycle);
        end
    end
end

end