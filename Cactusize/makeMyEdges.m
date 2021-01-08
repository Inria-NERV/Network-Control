function [s , t] = makeMyEdges(stems , cycles , cycles2StemsMatching , closestStemNodes , closestCycleNodes)
nStems = length(stems);
nCycles = length(cycles);
s = [];
t = [];
%% first build the stems and cycles 
% stems
for kStem = 1:nStems
    stem = stems{kStem};
    stemLength = length(stem);
    if stemLength > 1
        for kNode = 1:stemLength-1
            s = [s stem(kNode)];
            t = [t stem(kNode+1)];
        end
    else
        fprintf('stem with only one node !!!\n')
    end
end
% cycles
for kCycle = 1:nCycles
    cycle = cycles{kCycle};
    cycleLength = length(cycle);
    if cycleLength > 1
        for kNode = 1:cycleLength-1
            s = [s cycle(kNode)];
            t = [t cycle(kNode+1)];
        end
        s = [s cycle(end)];
        t = [t cycle(1)];
    else
        fprintf('cycle with only one node !!!\n')
    end
end
%% Then we deal with the distinctive edges
% existing edges 
[row,col] = find(cycles2StemsMatching == 1);
nEdgesDist = length(row);
for kE = 1:nEdgesDist
    s = [s closestStemNodes(row(kE) , col(kE))];
    t = [t closestCycleNodes(row(kE) , col(kE))];
end

end