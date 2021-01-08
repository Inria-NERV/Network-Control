function [stems,cycles] = extractStemCycles(A)
n = size(A,2);
nodes = 1:n;

[matched, unmatched, matchedBy ] = maximumMatching( A' );
Nd = length(unmatched);
%% get list of edges in the Matching
M = [];
for k=1:n
    if matchedBy(k)~= 0
        M = [M ;matchedBy(k) k];
    end
end

%% Get Stems
nStems = Nd;
stems = {};
nodesInStems = [];
for kStems = 1:nStems
    x = unmatched(kStems);
    hasMate = 1;
    stem = [x];
    nodesInStems = [nodesInStems x];
    stemLength = 1;
    while hasMate
        % find child in the matching
        child = find(matchedBy == x);
        stem = [stem child];
        nodesInStems = [nodesInStems child];
        x = child;
        hasMate = ismember(x, matchedBy);
    end
    stems{kStems} = stem;
end
%% Get Cycles
W = setdiff(nodes , nodesInStems); % nodes that are not in stems
nW = length(W);
visited = zeros(1,nW);
cycles = {};
kCycle = 0;
for kc = 1:nW
    if ~visited(kc)
        kCycle = kCycle +1;
        v = W(kc);
        x = v;
        visited(kc)=1;
        notEnd = 1;
        cycle = [v];
        while notEnd
            child = find(matchedBy == x); %get child
            cycle = [cycle child];
            kChild = find(W == child);
            visited(kChild)=1;
            grandChild = find(matchedBy == child);
            notEnd = grandChild ~= v ;
            x = child ; 
            %bleble = 8
        end
        cycles{kCycle} = cycle;
    end
end

end