function [stems , cycles , cycles2StemsMatching , closestStemNodes , closestCycleNodes , unmatchedCycles, ...
    cyclesConflict , stemsLengths,macroDistances , onlyStems] = makeBudsPhase1(stems , cycles , A)
nStems = length(stems);
nCycles = length(cycles);
stemsLengths = zeros(1,nStems);
% find shorstest stem
for kStem=1:nStems
    stemsLengths(1,kStem) = length(stems{kStem});
end
[M,shortestStem] = min(stemsLengths);

% check if the graph is connected
bins = conncomp(digraph(A) ,'Type',  'weak');
if ~isempty(find(bins ~=1))
    fprintf('the graph is not connected \n')
end

%% case 1 there are no cycles
onlyStems = 0;
if nCycles == 0
    fprintf('no cycles in the cacti\n')
    macroDistances = Inf(nStems, nCycles); %distance from stem i to cycle j
    closestStemNodes = zeros(nStems, nCycles); %closestStemNodes for stem i to cycle j
    closestCycleNodes = zeros(nStems, nCycles); %closestCycleNodes from stem i to cycle j
    cycles2StemsMatching = zeros(nStems,nCycles);
    visitedCycles = zeros(1,nCycles);
    % check for conflicts : cycle with distance 1 to many stems
    cyclesConflict = zeros(1,nCycles);
    unmatchedCycles = [];
    onlyStems = 1;
%% case 2 there are both
elseif nCycles ~= 0
    % first deal with special case where there are no Stems
    if nStems == 0 % this is the case 3: break one the cycles and turn it into a stem
        % Then do the same thing as case 2
        [ CC ] = controlCentrality( A, ones(size(A,1),1) );
        % find node with biggest control centrality and make it a Driver
        [~,chosenDriver] = max(CC);
        % go through cycles
        cycle2break = 0;
        for kCycle =1:nCycles
            if ismember(chosenDriver , cycles{kCycle})
                cycle2break = kCycle;
            end
        end
        % break the cycle
        [stems , newCycles , nStems, nCycles] = breaknArrangeCycles(cycles , cycle2break , chosenDriver);
        cycles = newCycles;
        if nCycles == 0
            fprintf('no cycles left in the cacti may be some bugs\n')
            onlyStems = 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = digraph(A');
    d = distances(G);
    macroDistances = Inf(nStems, nCycles); %distance from stem i to cycle j
    closestStemNodes = zeros(nStems, nCycles); %closestStemNodes for stem i to cycle j
    closestCycleNodes = zeros(nStems, nCycles); %closestCycleNodes from stem i to cycle j
    if onlyStems == 1
        cycles2StemsMatching = zeros(nStems,nCycles);
        visitedCycles = zeros(1,nCycles);
        % check for conflicts : cycle with distance 1 to many stems
        cyclesConflict = zeros(1,nCycles);
        unmatchedCycles = [];
        return
    end
    
    for kStem=1:nStems
        for kCycle = 1:nCycles
            [xDist , xStem , xCycle] = cycle2stemDist(stems{kStem} , cycles{kCycle} , d);
            macroDistances(kStem , kCycle) = xDist;
            closestStemNodes(kStem , kCycle) = xStem;
            closestCycleNodes(kStem , kCycle) = xCycle;
        end
    end
    cycles2StemsMatching = zeros(nStems,nCycles);
    visitedCycles = zeros(1,nCycles);
    % check for conflicts : cycle with distance 1 to many stems
    cyclesConflict = zeros(1,nCycles);
    if nStems > 1
        for kCycle = 1:nCycles
            cyclesConflict(1,kCycle) = isequal(macroDistances(:,kCycle) , ones(nStems , 1));
        end
    end
    % make buds
    for kStem=1:nStems
        % check for already existing buds
        closeCycles = find(macroDistances(kStem,:) ==1);
        nCloseCycle = length(closeCycles);
        for kCloseCycle = 1:nCloseCycle
            if cyclesConflict(1,kCloseCycle) == 1 % check for confict
                if ~visitedCycles(1,kCloseCycle)
                    % assign cycle to the shorstest stem
                    cycles2StemsMatching(shortestStem , kCloseCycle) = 1;
                    visitedCycles(1 , kCloseCycle) = 1;
                end
            else
                if ~visitedCycles(1,kCloseCycle)
                    % assign cycle to the shorstest stem
                    cycles2StemsMatching(kStem , kCloseCycle) = 1;
                    visitedCycles(1 , kCloseCycle) = 1;
                end
            end
        end
    end
    
    % find unmatched cycles
    unmatchedCycles = [];
    for kCycle = 1:nCycles
        if isempty(find(cycles2StemsMatching(:,kCycle)))
            unmatchedCycles = [unmatchedCycles kCycle];
        end
    end
end
end
