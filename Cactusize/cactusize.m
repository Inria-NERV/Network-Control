function [cacti , controlledCacti , cactiMatrix , BMatrix] = cactusize( A )
n = size(A,2);
fprintf('-----  extracting stems and cycles -----\n');

[stems , cycles] = extractStemCycles(A);
nStems = length(stems);
nCycles = length(cycles);
%% Make buds from existing connections
fprintf('-----  making buds 1 -----\n');

[stems , cycles , cycles2StemsMatching , closestStemNodes , closestCycleNodes , unmatchedCycles, ...
    cyclesConflict , stemsLengths,macroDistances , onlyStems] = makeBudsPhase1(stems , cycles , A);
nStems = length(stems);
nCycles = length(cycles);
%% Make buds with new connections
fprintf('-----  making buds 2 -----\n');

cycles2StemsMatching2 = makeBudsPhase2( cycles  , cycles2StemsMatching ,  ...
    unmatchedCycles  , stemsLengths , macroDistances);

%% make a graph structure with only the cacti and the isolated cycles
fprintf('----- constructing Cacti graph -----\n');
[s , t] = makeMyEdges(stems , cycles , cycles2StemsMatching , closestStemNodes , closestCycleNodes);
CactiInit = digraph(s,t);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Deal with  the Plotting stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make a graph structure for the complete controlled cacti
nodesNames = {};
for k =1:n
    nodesNames{k,1} = sprintf('%i',k) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build input matrix 
BMatrix = zeros(n , n);

always = 1; % still may be fix some bugs
if  always  % everything is normal: Then Nd = nStems
    Nd = nStems;             % And input nodes are connected to the biginning of each Stem
    controlledCacti = CactiInit;
    controlledCacti.Nodes.Name = nodesNames;
    for kD = 1:Nd
        controlledCacti = addedge(controlledCacti , sprintf('U%i',kD) , sprintf('%i',stems{kD}(1)) );
        BMatrix(stems{kD}(1) , stems{kD}(1)) = 1;
    end
    
    
    % Relate isolated cycles to corresponding Input Nodes
    [row,col] = find(cycles2StemsMatching2 == 2);
    nEdgesDist = length(row);
    ss = {};
    tt = {};
    for kE = 1:nEdgesDist
        ss{kE} = sprintf('U%i',row(kE));
        tt{kE} = sprintf('%i',closestCycleNodes(row(kE) , col(kE)));
        BMatrix(closestCycleNodes(row(kE) , col(kE)) , stems{row(kE)}(1)) = 1;
    end
    controlledCacti = addedge(controlledCacti , ss,tt );
    controlEdges = Nd + nEdgesDist;
end

figure;
plot(CactiInit);
figure;
layout = 'force';
pp =plot(controlledCacti , 'MarkerSize' , 10 ...
    , 'NodeFontWeight' , 'bold' , 'NodeFontSize' , 10 ,'Layout',layout);
pp.EdgeCData = [zeros(1,controlledCacti.numedges -controlEdges) ones(1,controlEdges)];
pp.NodeCData = [zeros(1,controlledCacti.numnodes - Nd) ones(1,Nd)];
colormap jet
cacti = CactiInit;
% adjacency matrix of cacti
cactiMatrix = full(adjacency(cacti)');

end