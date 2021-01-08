function cycles2StemsMatching2 = makeBudsPhase2( cycles  , cycles2StemsMatching ,  ...
      unmatchedCycles  , stemsLengths , macroDistances)
%% This function creates new connections between the unmatched cycles and the input nodes
nCycles = length(cycles);
[~,shortestStem] = min(stemsLengths);
nRemainingCycles = length(unmatchedCycles);
cycles2StemsMatching2 = cycles2StemsMatching;
for kC = 1:nRemainingCycles
    % find real cycle index
    kReal = unmatchedCycles(kC);
%     for kk = 1:nCycles
%         if isequal(cycles{kk} , unmatchedCycles(kC))
%             kReal = kk;
%         end
%     end
    % find closest stem~
    [bb,closeStem] = min(macroDistances(:,kReal)');
    % check if the cycle is isolated , then relate it to the shortest stem
    if bb == Inf
        cycles2StemsMatching2(shortestStem , kReal) = 2 ; 
        % we mark it with 2 to differentiate the new connection
    else
        cycles2StemsMatching2(closeStem , kReal) = 2 ; 
    end
end
end