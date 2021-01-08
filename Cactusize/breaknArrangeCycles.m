function [stems , newCycles , nStems, nCycles] = breaknArrangeCycles(cycles , cycle2break , chosenDriver)
% make stem   
myCycle = cycles{cycle2break};
% verification 
if ~ismember(chosenDriver , myCycle)
    fprintf('error when breaking cycle into stem\n')
end
driverInd = find(myCycle == chosenDriver);
stem = [myCycle(driverInd : end) myCycle(1: driverInd - 1) ];
stems = {};
stems{1} = stem;
nStems = 1;
% Now rearrange the remaining cycles
nCycles = length(cycles);
if nCycles == 1
    newCycles = {};
    nCycles = 0;
else
    newCycles = {};
    kk = 1;
    for k = 1:nCycles
        if k ~= cycle2break
            newCycles{kk} = cycles{k};
            kk = kk + 1;
        end
    end
nCycles = nCycles - 1;
end
end