function [rowInd , colInd , weights] = findMatrixDiff(origMat , transMat)
[rowInd,colInd] = find(origMat ~= transMat);
linearidx = sub2ind(size(origMat),rowInd,colInd);
linearMat = origMat(:);
weights = linearMat(linearidx);
end