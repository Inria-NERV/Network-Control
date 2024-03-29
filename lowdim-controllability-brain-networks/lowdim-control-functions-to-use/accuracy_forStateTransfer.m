function [accuracy , trajectory_length] = accuracy_forStateTransfer(matrix_path , x0_path, xf_path , target, drivers , r_dim )
%--------------------------------------------------------------------------
% Copyright (c) [2024] [Remy Ben messaoud]
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%--------------------------------------------------------------------------

% Function Header
% Author: Remy Ben Messaoud
% Function Name: accuracy_forStateTransfer
% Associated Work: "Low-dimensional controllability of brain networks"
% Reference: Messaoud, R. B., Du, V. L., Kaufmann, B. C., Couvy-Duchesne, B., Migliaccio, L., Bartolomeo, P., ... & Fallani, F. D. V. (2023). Low-dimensional controllability of brain networks. arXiv preprint arXiv:2311.11132.
% Description: Your function description

% Inputs,
% INPUT 1
% matrix_path = adjacency matrix of the network (nxn matlab matrix) (path
% to csv file)

% INPUT 2
% x0_path = nx1 vector = initial state (path to csv file)

% INPUT 3
% xf_path = nx1 vector = final state (path to csv file)


% INPUT 3
% drivers = vector of indices of nodes to to check as drivers and compute
% their control cnetrality ... default  drivers =1:n , check all nodes
% drivers = 1xnDrivers vector
% default : drivers =1:n compute control centreality for all nodes

%%%% optinal arguments
% INPUT 4
% r_dim = the low-dimension of the output control = number of eigenmaps of
% the targeted network to control... default r_dim=5;

% INPUT 4
% tf = time allowed for state transfer (default =1, does not change the
% trends) ... We cannot calculate the enrgy in the infinite time case


% Outputs,
% OUTPUT 1
% tarnsferCapacity = vector of size nDrivers x 1 that returns the capacity
% for each driver to achieve the state transfer  (inverse of the energy)
% !!!! There can be negatve values because of numerical errors, (negative values have to
% be discarded and should be zero in principle

% OUTPUT 2
% lowdim_transferCapacity = vector of size nDrivers x 1 that returns the capacity
% for each driver to achieve the state transfer  (inverse of the energy)
% but in the low-dimensional space (dim=5)
% !!!! Even here there can be some residual negative values because of numerical
% errors, (negative values have to
% be discarded and should be zero in principle

if nargin<6
    r_dim=5;
end

if nargin<5
    drivers=1:n;
end

if nargin<4
    target=1:n;
end

%% %%%%%%%%%%%%%%%%%%%%%
nDrivers = length(drivers) ;
targetSize = length(target);

%% load connectome
matrix = table2array(readtable(matrix_path));

% add upper part
matrix = matrix+ matrix';

%%% get only  parcels of interest

goodIdx = [1:100]; % this valid for Scaefer100 % change as needed

matrix = matrix(goodIdx , goodIdx);
n = size(matrix,2);


matrix = matrix - diag(diag(matrix));

matrix = matrix/max(max(matrix));

%% load initial and final states
removeStateMean = 1;

x0_raw = table2array(readtable(x0_path));

xf_raw = table2array(readtable(xf_path));


%%% remove mean state
if removeStateMean
    x0 = x0_raw - mean([x0_raw]);
    xf = xf_raw - mean([x0_raw]);
else
    x0 = x0_raw;
    xf = xf_raw ;
end

%% Control options
% constatns
T = 8; %in s
STEP = 0.01 ; %in s
t = 0:STEP:T;
nSamps = length(t);
R = eye(n,n);
Q = zeros(n,n);
rhoOpt = 10^(-6);


%% normalize matrix to A
lambdaMax = eigs(matrix , 1);
A =  matrix - 1.001 * lambdaMax* eye(n);
% the coef of 1.001 is to ensure Re(Lambda(A))<0 to be stable

%% isolate target network

%%%%% big C matrix for target Control %%%%%%%
Ctarget = zeros(targetSize , n);
for ktarget = 1 : targetSize
    Ctarget(ktarget , target(ktarget)) = 1;
end
targetNet = matrix(target , target);

%% get the Laplacian eigenmaps of the target network
orderEigenmaps = 1;


targetlaplac = diag(sum(targetNet)) - targetNet;
targetlaplac = (targetlaplac + targetlaplac')/2; % symmetrize just in case

[VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
lambdaNotSorted = diag(lambdaMatNotSorted);

[lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
V = VnotSorted(: , idxAsc);

if orderEigenmaps
    xfGFT = (V')*(Ctarget* xf );
    x0GFT = (V')*(Ctarget* x0);
    [lambdaMag , idxGFTxf] = sort(abs(abs(xfGFT) - abs(x0GFT))./(abs(xfGFT) + abs(x0GFT)) , 'descend');

    % [lambdaMag , idxGFTxf] = sort((abs(xfGFT) ), 'descend');

    orderedModesIdx =idxGFTxf; % descend
else
    orderedModesIdx =1:n;
end

%% loop over drivers
accuracy = zeros(nDrivers,1);
trajectory_length = zeros(nDrivers,1);


C = (V(:,orderedModesIdx(1 : r_dim))')*Ctarget ;

yf = C * xf;


for k = 1:nDrivers

    %% build input matrix B
    driversInd = drivers(k);
    Drivers = zeros(1,n);
    Drivers(driversInd) = 1;
    B = zeros(n);
    B = B - diag(diag(B)) + diag(Drivers);

    %% calculate trajectory with driver
    [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] = ...
    myOutputControlFunction(A, B, T, t, x0, xf, rhoOpt , R , Q , C , yf  );

    %% evaluate performance of driver
    accuracy(k,1) =  cosineSimilarity(ySimu(end,:)' , yf);
    trajectory_length(k,1) = trajLength(xSimu);

end


end

function acc = cosineSimilarity(x,y)
acc = ((x')*y)/(norm(x) * norm(y));
end

function L = trajLength(x , dimNorm)
if nargin < 2
    dimNorm = 2;
end
dimNorm = max([2 dimNorm]);

[nSamp, n] = size(x);
xPrev = x(1,:);
L=0;
for kSamp=2:nSamp
    L = L + norm(x(kSamp , :) - xPrev  , dimNorm);
    xPrev = x(kSamp,:);
end
% L = L/T;
end
