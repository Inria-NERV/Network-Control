function lowCC= low_dimensional_control_centrality(matrix , target, drivers, r_dim)
%--------------------------------------------------------------------------
% Copyright (c) [2023] [Remy Ben messaoud]
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
% Function Name: low_dimensional_control_centrality
% Associated Work: "Low-dimensional controllability of brain networks"
% Reference: Messaoud, R. B., Du, V. L., Kaufmann, B. C., Couvy-Duchesne, B., Migliaccio, L., Bartolomeo, P., ... & Fallani, F. D. V. (2023). Low-dimensional controllability of brain networks. arXiv preprint arXiv:2311.11132.
% Description: Your function description

% Inputs,
% INPUT 1
% matrix = adjacency matrix of the network

% INPUT 2
% target = vector of indices of nodes composing the target network target =
% default : target =1:n control the entire network
% target = 1xtargetsize vector

% INPUT 3
% drivers = vector of indices of nodes to to check as drivers and compute
% their control cnetrality ... default  drivers =1:n , check all nodes
% drivers = 1xnDrivers vector

% INPUT 4
% r_dim = the low-dimension of the output control = number of eigenmaps of
% the targeted network to control... default r_dim=5;

% Outputs,
% OUTPUT 1
% lowCC = vector of size 1xnDrivers that returns the low-dimensional
% control centrality for each driver

%% %%%%%%%%%%%%%%%%%%%%%
nDrivers = length(drivers) ;
n = size(A,1);

targetsize = length(target);

%% normalize matrix to A
lambdaMax = eigs(matrix , 1);
A =  matrix - 1.001 * lambdaMax* eye(n);
% the coef of 1.001 is to ensure Re(Lambda(A))<0 to be stable

%% isolate target network
%%%%% big C matrix for target Control %%%%%%%
bigC = zeros(targetSize , n);
for ktarget = 1 : targetSize
    bigC(ktarget , target(ktarget)) = 1;
end
targetnet = matrix(target , target);

%% get the Laplacian eigenmaps of the target network
targetlaplac = diag(sum(targetNet)) - targetNet;

[VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
lambdaNotSorted = diag(lambdaMatNotSorted);

[lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
V = VnotSorted(: , idxAsc);

%% loop over drivers
lowCC = zeros(1,nDrivers);
for k = 1:nDrivers

    %% build input matrix B
    driversInd = drivers(k);
    Drivers = zeros(1,n);
    Drivers(driversInd) = 1;
    B = zeros(n);
    B = B - diag(diag(B)) + diag(Drivers);


    %% build controllability Gramian
    sys = ss(A,B,eye(n) , zeros(n,n));
    W = gram(sys,'c');
    %% gramian for the target
    W = (W + W')/2;
    targetGram = (bigC * W *(bigC'));

    %% gramian for the eigenmaps of the target
    Ceig = V(:,orderedModesIdx(1 : r_dim))';

    Wbar = Ceig * targetGram * Ceig';

    targetGramSpec = (Wbar + Wbar')/2;
    lambdaSpecRed = eig(targetGramSpec);
    [~ , ismallabs] = min(abs(lambdaSpecRed));

    lowCC(k) = lambdaSpecRed(ismallabs);
end
end

