function [tarnsferCapacity , lowdim_transferCapacity] = controlMetric_forStateTransfer(matrix , x0 , xf ,drivers, r_dim , tf)
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
% Function Name: controlMetric_forStateTransfer
% Associated Work: "Low-dimensional controllability of brain networks"
% Reference: Messaoud, R. B., Du, V. L., Kaufmann, B. C., Couvy-Duchesne, B., Migliaccio, L., Bartolomeo, P., ... & Fallani, F. D. V. (2023). Low-dimensional controllability of brain networks. arXiv preprint arXiv:2311.11132.
% Description: Your function description

% Inputs,
% INPUT 1
% matrix = adjacency matrix of the network (nxn matlab matrix)

% INPUT 2
% x0 = nx1 vector = initial state

% INPUT 3
% x0 = nx1 vector = final state


% INPUT 3
% drivers = vector of indices of nodes to to check as drivers and compute
% their control cnetrality ... default  drivers =1:n , check all nodes
% drivers = 1xnDrivers vector
% default : drivers =1:n compute control centreality for all nodes


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

n = size(matrix,1);

if nargin<6
    tf=1;
end

if nargin<5
    r_dim=5;
end
if nargin<4
    drivers=1:n;
end


% We have an expression of energy in finite time but not ininfinte time the
% time parameter does not change the trends

%% %%%%%%%%%%%%%%%%%%%%%
nDrivers = length(drivers) ;


%% normalize matrix to A
lambdaMax = eigs(matrix , 1);
A =  matrix - 1.001 * lambdaMax* eye(n);
% the coef of 1.001 is to ensure Re(Lambda(A))<0 to be stable


%% get the Laplacian eigenmaps of the target network
targetlaplac = diag(sum(matrix)) - matrix;

[VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
lambdaNotSorted = diag(lambdaMatNotSorted);

[lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
V = VnotSorted(: , idxAsc);

orderedModesIdx = 1:n;

%% natural state evolution
evolveState = expm(A*tf)*x0; % natural 

%% loop over drivers
tarnsferCapacity = zeros(nDrivers,1);
lowdim_transferCapacity = zeros(nDrivers,1);

for k = 1:nDrivers

    %% build input matrix B
    driversInd = drivers(k);
    Drivers = zeros(1,n);
    Drivers(driversInd) = 1;
    B = zeros(n);
    B = B - diag(diag(B)) + diag(Drivers);


    %% build controllability Gramian
    sys = ss(A,B,eye(n) , zeros(n,n));
    opt = gramOptions('TimeIntervals',[0 tf]);
    W = gram(sys,'c' , opt);
    W = (W + W')/2;

    %% get transfer capacity : inverse of energy
    tarnsferCapacity(k) =1/( ...
        ((xf - evolveState)')* (W\(xf - evolveState)) ...
                            );

    %% gramian for the eigenmaps of the target (Low-dimensional)
    Ceig = V(:,orderedModesIdx(1 : r_dim))';

    Wbar = Ceig * W * Ceig';

    targetGramSpec = (Wbar + Wbar')/2; % symmetrize just in case
    

    lowdim_transferCapacity(k) =1/(...
        ((Ceig*xf - Ceig*evolveState)')* (targetGramSpec\(Ceig*xf - Ceig*evolveState))...
        );

end

end

