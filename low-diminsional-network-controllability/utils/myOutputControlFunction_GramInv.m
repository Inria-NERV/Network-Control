function [ySimu, xSimu, U ] = myOutputControlFunction_GramInv(A, B, T, t, x0, yf , C , W)
n =size(A,1);
nSamps = length(t);

%
sys = ss(A,B, eye(n), 0);
sysObs = ss(A,B,C,0);

if nargin < 9
    isStable = isstable(sys);
    if isStable == 1
        %%%  finite time gramian
        opt = gramOptions('TimeIntervals' , [0 T]);
        W = gram(sys , 'c' , opt);

        %%% infinite time gramian
        %     W = gram(A , B);

    else
        [W,Tr] = CtrGram(A,B);
    end
end

%% 
% [W,Tr] = CtrGram(A,B);
targetGram = (C*W*(C'));

% %% in loop
% ticT = tic;
% Winv = inv(W);
U = zeros(nSamps , n);
% X_opt = zeros(nSamps , n);
P_adj  = nan;

gramInvMat = (C')*(targetGram\(   C*(-expm((T)*A)*x0) + yf)    );

for kSamp = 1:nSamps
%     Wt = myGram(A, B, t(kSamp),STEP);
    U(kSamp , :) = ((B')*expm((T-t(kSamp))*A')*gramInvMat)' ;
end

U(isnan(U)) = 0;
U(U == Inf) = 0;



%% simulate

[ySimu,tt,xSimu] = lsim(sysObs,U/(1),t,x0);

end