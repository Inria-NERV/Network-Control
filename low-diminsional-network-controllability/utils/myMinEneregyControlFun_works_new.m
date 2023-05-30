function [ySimu, xSimu, U ] = myMinEneregyControlFun_works(A, B, T,STEP, t, x0, yf , C , W)
n =size(A,1);
nSamps = length(t);
% sys_xp = ss(A,B,eye(n),0);
% [X_opt, U_opt, n_err , P_adj] = optim_fun(A, T,STEP, B, x0, xf, rho, S);
% 
% U = reshape(U_opt , [n , nSamps])';
% 
% [Y,tt,xSimu] = lsim(sys_xp,U,t,x0);
% 
% [ xOpenLoopbis ] = open_loop_control( A, B, x0, U_opt , STEP);
% %[ xSelfDynamicBis ] = open_loop_control( A, B, x0, zeros(size(U_opt)) , STEP);

%xfSimu = squeeze( xOpenLoopbis(:,1,end)) ;
% W = myGram(A, B, T,STEP);

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

%     if isnan( Ut(1)) || Ut(1) == inf
%         U(kSamp , :) = zeros(1,n);
%     else
%         U(kSamp , :) = Ut;
%     end


%     Xt = (expm(t(kSamp)*A)*(x0 + Wt * W\(expm((-T)*A)*x0 - xf)))';
%     if isnan( Xt(1)) || Xt(1) == inf
%         X_opt(kSamp , :) = zeros(1,n);
%     else
%         X_opt(kSamp , :) = Xt;
%         %X_opt(kSamp , :) = zeros(1,n);
%     end
% end
% tacT = toc;
% fprintf('loop calcul = %d\n' , tacT-ticT);
% 
% %% try closed form for U 
% ticT = tic;

%% simulate

[ySimu,tt,xSimu] = lsim(sysObs,U/(1),t,x0);

end