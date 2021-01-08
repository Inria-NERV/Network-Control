function [X_opt, xSimu, U , P_adj] = myMinEneregyControlFun(A, B, T,STEP, t, x0, xf)
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
W = myGram(A, B, T,STEP);
% Winv = inv(W);
U = zeros(nSamps , n);
X_opt = zeros(nSamps , n);
P_adj  = nan;
for kSamp = 1:nSamps
    Wt = myGram(A, B, t(kSamp),STEP);
    Ut = ((B')*exp((T-t(kSamp))*A')*W\(-exp((T-t(kSamp))*A)*x0 + xf))' ;
    Xt = (exp(t(kSamp)*A)*(x0 + Wt * W\(exp((-T)*A)*x0 - xf)))';
    if isnan( Ut(1)) || Ut(1) == inf
        U(kSamp , :) = zeros(1,n);
    else
        U(kSamp , :) = Ut;
    end
    if isnan( Xt(1)) || Xt(1) == inf
        X_opt(kSamp , :) = zeros(1,n);
    else
        %X_opt(kSamp , :) = Xt;
        X_opt(kSamp , :) = zeros(1,n);
    end
end
sys_xp = ss(A,B,eye(n),0);
[Y,tt,xSimu] = lsim(sys_xp,U,t,x0);
end