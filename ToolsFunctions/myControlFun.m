function [X_opt, xSimu, U , P_adj] = myControlFun(A, B, T,STEP, t, x0, xf, rho, S)
n =size(A,1);
nSamps = length(t);
sys_xp = ss(A,B,eye(n),0);
[X_opt, U_opt, n_err , P_adj] = optim_fun(A, T,STEP, B, x0, xf, rho, S);

U = reshape(U_opt , [n , nSamps])';

[Y,tt,xSimu] = lsim(sys_xp,U,t,x0);

[ xOpenLoopbis ] = open_loop_control( A, B, x0, U_opt , STEP);
%[ xSelfDynamicBis ] = open_loop_control( A, B, x0, zeros(size(U_opt)) , STEP);

%xfSimu = squeeze( xOpenLoopbis(:,1,end)) ;

end