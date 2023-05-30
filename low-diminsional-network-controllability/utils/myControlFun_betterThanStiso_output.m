function [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] =...
    myControlFun_betterThanStiso_output(A, B, T,STEP, t, x0, xf, rho, S , R , Q , C , yf  )

% if nargin < 10
%     prec = 0;
% end

n =size(A,1);
% nSamps = length(t);
sys = ss(A,B,C,0);

%%%%%%%%%%%%%%%%%%%%%%%
invR = inv(R);
% Sbar = eye(n) - S;

Atilde = [A -B*invR*B'/(2*rho) ; -2*Q -A'];

M = expm(Atilde*T);
M11 = M(1:n,1:n);
M12 = M(1:n,n+1:end);
M21 = M(n+1:end,1:n);
M22 = M(n+1:end,n+1:end);

condAtilde = cond(Atilde);

% N = inv(Atilde) * (M-eye(size(Atilde)));
% c = N*[zeros(n);Q]*2*xf;
if norm(Q) == 0
    c = zeros(2*n,1);
else
    N = Atilde\(M-eye(size(Atilde)));
    c = N*[zeros(n);Q]*2*xf;
end

% c =  zeros(2*n, 1);
c1 = c(1:n);
c2 = c(n+1:end);

% M12 = M12;

%%% lagrange param
% outputSize = length(yf);
% v = (C')* lagrngParam * ones(outputSize , 1);
% A = M22 - 2*(C')*C*M12;
% b = (2 *(C')*C*M11 - M21)*x0 + (2 *(C')*C*c1 - c2) - 2 *(C')*yf + v  ;


touchyMat = M22 - 2*(C')*C*M12;
b = (+2 *(C')*C*M11 - M21)*x0 + (2 *(C')*C*c1 - c2) - 2 *(C')*yf   ;


p0 = pinv(touchyMat) * b;

n_err = norm(touchyMat*p0-b,1)/(norm(touchyMat,1)*norm(b,1));
% fprintf('relative error = %.3e \n' , relativeError);
% p0 = vpa(pinv(touchyMatrix)) * (-[S*M11;Sbar*M21]*x0 - [S*c1;Sbar*c2] + [S*xf;zeros(n,1)]);

% digits(32);
condPinv = cond(touchyMat);

% n_err = norm([S*M12;Sbar*M22]*p0 - (-[S*M11;Sbar*M21]*x0 - [S*c1;Sbar*c2] + [S*xf;zeros(n,1)])); % norm(error)

sys_xp = ss(Atilde,[zeros(n);2*S],eye(2*n),0);

% t = 0:STEP:T;

% U_const = [];
% while size(U_const,1) < length(t)
%     U_const = [U_const;2*xf'];
% end

%%% output case
U_const = [];
while size(U_const,1) < length(t)
    U_const = [U_const;2*((C')*yf)'];
end


[Y,tt,xp] = lsim(sys_xp,U_const,t,[x0;p0]);

% sys = ss(A,B*B'/(2*rho),eye(n),0);
% [Y,T,X] = lsim(sys,-xp(:,n+1:end),tt,x0);

U = [];
for i = 1:length(t)
    U(i,:) = (-(1/(2*rho))*invR*(B')*(xp(i,n+1:end)'))';
end
% Modif Here
% U_opt = reshape(U', [n,1,length(t)]);

% X_opt = xp(:,1:n);
P_adj = xp(:,n + 1:end);




%%%%%%%%%%%%%%%%%%%%%%

% U = reshape(U_opt , [n , nSamps])';

[ySimu,tt,xSimu] = lsim(sys,U,t,x0);

% [ xOpenLoopbis ] = open_loop_control( A, B, x0, U_opt , STEP);
%[ xSelfDynamicBis ] = open_loop_control( A, B, x0, zeros(size(U_opt)) , STEP);

%xfSimu = squeeze( xOpenLoopbis(:,1,end)) ;

end