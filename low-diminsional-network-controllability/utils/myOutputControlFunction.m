function [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] =...
    myOutputControlFunction(A, B, T, t, x0, xf, rho , R , Q , C , yf  )

% if nargin < 10
%     prec = 0;
% end

n =size(A,1);
% nSamps = length(t);
sys = ss(A,B,C,0);

%%%%%%%%%%%%%%%%%%%%%%%

Atilde = [A -B*(R\(B'))/(2*rho) ; -2*Q -A'];

M = expm(Atilde*T);
M11 = M(1:n,1:n);
M12 = M(1:n,n+1:end);
M21 = M(n+1:end,1:n);
M22 = M(n+1:end,n+1:end);

condAtilde = cond(Atilde);


if norm(Q) == 0
    c = zeros(2*n,1);
else
    N = Atilde\(M-eye(size(Atilde)));
    c = N*[zeros(n);Q]*2*xf;
end

% c =  zeros(2*n, 1);
c1 = c(1:n);
c2 = c(n+1:end);


touchyMat = M22 - 2*(C')*C*M12;
b = (+2 *(C')*C*M11 - M21)*x0 + (2 *(C')*C*c1 - c2) - 2 *(C')*yf   ;


p0 = pinv(touchyMat) * b;

n_err = norm(touchyMat*p0-b,1)/(norm(touchyMat,1)*norm(b,1));

condPinv = cond(touchyMat);


sys_xp = ss(Atilde,[zeros(n);zeros(n)],eye(2*n),0);

%%% output case
U_const = [];
while size(U_const,1) < length(t)
    U_const = [U_const;2*((C')*yf)'];
end


[Y,tt,xp] = lsim(sys_xp,U_const,t,[x0;p0]);


U = [];
for i = 1:length(t)
    U(i,:) = (-(1/(2*rho))*(R\(B'))*(xp(i,n+1:end)'))';
end

% X_opt = xp(:,1:n);
P_adj = xp(:,n + 1:end);


%%%%%%%%%%%%%%%%%%%%%%

[ySimu,tt,xSimu] = lsim(sys,U,t,x0);


end