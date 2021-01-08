function W = myGram(A , B , t , dt)
    n = size(A,1);
    W = zeros(n,n);
    for tau = t
        W = W + dt*(exp(tau*A)*(B*B')*exp(tau*A'));
    end
end