function normA = myNormalization(A,mode , c, n)
if strcmp(mode ,'basic')
    lambdaO = eig(A);
    normA = A /(max(abs(lambdaO)) + c) - eye(n);
elseif strcmp( mode ,'laplacian')
    Gtemp = digraph(A');
    deg = indegree(Gtemp)+outdegree(Gtemp);
    lambdaO = eig(A);
    normA = A /(max(abs(lambdaO)) + c) - diag(deg);
end
end