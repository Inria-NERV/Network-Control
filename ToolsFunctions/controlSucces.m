function s = controlSucces(x,xf,x0, eps)
if eps =='auto'
    eps = norm(xf - x0)/100;
end
s = (norm(x(end,:)' - xf) < eps);
end