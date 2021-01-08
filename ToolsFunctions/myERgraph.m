function A = myERgraph(n , p , selfLoops)
A = zeros(n);
if selfLoops == 0
    for i=1:n
        for j=1:n
            if  i ~= j
                x = rand();
                if x < p
                    A(i,j) = 1;
                end
            end
        end
    end
elseif selfLoops == 1
    for i=1:n
        for j=1:n
            x = rand();
            if x < p
                A(i,j) = 1;
            end
            
        end
    end
end   
end