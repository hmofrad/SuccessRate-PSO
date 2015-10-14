% Schwefel’s Problem 1.2
function fit=f3(x) 
    [G,D] = size(x);
    fit = zeros(G,1);
    for i=1:G
        for j=1:D
            for k=1:j
                fit(i) = fit(i) + x(i,k)^2;
            end
        end
    end
end