% f5 = Rosenbrock function
function fit = f5(x)
    [G,D] = size(x);
    fit = sum(100 .* (x(:,1:D-1) .^2 - x(:,2:D)) .^ 2  + (x(:,1:D-1)-1) .^2 , 2);
end