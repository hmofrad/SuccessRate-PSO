% Generalized Schwefel’s Problem 2.26
function fit=f8(x)
    [G D] = size(x);
    fit=420.9678-sum(x.*sin(sqrt(abs(x))),2);
    fit = 420.9678*D+fit;
end