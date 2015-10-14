% Ackley’s Function
function fit = f10(x)
    [G,D] = size(x);
    fit = sum(x.^2,2);
    fit = -20 .* exp(-0.2 .* sqrt(fit ./ D)) - exp(sum(cos(2 .* pi .* x) , 2) ./ D) + 20 + exp(1);
end