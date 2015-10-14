% G. Quartic Function i.e. Noise
function fit=f7(x) 
    [G,D] = size(x);
    fit = sum(x.^4,2) + rand(G,1);
end