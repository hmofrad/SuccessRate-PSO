% Generalized Rastrigin’s Function
function fit=f9(x)
    fit=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end