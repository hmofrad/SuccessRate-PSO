% Generalized Griewank Function
function fit=f11(x)
    [G,D]=size(x);
    fit=1;
    for i=1:D
        fit=fit.*cos(x(:,i)./sqrt(i));
    end
    fit=sum(x.^2,2)./4000-fit+1;
end