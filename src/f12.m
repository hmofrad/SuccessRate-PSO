% Generalized Penalized Functions 1
function fit = f13(x)
    [G,D] = size(x);
    y  = 1 + 1/4 .* (x + 1);
    c1 = sin(3 * pi * y(:,1)).^2;
    c2 = sum(((y(:,1:D-1)-1).^2) .* (1 + 10 .* sin(pi * y(:,2:D)).^2),2);
    c3 = (y(:,D)-1).^2;
    fit = pi/D * (c1 + c2 + c3);

    u = zeros (1,D);
    for i = 1:G
        for j = 1:D
            u(j) = U(x(i,j));
        end
        fit(i) = fit(i) + sum(u);
    end
    % U function
    function u = U(x)
        a = 10;
        k = 100;
        m = 4;
        if (x>a)
            u = k*(x - a)^m;
        elseif (x >= -a && x <= a)
            u = 0;
        elseif  (x < -a)
            u = k*(-x - a)^m;
        end
    end
end