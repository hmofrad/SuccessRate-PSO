% Calculate search direction from phi via
% a polar to Cartesian coordinate transform
% for one member of group using equation (1)
function [D] = Polar2Cartesian (phi)
    dim = length(phi) + 1;
    D = zeros (1,dim);       % D = [d1 ... dn];
    D(1) = prod(cos(phi),2); % calculate d1
    for i = 2:dim
        D(i) = sin(phi(i-1)) * prod(cos(phi(i:end))); % calculate di
    end
    D(end) = sin(phi(end));  % calculate dn
end
