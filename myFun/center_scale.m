function z = center_scale(X, mu, sigma)

% [] is a special case for std and mean, just handle it out here.
if isequal(X, [])
    z = []; 
    return; 
end

% Compute X's mean and sd, and standardize it
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus, X, mu);
z = bsxfun(@rdivide, z, sigma0);

end

