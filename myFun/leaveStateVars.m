function Y = leaveStateVars(Q, xq, vars, p, varargin)
    
% Useful variables
if length(vars) == 1 
    nVars = vars; 
else
    nVars = length(vars);
end
if nVars < 1
    error('No variables are present!');
end

if size(p, 1) == 1 
    nSamples = p; 
else
    nSamples = size(p, 1); 
end

if size(xq, 1) == 1
    n = xq; 
else
    n = size(xq, 1); 
end
    
    
%% Main

% Initialize matrix
Y = zeros(nVars*n, nSamples);

% Leave "StateVars"
for j = 1:nSamples
    temp = [];
    
    for i = 1 : nVars
        j1 = 1+(j-1)*n;
        j2 = j*n;
        temp = [temp; Q(i, j1:j2)'];
    end
    
    % Data matrix
    Y(:,j) = temp;      % (n x nVars) X (nSamples)           
end
    
    
end




