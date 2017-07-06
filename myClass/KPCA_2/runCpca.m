function runCpca(obj, varargin)
%% Inputs
n_args = length(varargin);
if n_args > 0
    % Run CPCA with this approximation order
    if varargin{1} > 0
        obj.pca_approximation_order = varargin{1};
    end
end

%% CPCA
% Evaluate the CPCA scores
obj.update_gamma();

% Recover data
obj.recoverCpca();

% Estimate errors
obj.getCpcaErrors();

end





