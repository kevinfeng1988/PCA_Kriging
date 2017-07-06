function runPca_hiddenFunction(obj, varargin)

% Centered and scaled data (mean and scaling factors)
if ~obj.is_local 
    obj.runCenterScale();
% Do not run if this is a local cluster, the property centered_scaled_data
% has probably been already set.
end

% Run PCA
if obj.is_local && obj.is_mesh_variable
    [obj.pca_modes, scores, obj.pca_eigenvalues, obj.correlation_matrix] = local_fun_PCA(obj.training_data);
else
    [obj.pca_modes, scores, obj.pca_eigenvalues, obj.correlation_matrix] = local_fun_PCA(obj.centered_scaled_data);
end

% Scores
obj.pca_scores = scores';

% Set approximation order if not set, change it if greater than the size of
% the PCA modes matrix
if isempty(obj.pca_approximation_order)
    obj.pca_approximation_order = size(obj.pca_modes, 2);
elseif obj.pca_approximation_order > size(obj.pca_modes, 2)
    obj.pca_approximation_order = size(obj.pca_modes, 2);
end

end


function [coeff, scores, latent, C] = local_fun_PCA(data, varargin)
%% Dimensions
[n_var, n_samples] = size(data);

%% Inputs
n_args = length(varargin);
% Center-scale?
cen = false; 
if n_args > 0
    cen = varargin{1};
end
% Choose algorithm
if n_var > n_samples
    alg = 'svd';
else
    alg = 'svd';
end
if n_args > 1
    alg = varargin{2};
end

%% Correlation matrix
try
    if n_var > n_samples
        C = cov(data);
    else
        C = cov(data');
    end
catch ME
    C = 0;
end

%% Apply PCA (Using Matlab's default toolbox)
[coeff, scores, latent] = pca(data', 'Centered',cen, 'Algorithm',alg);

end





