function [decoded_data, varargout] = pca_decode(obj, modes, scores, scaling_factors, mean_column, varargin)

% From the PCA score space to the data space
decoded_data = modes * scores;

% Relavant dimension
N = size(modes, 1); 

% If a further input is passed, this function understands that a centroid
% has to be added back into the data (this method is usually called with
% one additional input when working with a local cluster).

if ~isempty(varargin)
    if obj.is_mesh_variable
        centroid = varargin{1}(:)';
        % The input mean_column will be interpreted as the centroid
        M = repmat(centroid, size(decoded_data, 1), 1);
        decoded_data = M + decoded_data;
    else
        centroid = varargin{1}(:);
        % The input mean_column will be interpreted as the centroid
        M = repmat(centroid, 1, size(decoded_data, 2));
        decoded_data = M + decoded_data;
    end
end

% Once the centroid has been added back, we obtain the data that was
% clustered (usually, centered-scaled data). We need to
% uncenter-unscale the data in order to build back the original training
% data. Anyway, it is chosen not to do this operation now.

if ~obj.is_local
    decoded_data = obj.uncenterUnscale(decoded_data, mean_column, scaling_factors);
end


end


