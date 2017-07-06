function getKlcpcaPredictions(obj, varargin)

% Specific cluster
idx_clust = [];
if ~isempty(varargin)
    idx_clust = varargin{1};
end

% Choose according the clustered dimension of the data matrix
if ~obj.clustering_dimension
    % Work differently for input varargin{1} 
    if ~isempty(idx_clust)
        % If a specifi cluster is indicated, only work for that
        obj.local_pca{idx_clust}.getKcpcaPredictions();
        I = obj.local_nz_idx_clust{idx_clust};
        obj.klcpca_predictions(I,:) = obj.local_pca{idx_clust}.kcpca_predictions;
    else
        % Rebuild from all clusters
        for ii = 1 : length(obj.local_pca)
            obj.local_pca{ii}.getKcpcaPredictions();
            I = obj.local_nz_idx_clust{ii};
            obj.klcpca_predictions(I,:) = obj.local_pca{ii}.kcpca_predictions;
        end
    end
else
    % There already should be the correct sets of Local PCA scores
    % from which to decode data
    
end
    
end