function getKlpcaPredictions(obj, varargin)

% Specific cluster
idx_clust = [];
if ~isempty(varargin)
    idx_clust = varargin{1};
end

% Choose according to the clustered dimension of the data matrix
if ~obj.clustering_dimension
    % Work differently for input varargin{1} 
    if ~isempty(idx_clust)
        % If a specifi cluster is indicated, only work for that
        obj.local_pca{idx_clust}.getKpcaPredictions();
        I = obj.local_nz_idx_clust{idx_clust};
        obj.klpca_predictions(I,:) = obj.local_pca{idx_clust}.kpca_predictions;
    else
        % Rebuild from all clusters
        for ii = 1 : length(obj.local_pca)
            obj.local_pca{ii}.getKpcaPredictions();
            I = (obj.local_idx == ii);
            obj.klpca_predictions(I,:) = obj.local_pca{ii}.kpca_predictions;
        end
    end
else
    % Get useful sizes
    n_predictions = size(obj.prediction_points, 1);
    n_clusters = length(obj.local_pca);
    % Rebuild from all clusters
    dist = zeros(n_clusters, n_predictions);
    for ii = 1 : n_clusters
        obj.local_pca{ii}.getKpcaPredictions();
        % Get distance of each cluster from every prediction point
        d = getDistances(obj.local_pca{ii}.training_points, obj.prediction_points, obj.dist_on_cs_points);
        dist(ii,:) = d(:)';
    end
    % Rebuild for each prediction point
%     obj.klpca_predictions = combinePredictions(predictions, mse, kriged);
    dim_out = size(obj.local_pca{1}.kpca_predictions, 1);
    for jj = 1 : n_predictions
        % Weights evaluated from the meta-features
%         weights = 1 / dist(:,jj); % Get weights for each prediction
%         tot_w = sum(weights); % Sum of the weights
        % For each cluster, get the j-th prediction and weigh it
        % accordingly. The prediction for the j-th prediction point will be
        % a weighed average of the clusters predictions, with the distance
        % evaluated before as weights
        % For prediction jj, get one MSE per model/cluster
        mse_model = zeros(n_clusters, 1);
        for ii = 1 : n_clusters
            % MSE check
            temp = abs(obj.local_pca{ii}.mse(:,jj)) ./ (abs(obj.local_pca{ii}.kriged_pca_scores(:,jj)) + eps);
            mse_model(ii) = max( mean( temp, 1) ); % MSE for the ii-th model
        end
        % Weights inversely proportional to MSE
        weights = 1 ./ mse_model; 
        weights = weights ./ sum(weights); % Weights for the jj-th prediction
        sol = zeros(dim_out, n_clusters); % Columns are each model's prediction for jj-th prediction point 
        for ii = 1 : n_clusters
            sol(:,ii) = obj.local_pca{ii}.kpca_predictions(:,jj); % Load prediction for model kk
        end
        obj.klpca_predictions(:,jj) = sol * (weights); % Weighed average of the predictions
        % Use also the distances as meta-features
        weights_2 = 1 ./ dist(:,jj);
        weights_2 = weights_2 ./ sum(weights_2);
        obj.klpca_predictions(:,jj) = sol * (weights + weights_2) / 2;
        % We can also get the closest one, instead of combining
        [~, I] = max(weights_2);
        obj.klpca_predictions(:,jj) = sol(:,I);
    end
end

% They predicted the centered-scaled data
obj.klpca_predictions = obj.uncenterUnscale(... 
                        obj.klpca_predictions, ...
                        obj.mean_column, ...
                        obj.scaling_factors);
    
end



