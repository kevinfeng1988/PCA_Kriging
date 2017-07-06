function runClustering(obj, varargin)
%% Description
%{
    Routine that partitions data into clusters.
%}


%% Input
nin = length(varargin);
if nin > 0
    obj.number_of_clusters = varargin{1};
end

%% Check the number of clusters

if obj.number_of_clusters == 1
    obj.local_idx = ones(size(obj.training_data, obj.clustering_dimension + 1), 1);
    obj.local_pca = {};
    obj.local_nz_idx_clust = {};
    fprintf('\nLocal PCA was run, with only one cluster.\n');
    return; % Do not run the clustering optimization procedure
end

%% Clustering procedure
fprintf('\nClustering data... \n');

if obj.clustering_dimension && obj.lpca_clustering_options == 0
    obj.choosePartitioningCriteria(); % Partition variables
elseif ~obj.clustering_dimension && obj.lpca_clustering_options == 0
    obj.choosePartitioningCriteria(); % Partion observations
end

fprintf('\n...data clustered. \n');

obj.local_idx_stored{end+1} = obj.local_idx;
obj.number_of_clusters_stored(end+1) = numel(unique(obj.local_idx));

end




