function runLpcaKriging(obj, varargin)
%% Inputs

% Approximation order
nin = length(varargin);
if nin > 0 && ~isempty(varargin{1})
    obj.pca_approximation_order = varargin{1};
end

% Manually set the kriged pca scores
manset_scores = [];
if nin > 1 && ~isempty(varargin{2})
    manset_scores = varargin{2};
    if ~isa(manset_scores, 'cell')
        manset_scores = [];
    end
    if islogical(varargin{2}) && varargin{2}
        obj.get_Kriging_targets();
        manset_scores = true;
    end
end

fprintf('Kriging-LPCA will be performed.\n');

%% Kriging on the LPCA scores
for ii = 1 : length(obj.local_pca)
    % Set training points: depending what was clustered (variables or
    % observations), the property 'training_points' is set accordingly.
    if ~obj.clustering_dimension 
        % Clustered variable: same regression for all clusters
        obj.local_pca{ii}.training_points = obj.training_points; 
    else
        % Clustered observations: different regression per cluster
        I = ismember(obj.local_idx, ii);
        obj.local_pca{ii}.training_points = obj.training_points(I,:);
    end
    % In both cases, all clusters will try to make the same predictions
    obj.local_pca{ii}.prediction_points = obj.prediction_points;
    % Stacking?
    obj.local_pca{ii}.stacking = obj.stacking;
    obj.local_pca{ii}.trend_function = obj.trend_function;
    obj.local_pca{ii}.correlation_function = obj.correlation_function;
    obj.local_pca{ii}.includeMoreSamples = obj.includeMoreSamples;
    obj.local_pca{ii}.correlation_function_optimization = obj.correlation_function_optimization;
    
    % Be sure the LPCA object have info about the original data (WHY?)
%     if ~obj.clustering_dimension
%         obj.local_pca{ii}.original_data = ...
%                         obj.original_data(obj.local_nz_idx_clust{ii}, :);
%     else
%         obj.local_pca{ii}.original_data = ...
%                         obj.original_data(:,obj.local_nz_idx_clust{ii});
%     end
    
    % Run Kriging for the PCA scores inside the LPCA objects. The routine
    % that is run depends on two main events: (1) the kriged scores are
    % manually set; (2) clustering was performed on observations, not on
    % the variables. In the case (2), it is necessary to work differently.
    % Please remember that the property 'prediction_points' is the same for
    % every cluster, regardless of 'clustering_dimension'
    if isempty(manset_scores) && ~obj.clustering_dimension
        evalc('obj.local_pca{ii}.runPcaKriging();');
    elseif isempty(manset_scores) && obj.clustering_dimension 
        evalc('obj.local_pca{ii}.runPcaKriging();');
        % [TO BE CONTINUED]
        % Idea: perform stacking with the different Kriging models and/or
        % predictions.
        % For now: instead of stacking, weigh the predictions (of the
        % scores) based on the distance between prediction point and
        % average (or cluster's closest) training point.
    elseif islogical(manset_scores) && manset_scores
        evalc('obj.local_pca{ii}.runPcaKriging(obj.local_pca{ii}.pca_scores_kriging_targets);');
    elseif isa(manset_scores, 'cell')
        evalc('obj.local_pca{ii}.runPcaKriging(manset_scores{ii});');
    end
end

% Get predictions of LPCA
obj.getKlpcaPredictions();

fprintf('Kriging-LPCA was performed.\n');

% Estimate the errors
obj.getKlpcaPredictionsErrors();

end
