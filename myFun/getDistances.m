function [dist, varargout] = getDistances(training_points, prediction_points, varargin)

% Further inputs
n_args = length(varargin);
if n_args > 0
    cs_inputs = varargin{1};
    % Center-scale points
    if cs_inputs
        training_points = zscore(training_points, 0, 1);
        prediction_points = zscore(prediction_points, 0, 1);
    end
end

% Get distances of this set of training points to the j-th
% prediction point
n_samples = size(training_points, 1);
n_predictions = size(prediction_points, 1);
dist_t = zeros(n_samples, n_predictions);
for k = 1 : n_samples
    for jj = 1 : n_predictions
        dist_t(k,jj) = sqrt( norm( training_points(k,:) - prediction_points(jj,:) ) );
    end
end
% Get distances of this cluster from the prediction points, they
% will be used as meta-features
dist = zeros(n_predictions,1);
for jj = 1 : n_predictions
    dist(jj) = min(dist_t(:,jj));
end


end