function [point, varargout] = AdaptiveSampling(starting_points, Y, candidate_points, varargin)
%% Description
%{
Adaptive Sampling for basis improvements

INPUTS
    Y_full: data-set
    p: candidate points
    p0: starting points

OUTPUTS
    point: where to perform PCA

%}

%% Input check
% varargin size
nin = length(varargin);

% varargin check
opts = [];
if nin > 0
    opts = varargin{1};
end

% Rename input variables
p0 = starting_points;
p = candidate_points;

% Be sure that no value in p0 is in p
I = ismember(p,p0,'rows');
p(I) = [];


%% Relative Influence
Infl_rel = samplingInfluence(starting_points, Y);


%% Enrichment Potential
% Memory allocation
Pot = zeros(size(p,1),1);

% Potential of each candidate point
for i = 1 : size(Pot,1)
    
    % Evaluate the potential w.r.t. every starting point
    for j = 1 : size(p0,1)
        dist(j) = norm(p(i,:) - p0(j,:));
    end
    
    % Find the j for which temp is min
    [~, j] = min(dist);
    
    % Choose the j for which temp is max
    Pot(i) = dist(j) * Infl_rel(j);    
end


%% Candidate point to be chosen
% Get the index of the candidate point that maximizes the Pot
[~, ind] = max(Pot);


%% Output
% The candidate point
point = p(ind,:);

% Further outputs
if nargout > 1
    varargout{1} = Pot;
    if nargout > 2
        varargout{2} = Infl_rel;
    end
end


end
