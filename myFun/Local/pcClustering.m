function [idx, varargout] = pcClustering(Y,varargin)
%% Description
%{
Divides data into 2 clusters, according to the sign of the (first) PCA
mode. For now, it works only for K=2!
The number of clusters was called K, because of the relation to the
K-Means.
%}
%% Inputs
% nargin
nin = length(varargin);

% Size
[rows, cols] = size(Y);

% Number of clusters K
if nin > 0
    K = varargin{1};
end

%% PCA
% PCA options
opts.method=1; opts.con=0; opts.mycon=[]; opts.k=1; opts.cs=false;

% Create PCA object
obj = pcaData(Y',opts);

%% Cluster indicator vector
% T1st Principal direction
u = obj.modes(:,1); 

% 1st Principal Component
v = u' * Y'; v=v(:);

% Initialize the indicator vector
h = zeros(rows, 1);

% Create index position vectors
i = (v>0); j = (v<=0);

% Get the cluster indicator vector
if numel(i) == 0 || numel(j) ==0 % if no subdivision took place
    h = ones(rows, 1); % there is only one cluster
else
    h(i) = 1; h(j) = 2;
end

%% Output
% Cluster index indicator vector
idx = h;

% Optional outputs: Clustered data
if nargout > 1
    imax = numel(unique(idx));
    clustY=cell(imax,1);
    for i=1:imax
        clustY{i} = Y(idx==i,:);
    end
    varargout{1} = clustY;
end

end



