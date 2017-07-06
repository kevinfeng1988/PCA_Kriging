function [Infl_rel, varargout] = samplingInfluence(points, Y, varargin)
%% Description
%{
part of: "Adaptive Sampling for basis improvements"

Evaluation of the Influence of each snapshot: 
- on the PCA basis;
- on the regression model;

INPUTS
    Y: matrix of observations (column-vectors)
    points: input space locations of the observations

OUTPUTS
    Infl_rel: relative influence of each observation
%}


%% Input check
% Optional inputs
n_args = length(varargin);
reg_based = false; % Regression-based influence
if n_args > 0
    reg_based = varargin{1};
end

% Number of points
n_points = size(points,1);

%% PCA 
pca_model = pcaFun(Y); % Create a PCA structure

%% Snapshot influences on the modes
% Memory allocation
pod_nosnap = cell(n_points,1);
Infl = [];

% Create the matrix of snapshot influences on the modes
for j = 1 : n_points
    % Create the matrix with the missing snapshot
    Yj = Y;   
    Yj(:,j) = []; 
    % PCA with a missing snapshot
    pod_nosnap{j} = pcaFun(Yj); 
    % Infl(i,j): influence of snapshot j on the pca mode i
    for i = 1 : pod_nosnap{j}.k
        modeprod = abs( pca_model.modes(:,i)' * pod_nosnap{j}.modes(:,i) );
        Infl(i,j) = 1/(modeprod+eps) - 1; 
    end
end

%% Influence of the snapshots on the modal basis
% Memory allocation
Infl_s = zeros(length(pod_nosnap), 1);

% Influence of the snapshot j on the modal basis
for j = 1 : length(pod_nosnap)
    sv = (pod_nosnap{j}.eigenv).^(.5); sv = sv(:); % Singular value
    Infl_s(j,1) = sv' * Infl(:,j);
end

%% Output
% Relative influence of the jth snapshot on the modal basis
Infl_rel = Infl_s / norm(Infl_s, 1);

% Return output
if nargout > 1 
    varargout{1} = Infl;
end


end




function this = pcaFun(Y)


[this.Ycs, this.m, this.d] = zscore(Y, 0, 2);


% Apply PCA (Using Matlab's default toolbox)
[coeff, scores, latent] = pca(this.Ycs');


% Approximation order
this.k = size(coeff, 2);

% Modes
this.modes = coeff;        

% Scores
this.a = scores';           

% Eigenvalues
this.eigenv = latent;       

end



