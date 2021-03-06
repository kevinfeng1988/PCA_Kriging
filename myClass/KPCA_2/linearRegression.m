function [predictions, mse, k] = linearRegression(samples, values, prediction_points, varargin)
%% INPUTS:
%    x - independent variables.  These should be ordered with each
%        observation in a separate row and each variable in a column.
%    y - dependent variables.  These should be ordered with each observation
%        in a separate row and each variable in a column.
%   OPTIONAL INPUTS:
%    atol    - refinement tolerance
%    maxOrd  - maximum order for basis functions (default is 3)
%    nxtrial - number of (uniform) bins in x space to use (default is 4)
%    maxiter - maximum number of refinement iterations (default is 30)

tol = 1e-5;
maxOrd = 3;
nxtrial = 4;
maxiter = 50;

%% MARS
predictions = zeros( size(values,1), size(samples,1) );
mse = predictions * 0;
k = cell(size(values, 1),1);
parfor i = 1 : size(values, 1)
    k{i} = mars(samples, values(i,:)', tol, maxOrd, nxtrial, maxiter);
%     temp = k(prediction_points)';
    predictions(i,:) = k{i}(prediction_points)';
end


end


