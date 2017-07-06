function [predictions, mse, k] = performKriging(samples, values, prediction_points, trendFun, corrFun, varargin)
% Starting guess for the evaluation of the hyperparameters for Kriging
guess = 1e-3; % Manually set only in case it is not passed as argument
if ~isempty(varargin)
    guess = varargin{1};
end
morevalues = [];
moresamples = [];
if length(varargin) > 2
    morevalues = varargin{2};
    moresamples = varargin{3};
end

% NOTE: Different Krigig models (for example with different h-parameters or
% trend/correlation functions) might be mixed in order to have a better
% predictor [NOT IMPLEMENTED].
% --- Pseudo-code 1:
% if(there_is_more_than_one_model)
% k_i = ood(samples, values, guess_i, trendFun_i, corrFun_i);
% [predictions_i, mse] = k_i.predict(prediction_points);
% predictions = combine_models(predictions_1, predictions_2, ...);

% Kriging, training the model
if isempty(morevalues) || isempty(moresamples)
    k = ood(samples, values, guess, trendFun, corrFun); 
else
    k = ood(moresamples, morevalues, guess, trendFun, corrFun);
end
if ~isempty(prediction_points)
    % Predictions
    [predictions, mse] = k.predict(prediction_points);
    % NaN check
    predictions = predictions';     mse = mse';
    checkNaN = isnan(predictions);  predictions(checkNaN) = 0;
    checkNaN = isnan(mse);          mse(checkNaN) = 0;
else
    predictions = [];
    mse = [];
end

end



function k = ood(samples, values, guess, varargin)
% Default functions
trendFun = 'regpoly2';   % Trend function
corrFun = @corrmatern32; % Correlation function

% Inputs
nin = length(varargin);
if nin > 0 
    trendFun = varargin{1}; % User-provided trend function
    if nin > 1
        corrFun = varargin{2};  % User-provided correlation function
    end
end

% Check validity of the inputs
if ~strcmp(trendFun(1:end-1),'regpoly') && ~strcmp(trendFun,'')
    fprintf('\tRegression function was mispecified, the default one was set.\t');
    trendFun = 'regpoly2'; % If the trend fun was not properply specified, choose the default one
end

% Samples
values = values'; % Transpose VALUES (to be interfaced with the Kriging toolbox)
[nSamples, dim] = size(samples);

% Check if the chosen trend function is ok 
if strcmp(trendFun, 'regpoly2')
    p = .5 * (dim + 1) * (dim + 2);
    if nSamples < p
        trendFun = 'regpoly1';
    end
end
if strcmp(trendFun, 'regpoly1')
    p = dim + 1;
    if nSamples < p
        trendFun = 'regpoly0';
    end
end

% Corrgaussp needs one more hyperparameter
if strcmp(corrFun, 'corrgaussp')
   dim = dim + 1; 
end

isKriging = true;
if isKriging
    % Kriging options
    opts = Kriging.getDefaultOptions();
    opts.hpOptimizer = SQPLabOptimizer( dim, 1 );
    opts.regressionMaxOrder = dim; 

    % Hyperparameters 
    Val = 1e10;                            
    lb = zeros(1, dim);  ub = Val * ones(1, dim); % (1 x d) Lower and upper bounds
    theta0 = guess * ones(1, dim);  % (1 x d) Guess
    opts.hpBounds = [lb; ub];                 % (2 x d) Hyperparameter optimization bounds

    % Kriging
    evalc('k = Kriging( opts, theta0, trendFun, corrFun );');
    try
        evalc('k = k.fit(samples, values);');
    catch ME 
        evalc('k = k.fit(samples, values);');
    end
else
    % Blind Kriging options
    opts = BlindKriging.getDefaultOptions();
    opts.hpOptimizer = SQPLabOptimizer( dim, 1 );
    opts.regressionMaxOrder = dim; 

    % Hyperparameters 
    Val = 1e10;                            
    lb = zeros(1, dim);  ub = Val * ones(1, dim); % (1 x d) Lower and upper bounds
    theta0 = guess * ones(1, dim);  % (1 x d) Guess
    opts.hpBounds = [lb; ub];                 % (2 x d) Hyperparameter optimization bounds

    % Blind Kriging
    evalc('k = BlindKriging( opts, theta0, trendFun, corrFun );');
    try
        evalc('k = k.fit(samples, values);');
    catch ME 
        evalc('k = k.fit(samples, values);');
    end
end

end


