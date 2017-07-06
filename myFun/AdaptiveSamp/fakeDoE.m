function [p, Yas, p_miss, Pot, Infl_rel, varargout] = fakeDoE(p_orig, Y_orig, p0, nSamplesMax, varargin)
%% Description
%{
This function performs a fake Design of Experiment.
Given an 'original' dataset Y_orig, of M+W observations, this function
chooses the best M osbervations, with M provided by the user.

The term 'fake' is due to the fact that, once Adaptive Sampling has chosen
the new point to sample, the observation is added to the dataset, but this
observation is not the output of any simulation: we already have
the data!

This can be useful when performing validation on predicting models: we have
a dataset, we select the training observations (with the Adaptive Sampling
technique) and use the remainder for validation of the predictions.


INPUTS

p_orig:
Y_orig:
p0:
nSamplesMax:


OUTPUTS

p:
Yas:
p_miss;
Infle_rel:

%}



%% Inputs

% Size
[rows, cols] = size(Y_orig);

% Sorting
[r, c] = size(p0); % c: input space dimension; r: number of points
if c == 1
    p_orig = sort(p_orig);
    p0 = sort(p0);
end

% At least 3 training points have to be provided
if r < 3
    error('Not enough starting points.');
end

% Final number of samples (check)
if nSamplesMax > cols
    error('Too many samples demanded. nSampleMax > number of cols of Y_orig.');
end

% Consistency: p_orig and Y_orig
if cols ~= length(p_orig)
    error('Number of original points not equal to number of columns of Y_orig.');
end

% Membership of p0 in p_orig
if any(~ismember(p0,p_orig))
    error('Asked to sample a point for which there is no data.');
end


% Optional inputs
nin = length(varargin);

% Filter for CDP
cdp_filter = false;
if nin > 0
    cdp_filter = varargin{1};
end



%% Fake DoE

% Candidate points
cdp = p_orig;             
i = ismember(cdp, p0, 'rows');     
cdp(i,:) = [];

% The initial data set
Yas = sampleit(p0, p_orig, Y_orig);

Pot = [];
% Enrich the data set
while size(p0, 1) < nSamplesMax && ~isempty(cdp)

    % PCA Approximation Order
    opts.k = size(p0,1) - 1;
    if length(p0) < 3
        opts.k = 1;
    end

    % Get a new point
    tic;
    [newpoint, pot, ~] = AdaptiveSampling(p0, Yas, cdp);
    time_as = toc; fprintf('\nAS: %.2f s. ', time_as);

    % Update the training points
    p0 = [p0; newpoint]; p0 = sortrows(p0); 

    % Save Pot and Infl
    Pot = [Pot; max(pot)]; 

    % Update the candidate points
    i = ismember(cdp, p0, 'rows'); 
    cdp(i,:) = [];

    % Update the data set
    Yas = sampleit(p0, p_orig, Y_orig);
    
    % Faster version
    if cdp_filter
        n_max = floor(.2 * (size(p_orig,1) - size(p0,1)));
        cdp = discrete_lhs(p_orig, p0, n_max);
    end
    
    % Print status
    prog = 100 * (length(p0) / nSamplesMax);
    fprintf('\nProgress is at %.2f per cent.\n', prog);
end

% Training points chosen
p = sortrows(p0); 

% Influences
Infl_rel = samplingInfluence(p, Yas, opts);

% Get the prediction points
I = ~ismember(p_orig, p, 'rows');
p_miss = p_orig(I,:);



%% Outputs
% p, Yas, p_miss, Pot, Infl_rel



%% Varargout



end




