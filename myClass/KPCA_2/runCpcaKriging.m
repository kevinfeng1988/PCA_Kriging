function runCpcaKriging(obj, varargin)

nin = length(varargin);
if nin > 0 && ~isempty(varargin{1})
    obj.kriged_pca_scores = varargin{1};
else
    % Kriging on the CPCA scores
    fprintf('Kriging on the CPCA scores is in progress...\n');
    temp = 'cpca';
    test_data = []; % Test data for CPCA is too costly to evaluate
    evalc('obj.update_Kriging(obj.cpca_scores, temp, test_data);');
    fprintf('Kriging on the CPCA scores has terminated.\n');
end

% Terminate here if this is a LPCA object
if obj.is_local
    return
end

% Get the predictions
obj.getKcpcaPredictions();

% Estimate the errors
obj.getKcpcaPredictionsErrors();

end

