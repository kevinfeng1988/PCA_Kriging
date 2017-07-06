function runLcpcaKriging(obj, varargin)

nin = length(varargin);
if nin > 0
    obj.pca_approximation_order = varargin{1};
end

fprintf('Kriging-LCPCA will be performed.\n');

% Kriging on the LPCA scores
for ii = 1 : length(obj.local_pca)
    obj.local_pca{ii}.runCpcaKriging();
    
    % Get the predictions
    I = obj.local_nz_idx_clust{ii};
    obj.klcpca_predictions(I,:) = obj.local_pca{ii}.kcpca_predictions;
end

% Estimate the errors
obj.getKlcpcaPredictionsErrors();

end
