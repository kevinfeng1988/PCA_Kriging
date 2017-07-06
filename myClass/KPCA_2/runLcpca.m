function runLcpca(obj, varargin)
% Perform LCPCA
for ii = 1 : length(obj.local_pca)
    obj.local_pca{ii}.my_constraint = obj.my_constraint;
    obj.local_pca{ii}.runCpca();
end

% Recover data
obj.recoverLcpca();

% Estimate errors
obj.getLcpcaErrors();


end



