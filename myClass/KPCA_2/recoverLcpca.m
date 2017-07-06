function recoverLcpca(obj, varargin)

% Recover data from Local Cpca
for ii = 1 : length(obj.local_pca)
    % Recover data in one cluster
    obj.local_pca{ii}.recoverCpca(); 
    
    % Write the recovered info in the correct rows
    obj.local_recovered_data_cpca(obj.local_idx == ii) = ...
        obj.local_pca{ii}.cpca_recovered_data; 
end

end
