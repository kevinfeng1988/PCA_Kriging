function [nz_X_k, nz_idx_clust, k, varargout] = partitionVQ(obj, X, idx, varargin)

    [rows, columns] = size(X); % Samples, Variables

    threshold = columns;
    
    isremove = true; 
    if ~isempty(varargin)
        isremove = varargin{1};
    end
    if ~isa(isremove,'logical')
        isremove = false;
    end
    
%     k = numel(unique(idx));
    k = max(idx);

    idx_clust = cell(k, 1);
    n_points = zeros(k, 1); 

    for j = 1 : k
        idx_clust{j} = find(idx == j);
        n_points(j) = size(idx_clust{j}, 1);
    %%% Uncomment -> default file
        if isremove
            if (n_points(j) < threshold)
                fprintf('\nNo points in the cluster n. %d, cluster removed \n', j);
            end
        end
    end
    if isremove ~= 1
        nz_idx = find(n_points > -1);
    else
        nz_idx = find(n_points > threshold);  % Uncomment -> default file
    end
    k_new = size(nz_idx, 1);
    k = k_new;
    nz_X_k = cell(k, 1);
    nz_idx_clust = cell(k, 1);
    for j = 1 : k
        nz_X_k{j} = zeros(n_points(j), columns);
        nz_idx_clust{j} = idx_clust{nz_idx(j)};
        nz_X_k{j} = X(nz_idx_clust{j}, :);
    end
    
    if nargout > 3
        idx = zeros(rows,1);
        for j = 1 : k
             idx(nz_idx_clust{j}) = j;
        end
        varargout{1} = idx;
    end
end



 
       