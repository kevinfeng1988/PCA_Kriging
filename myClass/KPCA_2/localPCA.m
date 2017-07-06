function [idx, uncentered_nz_X_k, varargout] = localPCA(obj, scal_X, n_eigs, k, varargin)

nin = length(varargin);


% Dimensions:
[rows, columns] = size(scal_X);

% The input argument X that is passed is already the centered-scaled data,
% see the method choosePartitioningCriteria().

%%%%%%%%%%% Convergence indicators initialization %%%%%%%%%%%
convergence = 0;    
iter = 0;                   
iter_max = 1e4;
eps_rec = 1.0;      
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;    
r_tol = 1.0e-08;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1) INITIALIZATION
    
% If something went wrong, use the old code (put aside)
% if false && (length(unique(idx)) ~= k)
%     % Initialize centroids
%     C_int = linspace(1, rows, k+2);     C = scal_X(round(C_int(2:k+1)), :);  
%     
%     % Initialize clusters 
%     nz_X_k = cell(k,1); 
%     sp = floor(rows/k);
%     for j = 1 : k - 1
%         nz_X_k{j} = scal_X(1 + sp*(j-1):j*sp, :);
%     end
%     nz_X_k{k} = scal_X(1 + sp*(k-1):end, :);
% end

% Centroids initialization
if nin > 0 && ~isempty(varargin{1})
    % User-supplied centroids initialization
    idx_0 = varargin{1}; 
    % Get clusters
    nz_X_k = cell(k, 1); 
    for j = 1 : k
        nz_X_k{j} = scal_X(idx_0 == j, :);
    end
    % Get cluster centroids
    C = zeros(k, columns);        
    for j = 1 : k
        C(j, :) = mean(nz_X_k{j}, 1);
    end
else
    % KMEANS initialization
    fprintf('\nCentroids and cluster assignments initialization. \n');
    [idx, C] = kmeans(scal_X, k, 'MaxIter', 1e8); % Use kmeans to initialize 
    k = length(unique(idx)); % Be sure k clusters were found 
    % Initialize clusters
    nz_X_k = cell(k,1); 
    for j = 1 : k
        nz_X_k{j} = scal_X(idx == j, :);
    end
end

% Initialization of eigenvectors
[eigvec, ~, ~] = performlocalpca(nz_X_k, n_eigs, C);
% eigvec = cell(k,1); 
% n_eigs_max = n_eigs;
% parfor j = 1 : k
%     eigvec{j} = eye(columns, n_eigs);
% end


%% 2) PARTITION
eps_rec_var = Inf;

while (convergence == 0 && iter < iter_max) && (k ~= 1)
    
    C_convergence = 0;      eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);
    date
    datestr(now, 'HH:MM:SS')

    sq_rec_err = zeros(rows, k);

% Choose the metric
    parfor j = 1 : k
        C_mat = repmat(C(j, :), rows, 1); 
        
        % Squared mean reconstruction error
        rec_err_os = scal_X - C_mat - (scal_X - C_mat) * eigvec{j} * eigvec{j}';
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2); 
    end
    
    % Evalutate the recovered minimum error
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    
    % Partition the data into clusters
    if rows > columns
        isremove = true;
    else
        isremove = false;
    end
    isremove = false; % ???
    [nz_X_k, nz_idx_clust, k] = obj.partitionVQ(scal_X, idx, isremove);

    fprintf('\nThe current number of clusters is %d.\n', length(nz_X_k));
    
    % Evaluate the relative recontruction errors in each cluster
    rec_err_min_rel_k = cell(k, 1);
    for j = 1 : k
        rec_err_min_rel_k{j} = rec_err_min_rel(nz_idx_clust{j}, 1);
    end
    
    % Evaluate the mean error in each cluster
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);
    for j = 1 : k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_X_k{j}, 1);
    end
    fprintf(...
        '\nGlobal mean recontruction error at iteration n. %d equal to %d',...
        iter, eps_rec_new);

        

%% 3) EVALUATE NEW CLUSTERS' CENTROIDS
    C_new = zeros(k, columns);        
    for j = 1 : k
        C_new(j, :) = mean(nz_X_k{j}, 1);
    end
    eps_old = eps_rec_var;
    eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
    fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);
%     if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) ...
%             && (n_eigs < n_eigs_max)) 
%         n_eigs = n_eigs + 1;
%         fprintf('\n Clusters dimension increased to %d \n', n_eigs);
%     end
    
    % Save the best idx yet
    if eps_old < eps_rec_var
        bestidx = idx;
    end

    % Judge convergence: clusters centroid and relative reconstruction
    % error
    if (eps_rec_var < r_tol)
        eps_rec_convergence = 1;
    end
    if (size(C) == size(C_new))
        C_var = abs((C_new - C) ./ (C_new + a_tol));
        if (C_var(:, :) < r_tol)
            C_convergence = 1;
        end
    end
    if ((iter > 1) && (C_convergence == 1) && (eps_rec_convergence == 1))
        convergence = 1;
        fprintf('\nConvergence reached in %d iteartion \n', iter);
    end

    % Update recontruction error and cluster centroids
    C = C_new;
    eps_rec = eps_rec_new;


%% 4) PERFORM LOCAL PCA
    [eigvec, ~, ~] = performlocalpca(nz_X_k, n_eigs, C);
    
    
    iter = iter + 1;
end


%% OUTPUT

% Return the clusters
if k == 1
    idx = ones(rows,1);
    [nz_X_k, nz_idx_clust, k] = obj.partitionVQ(scal_X, idx);
else
    [nz_X_k, nz_idx_clust, k] = obj.partitionVQ(scal_X, idx, isremove);
end


% Message about convergence
if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end


% Return the data
uncentered_nz_X_k = cell(k, 1);
for j = 1 : k
    uncentered_nz_X_k{j} = scal_X(nz_idx_clust{j}, :);
end

fprintf('\n\n');

% Optional outputs
if nargout > 2
    varargout{1}.nz_X_k = nz_X_k;
    varargout{1}.k = k;
    varargout{1}.nz_idx_clust = nz_idx_clust;
end

close all hidden


end




% Centroidal-Voronoi-Tessellations-based PCA
function delta = cvtpca(firstIn, varargin)
%{
This function is useful when clustering the data (before performing PCA).
It gives a metric to evaluate the distance between a space-vector and a
cluster's centroid.
%}


% INPUTS CHECK:
h = whos('firstIn');
if strcmp(h.class,'pcaData') || strcmp(h.class,'podData')
    X = firstIn.Y;      
    eigvec = firstIn.modes;
    gamma = firstIn.d;
else
    X = firstIn;
    if nargin < 2
        error('You must supply the second input: eigenvectors matrix.\n');
    end
    eigvec = varargin{1};
    if nargin > 2
        gamma = varargin{2};
    else
        gamma = ones(size(X, 1), 1);
    end
end
gamma = gamma(:);


% PRE-PROCESSING:
[n, m] = size(X);   
k = size(eigvec, 2);
delta = zeros(n, 1);

if max(size(gamma)) == 1
    eigvec_scaled = eigvec * gamma;
else
    eigvec_scaled = eigvec .* repmat(gamma, 1, k);
end


% MAIN:
tol = 1e-22;
for i = 1 : n
    x = X(i, :); 
    delta(i, 1) = 1 - ( 1 / (norm(x)^2 + tol) ) * ( norm(x * eigvec_scaled)^2 );
end
    
end





% Perform Local PCA
function [eigvec, n_eig, gamma, varargout] = performlocalpca(nz_X_k, n_eigs, C)

    % Number of clusters
    k = max(size(nz_X_k));
    
    % Initialization of cell arrays
    eigvec = cell(k, 1);
    n_eig = cell(k, 1);
    gamma = cell(k, 1);
    pod = cell(k,1);
    

    % Apply PCA (Using Matlab's default toolbox)
    parfor j = 1 : k
        % Load each cluster
        thisx = (nz_X_k{j} - repmat(C(j,:), size(nz_X_k{j},1), 1))';
        
        % Apply PCA
        [coeff, ~, ~] = pca(thisx'); 
        
        if n_eigs > size(coeff,2)
            this_n_eigs = size(coeff, 2);
        else
            this_n_eigs = n_eigs;
        end

        % Outputs
        n_eig{j} = n_eigs;
        eigvec{j} = coeff(:, 1:this_n_eigs);             
        gamma{j} = 1;
    end
          

    % Additional Outputs
    if nargout > 3
        varargout{1} = pod;
    end
    
end



