function [idx, uncentered_nz_X_k, varargout] = localPCA(X, n_eigs, varargin)
%% Input Check:
% 1:metric; 2:pcascorespace; 3:clnum
nin = length(varargin);

%% Dimensions:
[rows, columns] = size(X);
n_eigs_max = columns;

%% Names of the metric:
metricNames{1} = 'Euclidean'; 
metricNames{2} = 'PCA-based Reconstruction Error';
metricNames{3} = 'Not ready'; 
metricNames{4} = 'CVT-PCA'; % Centroidal Voronoi tessellations based PCA

%% CENTER AND SCALE THE DATA
X_old = X;
alongdim = 2;     
[scal_X, X_ave, X_gamma] = zscore(X, 0, alongdim);
scal_Xold = zscore(X_old, 0, alongdim);
if (alongdim==1)
    X_ave = repmat(X_ave, rows, 1);
elseif (alongdim==2)
    X_ave = repmat(X_ave, 1, columns);
end

%% Do you want to perform the partition on the PCA-space?
question = '\nPerform partitioning on the PCA scores space? [y/n]   ';
i = [];
if nin>1
    i = varargin{2};
end
if isempty(i)
    i = input(question);
end
if (i==1) || (i==true) || any(i==['y','Y'])
    ispcaspace = true; 
    opts.mycon=[]; opts.con=false; opts.method=1; 
    opts.k=3;%round(size(X,2));   
    temp=pcaData(X',opts);
    X = temp.a';
else
    ispcaspace = false;
end

%% Inputs check (the criterion can be chosen by providing it as input):
if ispcaspace
    metric=1;
else
    if (nin>0) && ~isempty(varargin{1})
        metric = varargin{1}; %get it from the input, if it exists
    else
        s='\nChoose the Partition criterion: \n';
        for i=1:length(metricNames)
            s = [s,' ',num2str(i),' - ',metricNames{i},'; \n'];
        end
        metric = input(s); %or ask the user to provide it
    end
end
critName = metricNames{metric};
fprintf('\nPartition criterion: %s\n',critName);

% F=X(:,1); % needed

%% INPUTS
% Choose the number od clusters
if nin>2 && ~isempty(varargin{3})
    k = varargin{3};
else
    k = input('\nSpecify the number of clusters \n');
end
% Convergence indicators initialization
convergence = 0;    iter = 0;                   
iter_max = 1200;
eps_rec = 1.0;      eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;    r_tol = 1.0e-05;

%% 1) INITIALIZATION
% Centroids initialization
C_int = linspace(1, rows, k+2);     C = scal_X(round(C_int(2:k+1)), :);

nz_X_k = cell(k,1); 
sp=floor(rows/k);
for j=1:k-1
    nz_X_k{j} = scal_X(1+sp*(j-1):j*sp,:);
end
nz_X_k{k} = scal_X(1+sp*(k-1):end,:);

if (strcmp(critName,metricNames{3}))
    C = zeros(k, columns);        
    for j=1:k
        C(j, :) = mean(nz_X_k{j}, 1);
    end
end

% Initialization of eigenvectors
[eigvec,n_eig,gamma] = performlocalpca(nz_X_k, n_eigs);

% eigvec=cell(k,1); gamma=cell(k,1);
% for j=1:k
%     eigvec{j}=eye(columns);
%     gamma{j}=1;
% end

%% 2) PARTITION

eps_rec_var = Inf;

while ((convergence == 0) && (iter < iter_max)) && (k~=1)
    
    C_convergence = 0;      eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);  

    sq_rec_err = zeros(rows, k);

%% Metric
    for j = 1 : k
        C_mat = repmat(C(j, :), rows, 1); 
        D = diag(gamma{j}); %D2=sparse(diag(gamma2{j}));
        
        if (metric==1)
            % Euclidean distance
            rec_err_os = (scal_X - C_mat);
            sq_rec_err(:, j) = sum(rec_err_os.^2, 2); 
            
        elseif (metric==2)
            % Squared mean reconstruction error
            PI = ( eigvec{j} * eigvec{j}' );
            rec_err_os = scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D;
            sq_rec_err(:, j) = sum(rec_err_os.^2, 2); 
            
        elseif (metric==3) 
            % Orthogonality
            c = C(j,:); normc = norm(c);
            for i = 1 : rows
                x = scal_X(i,:);
                normx = norm(x);
                temp = abs(x * c'); 
                dist = 1;
                if temp > eps
                    dist = 1 - temp / ( normx * normc );
                end
                distance(i,1) = mean(dist);
            end
            sq_rec_err(:,j) = distance;
            
        elseif (metric==4)
            % CVT-PCA
            sq_rec_err(:, j)=cvtpca(scal_X,eigvec{j},gamma{j});
            
        end
    end
    
    %%
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    
    % Partition the data into clusters
    [nz_X_k, nz_idx_clust, k] = partitionVQ(scal_X, idx, k);

   
    fprintf('\nClusters dimension \n');
    disp(nz_X_k);
    
    fprintf('\nThe current number of clusters is %d.\n',length(nz_X_k));
    
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
    fprintf('\nGlobal mean recontruction error at iteration n. %d equal to %d', iter, eps_rec_new);
%     fprintf('\nLocal mean recontruction error at iteration n. %d\n', iter); 
%     for j = 1 : k
%         fprintf('%d \t', eps_rec_new_clust(j));
%     end
%     fprintf('\n');
        
    %% 3) EVALUATE NEW CLUSTERS' CENTROIDS
    C_new = zeros(k, columns);        
    for j = 1 : k
        C_new(j, :) = mean(nz_X_k{j}, 1);
    end
    eps_old = eps_rec_var;
    eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
    fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);
    if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) ...
            && (n_eigs < n_eigs_max)) 
        n_eigs = n_eigs + 1;
        fprintf('\n Cluster %d dimension increased to %d \n', j,  n_eigs);
    end
    
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

    % 4) PERFORM LOCAL PCA
    if metric == 3
        nz_X_k_Cmat = cell(k, 1);
        for j = 1 : k
            thissize = size(nz_X_k{j} ,1);
            C_mat = repmat(C(j, :), thissize, 1);
            nz_X_k_Cmat{j} = nz_X_k{j} - C_mat;
        end
        [eigvec,n_eig,gamma]=performlocalpca(nz_X_k_Cmat,n_eigs);
    else
        [eigvec,n_eig,gamma]=performlocalpca(nz_X_k,n_eigs);
    end
    iter = iter + 1;
    
end

if (k == 1)
    idx = ones(rows,1);
    [nz_X_k, nz_idx_clust, k] = partitionVQ(scal_Xold, idx, k);
else
    [nz_X_k, nz_idx_clust, k] = partitionVQ(scal_Xold, idx, k);
end

if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end

% 
uncentered_nz_X_k = cell(k, 1);
for j = 1 : k
    uncentered_nz_X_k{j} = X_old(nz_idx_clust{j}, :);
end


% VARARGOUT (ADDED BY ME)
    if (nargout>2)
        varargout{1}.nz_X_k=nz_X_k;
        varargout{1}.k=k;
%         varargout{1}.eigvec=eigvec;
%         varargout{1}.
    end

close all hidden
% cd ..

end



%% Centroidal-Voronoi-Tessellations-based PCA
function delta=cvtpca(firstIn,varargin)
%This function is useful when clustering the data (before performing PCA).
%It gives a metric to evaluate the distance between a space-vector and a
%cluster's centroid.

% INPUTS CHECK:
h=whos('firstIn');
if strcmp(h.class,'pcaData') || strcmp(h.class,'podData')
    X = firstIn.Y;      
    eigvec = firstIn.modes;
    gamma = firstIn.d;
else
    X = firstIn;
    if (nargin<2)
        error('You must supply the second input: eigenvectors matrix.\n');
    end
    eigvec = varargin{1};
    if (nargin>2)
        gamma = varargin{2};
    else
        gamma = ones(size(X,1),1);
    end
end
gamma=gamma(:);

% PRE-PROCESSING:
[n, m] = size(X);   k = size(eigvec,2);
delta=zeros(n,1);

if max(size(gamma)) == 1
    eigvec_scaled = eigvec*gamma;
else
    eigvec_scaled = eigvec.*repmat(gamma,1,k);
end

% MAIN:
tol=1e-22;
for i=1:n
    x=X(i,:); 
    delta(i,1) = 1 - ( 1 / (norm(x)^2+tol) ) * ( norm(x*eigvec_scaled)^2 );
end
    
end

%% Perform Local PCA
function [eigvec,n_eig,gamma,varargout]=performlocalpca(nz_X_k,n_eigs)

    % Number of clusters
    k = max(size(nz_X_k));
    
    % Initialization of cell arrays
    eigvec = cell(k, 1);
    n_eig = cell(k, 1);
    gamma = cell(k, 1);
    pod = cell(k,1);
    
    % PCA opts
    pcaOpts.mycon=[]; pcaOpts.con=false;
    pcaOpts.cs=false;
    pcaOpts.method=1;
    pcaOpts.k=n_eigs;

    for j = 1 : k
        % Load the each cluster
        thisx=nz_X_k{j}';

        % Build the pcaData object
        pod{j}=pcaData(thisx,pcaOpts); 

        % Outputs
        eigvec{j}=pod{j}.modes;
        n_eig{j}=pod{j}.k;             
        gamma{j}=pod{j}.d;
    end
    
    % Additional Outputs
    if nargout>3
        varargout{1}=pod;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%             'Distance' ? Distance measure
%             'sqeuclidean' (default) | 'cityblock' | 'cosine' | 'correlation' | 'hamming'
%             [idx,c,sumd,sq_rec_err] = kmeans(scal_X,k,'Distance','sqeuclidean',...
%                 'EmptyAction','drop','MaxIter',1000);
    