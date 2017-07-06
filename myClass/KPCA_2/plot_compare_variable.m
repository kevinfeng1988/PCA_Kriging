function plot_compare_variable(obj, var, point)
%% Input

% Input: var
if ~isa(var,'char')
    var = obj.variable_names{var};
end


%% Main

% Initialize legend and title
legend_ = {'Data'};
title_ = var;

% Get vectors to plot
y_data = obj.get_variable(point, var, 'data', obj.is_mesh_variable);
y_pca = obj.get_variable(point, var, 'pca', obj.is_mesh_variable);
y_cpca = obj.get_variable(point, var, 'cpca', obj.is_mesh_variable);
y_lpca = obj.get_variable(point, var, 'lpca', obj.is_mesh_variable);
y_lcpca = obj.get_variable(point, var, 'lcpca', obj.is_mesh_variable);

% Get x-axis data
if size(obj.mesh, 2) == 1
    x = obj.mesh;
else
    x = 1:length(y_data); 
end

% Plot
plot(x, y_data, '-'); hold on;
if ~isempty(y_pca)
    plot(x, y_pca, '--');
    legend_{end+1} = 'PCA';
end
if ~isempty(y_cpca)
    plot(x, y_cpca, '--');
    legend_{end+1} = 'CPCA'; % Append to the legend
end
if ~isempty(y_lpca)
    plot(x, y_lpca, '--');
    legend_{end+1} = 'LPCA'; % Append to the legend
end
if ~isempty(y_lcpca)
    plot(x, y_lcpca, '--');
    legend_{end+1} = 'LCPCA'; % Append to the legend
end

% Plot settings
grid on;
legend(legend_);
title(title_);


end



