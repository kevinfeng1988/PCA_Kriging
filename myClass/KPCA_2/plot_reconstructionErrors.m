function [varargout] = plot_reconstructionErrors(obj, varargin)

% Check the points are there 
if isempty(obj.training_points) || isempty(obj.prediction_points)
    fprintf('No training or prediction points are available.\n');
    return;
end

% Check the dimension of the input space
if obj.is_mesh_variable
    x_train = obj.training_points;
    x_test = obj.prediction_points;
    dim = size(x_train, 2);
else
    x_train = unique( obj.training_points(:,size(obj.mesh,2)+1:end) );
    x_test = unique( obj.prediction_points(:,size(obj.mesh,2)+1:end) );
    dim = size(obj.training_points, 2) - size(obj.mesh,2);
end
if dim > 3
    fprintf('Dimension of input space is > 3.\nNo plot possible.\n');
    return;    
end

% Input
model = 'pca';
if ~isempty(varargin)
    model = varargin{1};
    I = strcmp({'pca', 'lpca', 'cpca', 'lcpca', 'pca_var', 'lpca_var', 'cpca_var', 'lcpca_var'}, model);
    if sum(I) < 1
        error('Model not valid.')
    end
end

% Get correct errors
if strcmp(model, 'pca')
    y_err = obj.pca_reconstruction_error_observations;
elseif strcmp(model, 'lpca')
    y_err = obj.local_pca_reconstruction_error_observations;
elseif strcmp(model, 'cpca')
    y_err = obj.cpca_reconstruction_error_observations;
elseif strcmp(model, 'lcpca')
    y_err = obj.local_cpca_reconstruction_error_variables;
elseif strcmp(model, 'pca_var')
    y_err = obj.pca_reconstruction_error_variables;
elseif strcmp(model, 'lpca_var')
    y_err = obj.local_pca_reconstruction_error_variables;
elseif strcmp(model, 'cpca_var')
    y_err = obj.cpca_reconstruction_error_variables;
elseif strcmp(model, 'lcpca_var')
    y_err = obj.local_cpca_reconstruction_error_variables;
end

% If present, run only for one specific variable
var = 'all';
if length(varargin) > 1
    var = varargin{2};
    y_err = obj.get_variable_errors(false, var, model, obj.is_mesh_variable);
end

% Number of points 
m = size(x_train, 1);

% Symbols
s1 = 'ko';
s2 = 'bx';

% Plot variables?
is_var = strcmp(model(end-2:end), 'var');

% Start new figure
figure();
if ~is_var
    % Input space with dimension = 1
    if dim == 1
        [x, I] = sortrows([x_train; x_test]);
        y = [y_err; 0*x_test]; y = y(I);
        bar(x, 100 * y); 
        grid on; hold on;
        ylabel('Error [%]');
    end

    % Input space with dimension = 2
    if dim == 2
        plot(x_train(:,1), x_train(:,2), 'ko');
        hold on; grid on;
        plot(x_test(:,1), x_test(:,2), s2);
        hold on; grid on;
        % Plot errors
        for ii = 1 : m
            x = x_train(ii,1);
            y = x_train(ii,2);
            txt = [' ', num2str(100 * y_err(ii),3), '%'];
            text(x,y,txt);
        end
    end

    % Input space with dimension = 3
    if dim == 3 
        plot3(x_train(:,1), x_train(:,2), x_train(:,3), 'ko');
        hold on; grid on;
        plot3(x_test(:,1), x_test(:,2), x_test(:,3), s2);
        hold on; grid on;
        % Plot errors
        for ii = 1 : m
            x = x_train(ii,1);
            y = x_train(ii,2);
            z = x_train(ii,3);
            txt = [' ', num2str(100 * y_err(ii),3), '%'];
            text(x,y,z,txt);
        end
    end
else
    bar(1:length(y_err), 100 * y_err); 
    grid on; hold on;
    ylabel('Error [%]');
end

% Legend and title
if strcmp(var, 'all')
    tit = ['Reconstruction errors: ', model];
else
    tit = ['Reconstruction errors (', var, '): ', model];
end
if strcmp(model(1), 'l')
    tit = [tit, ' ', num2str(obj.number_of_clusters)];
end
title(tit);
hold off;

end

