function plot_sampled_input_space(obj, varargin)

% Check the points are there 
if isempty(obj.training_points) || isempty(obj.prediction_points)
    fprintf('No training or prediction points are available.\n');
    return;
end

% Check the dimension of the input space
dim = size(obj.training_points, 2);
if dim > 3
    fprintf('Dimension of input space is > 3.\nNo plot possible.\n');
    return;    
end

% Input space with dimension = 1
if dim == 1
    y_vec = obj.training_points * 0 + 1;
    plot(obj.training_points, y_vec, 'ko');
    hold on; grid on;
    y_vec = obj.prediction_points * 0 + 1;
    plot(obj.prediction_points, y_vec, '*');
end

% Input space with dimension = 2
if dim == 2
    plot(obj.training_points(:,1), obj.training_points(:,2), 'ko');
    hold on; grid on;
    plot(obj.prediction_points(:,1), obj.prediction_points(:,2), '*');
end

% Input space with dimension = 3
if dim == 3 
    plot3(obj.training_points(:,1), obj.training_points(:,2), obj.training_points(:,3), 'ko');
    hold on; grid on;
    plot3(obj.prediction_points(:,1), obj.prediction_points(:,2), obj.prediction_points(:,3), '*');
end

% Legend and title
legend('Train. Points', 'Pred. Points');
title('Sampled input space');

end
