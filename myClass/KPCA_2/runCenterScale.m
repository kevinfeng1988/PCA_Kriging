function runCenterScale(obj, varargin)

% Main
if ~obj.center_scale
    obj.centered_scaled_data = obj.training_data;
    obj.mean_column = 0; 
    obj.scaling_factors = ones(size(obj.training_data,1), 1);
else
    [obj.centered_scaled_data, obj.mean_column, obj.scaling_factors] = ...
        zscore(obj.training_data, 0, 2); 
end

end