function [varargout] = get_variable_errors(obj, logical, var, model, varargin)

% Provide this info
if ~isempty(varargin)
    is_mesh_var = varargin{1};
else
    is_mesh_var = true;
end

% Get predicted/reconstructed variable
if logical
    Y = obj.get_variable(obj.prediction_points, var, model, is_mesh_var);
    Y_original = obj.get_variable(obj.prediction_points, var, 'data', is_mesh_var);
else
    Y = obj.get_variable(obj.training_points, var, model, is_mesh_var);
    Y_original = obj.get_variable(obj.training_points, var, 'data', is_mesh_var);
end

% Get error
varargout{1} = obj.get_errors(Y_original, Y);

end


