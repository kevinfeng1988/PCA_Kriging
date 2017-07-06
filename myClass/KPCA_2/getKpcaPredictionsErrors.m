function getKpcaPredictionsErrors(obj, varargin)

% If there is no test data, return
if isempty(obj.original_data)
    return
end

if obj.is_mesh_variable
    % Get errors
    obj.kpca_prediction_error_observations = ...
        obj.get_errors(obj.original_data, obj.kpca_predictions);
    obj.kpca_prediction_error_variables = ...
        obj.get_errors(obj.original_data', obj.kpca_predictions');
    obj.kpca_prediction_error_variables = ...
        aveVar(obj.kpca_prediction_error_variables, obj.variable_names, obj.mesh);
else
    % Rebuild the Y matrices
    n_samples = size(obj.kpca_predictions,2) / size(obj.mesh,1);
    Y_original = leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, n_samples);
    Y_predicted = leaveStateVars(obj.kpca_predictions, obj.mesh, obj.variable_names, n_samples);
    % Get errors
    obj.kpca_prediction_error_observations = obj.get_errors(Y_original, Y_predicted);
    obj.kpca_prediction_error_variables = obj.get_errors(Y_original', Y_predicted');
    obj.kpca_prediction_error_variables = ...
        aveVar(obj.kpca_prediction_error_variables, obj.variable_names, obj.mesh);
end

end




