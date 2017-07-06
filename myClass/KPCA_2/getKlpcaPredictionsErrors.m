function getKlpcaPredictionsErrors(obj, varargin)

% If there is no test data, return
if isempty(obj.original_data)
    return
end

if obj.is_mesh_variable
    % Get errors
    obj.klpca_prediction_error_observations = ...
        obj.get_errors(obj.original_data, obj.klpca_predictions);
    obj.klpca_prediction_error_variables = ...
        obj.get_errors(obj.original_data', obj.klpca_predictions');
    obj.klpca_prediction_error_variables = ...
        aveVar(obj.klpca_prediction_error_variables, obj.variable_names, obj.mesh);
else
    % Rebuild the Y matrices
    n_samples = size(obj.klpca_predictions,2) / size(obj.mesh,1);
    Y_original = leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, n_samples);
    Y_predicted = leaveStateVars(obj.klpca_predictions, obj.mesh, obj.variable_names, n_samples);
    % Get errors
    obj.klpca_prediction_error_observations = obj.get_errors(Y_original, Y_predicted);
    obj.klpca_prediction_error_variables = obj.get_errors(Y_original', Y_predicted');
    obj.klpca_prediction_error_variables = ...
        aveVar(obj.klpca_prediction_error_variables, obj.variable_names, obj.mesh);
end

end