function getKlcpcaPredictionsErrors(obj, varargin)

if obj.is_mesh_variable
    % Get errors
    obj.klcpca_prediction_error_observations = ...
        obj.get_errors(obj.original_data, obj.klcpca_predictions);
    obj.klcpca_prediction_error_variables = ...
        obj.get_errors(obj.original_data', obj.klcpca_predictions');
    obj.klcpca_prediction_error_variables = ...
        aveVar(obj.klcpca_prediction_error_variables, obj.variable_names, obj.mesh);
else
    % Rebuild the Y matrices
    n_samples = size(obj.klpca_predictions,2) / size(obj.mesh,1);
    Y_original = leaveStateVars(obj.original_data, obj.mesh, obj.variable_names, n_samples);
    Y_predicted = leaveStateVars(obj.klcpca_predictions, obj.mesh, obj.variable_names, n_samples);
    % Get errors
    obj.klcpca_prediction_error_observations = obj.get_errors(Y_original, Y_predicted);
    obj.klcpca_prediction_error_variables = obj.get_errors(Y_original', Y_predicted');
end

end