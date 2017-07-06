function getLpcaErrors(obj, varargin)

if obj.is_mesh_variable
    % Get errors
    obj.local_pca_reconstruction_error_observations = ...
        obj.get_errors(obj.training_data, obj.local_recovered_data_pca); 
    obj.local_pca_reconstruction_error_variables = ...
        obj.get_errors(obj.training_data', obj.local_recovered_data_pca');
    obj.local_pca_reconstruction_error_variables = ...
        aveVar(obj.local_pca_reconstruction_error_variables, obj.variable_names, obj.mesh);
else
    % Rebuild the Y matrices
    n_samples = size(obj.training_data, 2) / size(obj.mesh,1);
    Y_original = leaveStateVars(obj.training_data, obj.mesh, obj.variable_names, n_samples);
    Y_predicted = leaveStateVars(obj.local_recovered_data_pca, obj.mesh, obj.variable_names, n_samples);
    % Get errors
    obj.local_pca_reconstruction_error_variables = obj.get_errors(Y_original', Y_predicted');
    obj.local_pca_reconstruction_error_observations = obj.get_errors(Y_original, Y_predicted);
    obj.local_pca_reconstruction_error_variables = ...
        aveVar(obj.local_pca_reconstruction_error_variables, obj.variable_names, obj.mesh);
end

end