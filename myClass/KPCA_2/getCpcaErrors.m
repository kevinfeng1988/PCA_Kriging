function getCpcaErrors(obj, varargin)

obj.cpca_reconstruction_error_observations = ...
    obj.get_errors(obj.training_data, obj.cpca_recovered_data); 

obj.cpca_reconstruction_error_variables = ...
    obj.get_errors(obj.training_data, obj.cpca_recovered_data);

end