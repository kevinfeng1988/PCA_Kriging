function getLcpcaErrors(obj, varargin)

obj.local_cpca_reconstruction_error_observations = ...
    obj.get_errors(obj.training_data, obj.local_recovered_data_cpca); 

obj.local_cpca_reconstruction_error_variables = ...
    obj.get_errors(obj.training_data', obj.local_recovered_data_cpca');

end