function getDirectKrigingErrors(obj, varargin)

obj.direct_prediction_error_observations = ...
    obj.get_errors(obj.original_data, obj.kriged_direct_data); 

obj.direct_prediction_error_variables = ...
    obj.get_errors(obj.original_data', obj.kriged_direct_data');

end

