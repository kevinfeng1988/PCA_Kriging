function getPredictionErrors(obj, varargin)

if ~isempty(obj.kriged_pca_scores)
    obj.getKpcaPredictionsErrors();
end
if ~isempty(obj.kriged_cpca_scores)
    obj.getKcpcaPredictionsErrors();
end
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.kriged_pca_scores)
    obj.getKlpcaPredictionsErrors();
end
if ~isempty(obj.local_pca) && ~isempty(obj.local_pca{1}.kriged_cpca_scores)
    obj.getKlcpcaPredictionsErrors();
end

end