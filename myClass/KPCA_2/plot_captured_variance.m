function plot_captured_variance(obj, varargin)

% Get the cumulative energy
y = [0; cumsum(obj.pca_eigenvalues)/sum(obj.pca_eigenvalues)];

% Plot
plot(0:length(y)-1, 100 * y); grid on;
hold on;

% PCA energy
obj.pca_energy = y(obj.pca_approximation_order+1);

% Mark the current captured energy
plot(obj.pca_approximation_order+1, 100 * obj.pca_energy, '*');

% Axis labels
set(get(gca, 'XLabel'), 'String', 'Number of PCs');
set(get(gca, 'YLabel'), 'String', 'Cumulative variance [%]');

% Axis limits
xlim([-1, length(y)]);
ylim([0, 105]);

% Title 
set(get(gca, 'Title'), 'String', 'Captured variance');

hold off;

end
