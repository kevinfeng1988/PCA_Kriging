function plot_eigenspectrum(obj, varargin)

% Get the cumulative energy
y = obj.pca_eigenvalues;

% Plot
semilogy(y); grid on;

% Axis labels
xlabel('Number of eigenvalues');
ylabel('Value [-]');

% Title 
title('Eigenvalue spectrum');

end