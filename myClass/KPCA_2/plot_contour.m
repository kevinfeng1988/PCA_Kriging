function [varargout] = plot_contour(obj, variable, point, varargin)
% DESCRIPTION: Works only for Zhy Li's case.


if size(obj.mesh,2) ~= 2
    error('It looks like there is no 2D mesh.');
end

if size(point,1) > 1
    error('Too many points provided.');
end

if ~isempty(varargin)
    model = varargin{1};
else
    model = 'data';
end

% Get vectors to plot
if ischar(variable)
    y_data = obj.get_variable(point, variable, model);
else
    error(' ');
end

method = 'natural';
cont = false;

% Get unique coordinates
x = unique(obj.mesh(:,1));
y = unique(obj.mesh(:,2));

% Make the mesh more dense 
deltax = abs(x(2) - x(1));
deltay = abs(y(2) - y(1));
delta = min([deltax, deltay]) / 3.0;
x = min((obj.mesh(:,1))) : delta : max((obj.mesh(:,1)));
y = min((obj.mesh(:,2))) : delta : max((obj.mesh(:,2)));

% Plot original data
[xx, yy] = meshgrid(x, y);
F = scatteredInterpolant(obj.mesh(:,1), obj.mesh(:,2), y_data, method);
zz = F(xx, yy);
figure();
if cont
    y_min = min(y_data); y_max = max(y_data); % Range
    v = y_min + (0 : .05 : 1) * (y_max - y_min); % Color values
    [~, h, ~] = contourf(xx, yy, zz, v); 
    set(h(:), 'LineStyle', 'none');
else
    imshow(zz, [], 'XData', x, 'YData', y, 'Colormap', jet);
    axis xy; % Flip the figure
    axis on; 
end
xlabel('x [m]');
ylabel('y [m]');
title([variable, ' - ', model]);

if nargout > 0
    varargout{1} = zz;
end


end



% Local functions
function idx = get_var_dom(A, n, var_i)

idx = A(1+n*(var_i-1):n*var_i, :);

end






