function [varargout] = plot_contour(mesh, y_data, varargin)
% DESCRIPTION: Works only for Zhy Li's case.


cont = false;

% Get unique coordinates
x = unique(mesh(:,1));
y = unique(mesh(:,2));

% Make the mesh more dense 
deltax = abs(x(2) - x(1));
deltay = abs(y(2) - y(1));
delta = min([deltax, deltay]) / 3.0;
x = mesh(1,1) : delta : mesh(end,1);
y = mesh(1,2) : delta : mesh(end,2);

% Plot original data
[xx, yy] = meshgrid(x, y);
F = scatteredInterpolant(mesh(:,1), mesh(:,2), y_data);
zz = F(xx, yy);
figure();
if cont
    contourf(xx, yy, zz);
else
    image_plot = mat2gray(zz); 
%     range_in = [min(min(image_plot)), max(max(image_plot))];
%     range_out = [min(y_data), max(y_data)];
%     image_plot = imadjust(image_plot, range_in, range_out);
    imshow(image_plot, 'XData', x, 'YData', y, 'Colormap', jet); axis xy;
    xticks = x(1) : x(2)-x(1) : x(end);
    xticklabels = x(1) : x(2)-x(1) : x(end);
    set(gca, 'xtick', xticks);
    set(gca, 'xticklabel', xticklabels);
end

xlabel('x [m]');
ylabel('y [m]');


end