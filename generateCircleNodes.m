function [interior_nodes, boundary_nodes] = generateCircleNodes(center, radius, num_interior, num_boundary)
%GENERATECIRCLENODES Generates random interior and regular boundary nodes for a circle.
%   [INTERIOR_NODES, BOUNDARY_NODES] = GENERATECIRCLENODES(CENTER, RADIUS, NUM_INTERIOR, NUM_BOUNDARY)
%   generates:
%       - NUM_INTERIOR points uniformly distributed inside a circle.
%       - NUM_BOUNDARY points evenly spaced on the circle's circumference.
%
%   Inputs:
%       center: [x0, y0], coordinates of the circle's center (default: [0, 0])
%       radius: R, radius of the circle (default: 1)
%       num_interior: N, number of interior nodes
%       num_boundary: M, number of boundary nodes
%
%   Outputs:
%       interior_nodes: Nx2 matrix of [x, y] coordinates for interior points
%       boundary_nodes: Mx2 matrix of [x, y] coordinates for boundary points

    % Set default values if not provided
    if nargin < 1
        center = [0, 0];
    end
    if nargin < 2
        radius = 1;
    end
    if nargin < 3
        num_interior = 100;
    end
    if nargin < 4
        num_boundary = 50;
    end

    x0 = center(1);
    y0 = center(2);

    %% 1. Generate Boundary Nodes (Regularly spaced on circumference)
    % Create equally spaced angles
    theta = linspace(0, 2*pi, num_boundary + 1)';
    theta(end) = []; % Remove the duplicate 2*pi point (same as 0)
    
    % Calculate boundary coordinates
    x_boundary = x0 + radius * cos(theta);
    y_boundary = y0 + radius * sin(theta);
    boundary_nodes = [x_boundary, y_boundary];

    %% 2. Generate Interior Nodes (Uniformly random inside the circle)
    % Pre-allocate array for efficiency
    interior_nodes = zeros(num_interior, 2);
    points_generated = 0;
    
    % Use rejection sampling to generate points uniformly in the circle
    while points_generated < num_interior
        % p = haltonset(1,'Skip',1e3,'Leap',1e2);
        % X0 = net(p,1);
        % Generate random points in the bounding square [-R, R] x [-R, R]
        x_candidate = x0 + radius * (2 * unifrnd(-1,1) - 1); % x in [x0-R, x0+R]
        y_candidate = y0 + radius * (2 * unifrnd(-1,1) - 1); % y in [y0-R, y0+R]
        
        % Check if the point is inside the circle
        % (x - x0)^2 + (y - y0)^2 <= R^2
        if (x_candidate - x0)^2 + (y_candidate - y0)^2 <= radius^2
            points_generated = points_generated + 1;
            interior_nodes(points_generated, :) = [x_candidate, y_candidate];
        end
    end

    %% 3. Visualize the Results (Optional)
    figure;
    hold on;
    
    % Plot interior nodes
    scatter(interior_nodes(:, 1), interior_nodes(:, 2), 40, 'b', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    % Plot boundary nodes
    scatter(boundary_nodes(:, 1), boundary_nodes(:, 2), 60, 'r', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Plot the circle for reference
    viscircles(center, radius, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 1);
    
    % Format the plot
    axis equal;
    grid on;
    box on;
    title(sprintf('Circle Nodes\nInterior: %d, Boundary: %d', num_interior, num_boundary));
    legend('Interior Nodes', 'Boundary Nodes', 'Circle Boundary', 'Location', 'best');
    xlabel('X');
    ylabel('Y');
    hold off;
end