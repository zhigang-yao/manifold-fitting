% 1. Randomly generate 10000 points on Möbius strip
num_points = 10000;
R = 1; % Radius of the strip
w = 0.2; % Width of the strip

theta = 2*pi*rand(num_points, 1);
phi = -w + 2*w*rand(num_points, 1);

x = (R + phi.*cos(theta/2)).*cos(theta);
y = (R + phi.*cos(theta/2)).*sin(theta);
z = phi.*sin(theta/2);

points = [x y z];

% Plot the points
figure;
scatter3(points(:,1), points(:,2), points(:,3), 'filled');

% 2. Construct knn adjacency graph
knn_graph = knnsearch(points, points, 'K', 11, 'Distance', 'euclidean');
knn_graph = knn_graph(:, 2:end); % Exclude self

% Get first point
start_point = points(1, :);

% Initialize distances to infinity
distances = inf(num_points, 1);

% Distance from start point to itself is 0
distances(1) = 0;

% Set of unvisited points
unvisited = 1:num_points;

% While there are still unvisited points
while ~isempty(unvisited)
    % Find closest point with minimum distance
    [~, idx] = min(distances(unvisited));
    closest = unvisited(idx);

    % Mark closest point as visited
    unvisited(idx) = [];

    % Update distances for unvisited neighbors
    for i = knn_graph(closest, :)
        if ismember(i, unvisited)
            new_dist = distances(closest) + norm(points(closest, :) - points(i, :));
            if new_dist < distances(i)
                distances(i) = new_dist;
            end
        end
    end
end

% Plot the points with color indicating distance
figure;
scatter3(points(:,1), points(:,2), points(:,3), [], distances, 'filled');
colorbar;
