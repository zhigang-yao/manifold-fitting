function G = kgraph(points,k)
% 2. Construct knn adjacency graph
G = knnsearch(points, points, 'K', k, 'Distance', 'euclidean');
G = G(:, 2:end); % Exclude self