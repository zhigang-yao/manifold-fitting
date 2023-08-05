k = 20; idx_point = 1;

NumSample = 10000;
samples = randn(3, NumSample);
samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(3, NumSample);

sigma = 0.06;
points = manfit_ours(samples', sigma, samples');

G = kgraph(points,k);
distances = Geodesic_Distance(points,G, idx_point);

% Plot the points with color indicating distance
figure;
scatter3(points(:,1), points(:,2), points(:,3), [], distances, 'filled');
colorbar;