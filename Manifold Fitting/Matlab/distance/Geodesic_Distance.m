function distances = Geodesic_Distance(points,G, idx_point)


distances = inf(size(points,1), 1);

distances(idx_point) = 0;

unvisited = 1:size(points,1);

while ~isempty(unvisited)
    [~, idx] = min(distances(unvisited));
    closest = unvisited(idx);

    unvisited(idx) = [];

    for i = G(closest, :)
        if ismember(i, unvisited)
            new_dist = distances(closest) + norm(points(closest, :) - points(i, :));
            if new_dist < distances(i)
                distances(i) = new_dist;
            end
        end
    end
end