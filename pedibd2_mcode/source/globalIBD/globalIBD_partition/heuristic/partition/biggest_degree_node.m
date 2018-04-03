
function [index, max_neighbor] = biggest_degree_node(relationship, domain, range)

max_neighbor = 0;
neighbor_count = 0;
num = length(relationship);
num_neighbor(1:num) = 0;
for i = domain
    neighbor = (find(relationship(i,:)==1));
    num_neighbor(i) = length(intersect(neighbor, range));
    if( num_neighbor(i) > neighbor_count )
        neighbor_count = num_neighbor(i);
        max_neighbor = neighbor;
        index = i;
    end
end

end