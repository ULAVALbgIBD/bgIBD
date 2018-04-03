function [clique] = generate_clique(pairs, priority)

% preprocessing to make ibd2 consistent

ibd1 = (pairs == 1);
ibd2 = (pairs == 2);
max_neighbor = [];
neighbor_count = 0;
index = 0;

for i = 1:length(ibd1(1,:))
    neighbor = (find(ibd1(i,:)==1));
    if( length(neighbor) > neighbor_count )
        neighbor_count = length(neighbor);
        max_neighbor = neighbor;
        index = i;
    end
end



temp = priority(index,:);
[temp1, map] = sort(temp, 'descend');
mask(1:length(map)) = 0;
mask(max_neighbor) = 1;
count = 0;
for i = 1:length(map)
    if( mask(map(i)) == 1 )
        count = count + 1;
        max_neighbor(count) = map(i);
    end
end



clique = index;

for i = 1:length(max_neighbor)
    full = 1;
    for j = 1:length(clique)
        if( ibd1(clique(j),max_neighbor(i)) ~= 1 )
            full = 0;
        end
    end
    if( full == 1 )
        clique = union(clique,max_neighbor(i));
    end
end

clique2 = setdiff(find(ibd1(index,:)==1), clique);

end
