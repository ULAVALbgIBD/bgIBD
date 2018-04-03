

function lk_result = linkage_compute(family, inner1, inner2, max_path)


% relative mapping in id, temp(:,2:4) 
% so two individuals should be changed to their local id in a family

parents = family(:, 2:4);

lk_ancestor = ancestor_matrix(parents);

lk_result = inheritance_path(inner1, inner2, parents, lk_ancestor, max_path);

end

%%



