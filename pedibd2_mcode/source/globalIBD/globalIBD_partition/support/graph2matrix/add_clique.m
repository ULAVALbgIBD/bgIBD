
function [clique_assignment rows] = add_clique(clique_assignment, clique_seed, clique_expand)


        
    clique = [clique_seed, clique_expand];
    temp = sum(clique_assignment(:, clique), 1);
    full_assigned = clique(temp == 2);
    if( isempty(full_assigned) )
        [clique_assignment, rows] = assign_clique(clique_assignment, clique);
    else
        % if part of clique fully assigned, assign the remaining
        [clique_assignment, rows] = assign_clique(clique_assignment, setdiff(clique, full_assigned));
    end
    

    
end

