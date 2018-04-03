
% set partition based on neighborhood voting

function [two_cliques, new_pairs] = generate_2cliques(pairs, kinship2, posteriorIBD1, posterior, triples, clique_assignment, tolerance, range)


% preprocessing to make ibd2 consistent

ibd1 = (pairs == 1);
ibd2 = (pairs == 2);
ibd12 = (pairs >= 1);

index = 0;  %node with the most neighbors
[num, c] = size(pairs);
if( num <= 0 || num ~= c )
    disp('error in input pairs');
end

new_pairs = pairs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = sum(clique_assignment, 1);
if( ~any((temp==0)) )
    % none unassigned individual
    index = 0;
    two_cliques = [];  
    return;
end

%pick the node with maximum number of neighbors
[index, neighbor12] = biggest_degree_node(ibd12, find(temp==0), 1:num);

[index, final_p_clique, final_m_clique, final_pm_clique] = segregate_2homologous(index, pairs, triples, kinship2, posterior, posteriorIBD1, tolerance);


if( isempty(final_pm_clique) )
    disp('error in clique partition seeding');
end

% for small neighborhood partion, index might have been reassigned
two_cliques.rep_p = -index;
two_cliques.rep_m = index;

two_cliques.p_clique = final_p_clique; 
two_cliques.m_clique = final_m_clique;
two_cliques.pm_clique = final_pm_clique;


end













