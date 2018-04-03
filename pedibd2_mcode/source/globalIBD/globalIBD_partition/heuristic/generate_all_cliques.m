function [clique_config] = generate_all_cliques(pairs, kinship2, posteriorIBD1, posterior, triples, verify, debug, range)



new_pairs = reconcile_pairs(pairs, posterior, range);


clique_assignment(1:length(pairs)) = 0;
rep_list(1) = 0;

while(1)
    [pm_clique, clique_assignment, rep_list, new_pairs] = assign_one_component(new_pairs, kinship2, posteriorIBD1, posterior, triples, clique_assignment, rep_list, verify, debug, range);
    if( isempty(pm_clique) )
        break;
    end
end



temp = sum( clique_assignment, 1 );
len = length(find(temp ~= 2));


if( len > 0 )
    disp(['inconsistent cliques: ', num2str(len)]);
end

if( debug )
    if( isempty(verify) )
        display_cliques(clique_assignment, pairs, range);
    else
        subplot(1,2,1);        
        title('graph to viterbi');
        display_cliques(clique_assignment, pairs, range);
        subplot(1,2,2);
        title('graph to likelihood');
        display_cliques(clique_assignment, allele2pair(verify), range);
    end
end


clique_config.rep_list = rep_list;
clique_config.assignment = clique_assignment;
if( rep_list(1) ~= 0 )
    clique_config.num_cliques = length(rep_list);
else
    clique_config.num_cliques = 0;
end
if( length(clique_assignment(:,1)) ~= clique_config.num_cliques + 1 )
    disp('error in clique assignment');
end

% rep_list may reflect which parental allele of index is assigned
% for individuals receiving this assignment, parental origin is not
% available
% refine this to faciliate imprinting analysis

end











