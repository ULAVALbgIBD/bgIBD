function [index final_p_clique, final_m_clique, final_pm_clique] =  segregate_2homologous(index, pairs, triples, kinship2, posterior, posteriorIBD1, tolerance)

[num, c] = size(pairs);
if( num <= 0 || num ~= c )
    disp('error in input pairs');
end

ibd1 = (pairs == 1);
ibd2 = (pairs == 2);
ibd12 = (pairs >= 1);

neighbor1 = find(ibd1(index,:) == 1);
neighbor2 = find(ibd2(index,:) == 1);
neighbor12 = find(ibd12(index,:) == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p_clique_forced, m_clique_forced] = initial_partition(index, neighbor1, kinship2, posteriorIBD1);

valid_triples = triples.valid_triples;
triple_viewpair = triples.triple_viewpair;
for i = 1:num
    for j = 1:num
        triple_index = valid_triples(index,i,j);
        if( ibd1(index,i) && ibd1(index,j) && ibd1(i,j) )
            if( triple_index > 0 )
                if( triple_viewpair(index,i,j) > 0 )
                    ibd12(i,j) = 1;
                end
                if( triple_viewpair(index,i,j) < 0 )
                    ibd12(i,j) = 0;
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p_clique_expanded, m_clique_expanded, vote_both, vote_neither] = vote_partition(p_clique_forced, m_clique_forced, neighbor12, ibd12);


if( isempty(p_clique_expanded) && isempty(m_clique_expanded) )
    [final_p_clique, final_m_clique, final_pm_clique] = prune_2cliques_complex(vote_both, kinship2, posteriorIBD1, posterior, pairs);
else
    [final_p_clique, final_m_clique, final_pm_clique] = prune_2cliques(index, neighbor12, ibd12, p_clique_forced, p_clique_expanded, m_clique_forced, m_clique_expanded, vote_both, posterior, tolerance);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( isempty(final_pm_clique) )
    disp('error in clique partition seeding');
end

% for small neighborhood partion, index might have been reassigned

end