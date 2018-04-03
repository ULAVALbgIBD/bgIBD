function [final_p_clique, final_m_clique, final_pm_clique] = prune_2cliques_complex(vote_both, kinship2, posteriorIBD1, posterior, pairs)


final_p_clique = [];
final_m_clique = [];
final_pm_clique = [];
len = length(vote_both);

if( len == 1 )
    final_pm_clique = vote_both;
    return;
else
    final_pm_clique = vote_both(1);
    %aribitrarily assign first
end

p_clique_forced = [];
m_clique_forced = [];
pm_clique_forced = [];


for i = 1:len
    index = vote_both(i);
    [p_clique_forced, m_clique_forced] = initial_partition(index, vote_both, kinship2, posteriorIBD1);
    if( ~isempty(p_clique_forced) || ~isempty(m_clique_forced) )
        pm_clique_forced = index;
        final_pm_clique = [];
        break;
    end
end

for j = 1:length(p_clique_forced)
    p_node = p_clique_forced(j);
    p_expand = vote_both(find(pairs(p_node,vote_both)==2));
    final_p_clique = [final_p_clique, p_expand];
end
final_p_clique = unique(final_p_clique);
% remove duplicates

for j = 1:length(m_clique_forced)
    m_node = m_clique_forced(j);
    m_expand = vote_both(find(pairs(m_node,vote_both)==2));
    final_m_clique = [final_m_clique, m_expand];
end
final_m_clique = unique(final_m_clique);

for j = 1:length(pm_clique_forced)
    pm_node = pm_clique_forced(j);
    pm_expand = vote_both(find(pairs(pm_node,vote_both)==2));
    final_pm_clique = [final_pm_clique, pm_expand];
end
final_pm_clique = unique(final_pm_clique);

remainder = setdiff(vote_both, [final_p_clique, final_m_clique, final_pm_clique]);
if( ~isempty(remainder) )
    if( isempty(final_p_clique) )
        final_p_clique = remainder;
    else
        if( isempty(final_m_clique) )
            final_m_clique = remainder;
        else
            disp('small clique inconsistency');
        end
    end
end

end
























