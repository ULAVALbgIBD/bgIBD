function [pm_clique, clique_assignment, rep_list, pairs] = assign_one_component(pairs, kinship2, posteriorIBD1, posterior, triples, clique_assignment, rep_list, verify, debug, range)

% do not change tolerance for dense markers
tolerance = 0.5;
[two_cliques, pairs] = generate_2cliques(pairs, kinship2, posteriorIBD1, posterior, triples, clique_assignment, tolerance, range);

if( isempty(two_cliques) )
    % all individuals are fully assigned
    pm_clique = [];
    return;
else
    pm_clique = two_cliques.pm_clique;
end

[clique_assignment, rows] = add_clique(clique_assignment, two_cliques.pm_clique, two_cliques.p_clique);
rep_list(rows) = two_cliques.rep_p;
[clique_assignment, rows] = add_clique(clique_assignment, two_cliques.pm_clique, two_cliques.m_clique);
rep_list(rows) = two_cliques.rep_m;

if( debug )
    p_clique = two_cliques.p_clique;
    m_clique = two_cliques.m_clique;
    pm_clique = two_cliques.pm_clique;
    sort([p_clique, pm_clique])
    range(sort([p_clique, pm_clique]))
    pairs(sort([p_clique, pm_clique]),sort([p_clique, pm_clique]))
    sort([m_clique, pm_clique])
    range(sort([m_clique, pm_clique]))
    pairs(sort([m_clique, pm_clique]),sort([m_clique, pm_clique]))
end


i = 0;
while(1)
    if( debug )
        i = i + 1
    end
    [one_clique, pairs] = generate_1clique(pairs, kinship2, posteriorIBD1, posterior, triples, clique_assignment, tolerance, verify, debug, range);
    if( isempty(one_clique) )
        break;
    end
    [clique_assignment, rows] = add_clique(clique_assignment, one_clique.pm_clique, one_clique.new_clique);
    rep_list(rows) = one_clique.rep;
    if( debug )
        display_cliques(clique_assignment, pairs, pairs, range);
    end
end

if( debug )
    clique_assignment
end
