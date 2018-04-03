

function [one_clique, new_pairs] = generate_1clique(pairs, kinship2, posteriorIBD1, posterior, triples, clique_assignment, tolerance, verify, debug, range)

% preprocessing to make ibd2 consistent

ibd1 = (pairs == 1);
ibd2 = (pairs == 2);
ibd12 = (pairs >= 1);

index = 0;  %node with the most neighbors
num = length(pairs(1,:));

new_pairs = pairs;
rep = index;

temp = sum(clique_assignment, 1);
if( isempty(find(temp==1)) )
    one_clique = [];
    index = 0;
    pm_clique = [];
    clique = [];
    rep = [];
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pick the node with maximum number of neighbors
% pick large cliques first with good evidence
[index, neighbor12] = biggest_degree_node(ibd12, find(temp==1), find(temp>0));
neighbor1 = find(ibd1(index,:) == 1);
neighbor2 = find(ibd2(index,:) == 1);
neighbor12 = find(ibd12(index,:) == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% index is the biggest node to be expanded
% will remove more ambiguity in this step
% consider ibd1 relationship to index only
% ibd2 relationship cannot be biased towards each clique 
[p_clique_forced, m_clique_forced] = initial_partition(index, neighbor1, kinship2, posteriorIBD1);


pre_clique = find(clique_assignment(find(clique_assignment(:,index) == 1), :) == 1);

pm_clique = neighbor2;
% pm_clique must be included in the pre_clique, otherwise errors



p_len = length(p_clique_forced);
m_len = length(m_clique_forced);
p_shared = length(intersect(p_clique_forced, pre_clique));
m_shared = length(intersect(m_clique_forced, pre_clique));
pm_len = length(pm_clique);
pm_shared = length(intersect(pm_clique, pre_clique));



if( debug )
    if( p_shared ~= 0 && p_shared ~= p_len )
        disp('possible inconsistency in global IBD');
    end
    if( m_shared ~= 0 && m_shared ~= m_len )
        disp('possible inconsistency in global IBD');
    end
    if( pm_shared ~= pm_len )
        disp('possible inconsistency in global IBD');
    end

    if( p_len ~= 0 && m_len ~= 0 )
        if( p_shared == 0 && m_shared == 0 )
            disp('possible inconsistency in global IBD');
        end
        if( p_shared ~= 0 && m_shared ~= 0 )
            disp('possible inconsistency in global IBD');
        end
    end
end


if( p_len ~= 0 && m_len ~= 0 )
% if not fully overlap with one, there is also inconsistency

    if( p_shared >= m_shared )
        m_clique = setdiff(neighbor12, [pre_clique, p_clique_forced, pm_clique]);
        p_clique = setdiff(pre_clique, pm_clique);
        clique = m_clique;
        rep = index;
    end

    if( m_shared > p_shared )
        p_clique = setdiff(neighbor12, [pre_clique, m_clique_forced, pm_clique]);
        m_clique = setdiff(pre_clique, pm_clique);
        clique = p_clique;
        rep = -index;
    end
    
end


if( p_len == 0 && m_len == 0 )
    % arbitrarily assign
    m_clique = setdiff(neighbor12, [pre_clique, p_clique_forced, pm_clique]);
    p_clique = setdiff(pre_clique, pm_clique);
    clique = m_clique;   
    rep = index;
end

if( p_len > 0 && m_len == 0 )
    if( p_shared == 0 )
        p_clique = setdiff(neighbor12, [pre_clique, m_clique_forced, pm_clique]);
        m_clique = setdiff(pre_clique, pm_clique);
        clique = p_clique;
        rep = -index;
    else
        m_clique = setdiff(neighbor12, [pre_clique, p_clique_forced, pm_clique]);
        p_clique = setdiff(pre_clique, pm_clique);
        clique = m_clique;
        rep = index;
    end
end

if( p_len == 0 && m_len > 0 )
    if( m_shared == 0 )
        m_clique = setdiff(neighbor12, [pre_clique, p_clique_forced, pm_clique]);
        p_clique = setdiff(pre_clique, pm_clique);
        clique = m_clique;  
        rep = index;
    else
        p_clique = setdiff(neighbor12, [pre_clique, m_clique_forced, pm_clique]);
        m_clique = setdiff(pre_clique, pm_clique);
        clique = p_clique;
        rep = -index;
    end
end


clique = prune_1clique(pm_clique, clique, ibd12, posterior, tolerance);

one_clique.pm_clique = pm_clique;
one_clique.new_clique = clique;
one_clique.rep = rep;

if( debug )

    sort([p_clique, pm_clique])
    range(sort([p_clique, pm_clique]))
    pairs(sort([p_clique, pm_clique]),sort([p_clique, pm_clique]))
    sort([m_clique, pm_clique])
    range(sort([m_clique, pm_clique]))
    pairs(sort([m_clique, pm_clique]),sort([m_clique, pm_clique]))

end




end













