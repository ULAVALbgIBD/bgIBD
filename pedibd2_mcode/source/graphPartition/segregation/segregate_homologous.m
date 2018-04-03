function [final_p_clique, final_m_clique, final_pm_clique, error] ...
    =  segregate_homologous ...
    (index, ...
    pairs, ...
    triples, ...
    kinship2, singleLINEAL, ...
    posterior, posteriorIBD1, tolerance)

final_p_clique = [];
final_m_clique = [];
final_pm_clique = [];
error = 0;
[nGENO, c] = size(pairs);
if( nGENO <= 0 || nGENO ~= c )
    error = 1;
    disp('error in input pairs');
    return;
end


% logical matrix
ibd1 = (pairs == 1);
ibd2 = (pairs == 2);
ibd12 = (pairs >= 1);

neighbor1 = find(ibd1(index,:));
neighbor2 = find(ibd2(index,:));
neighbor12 = find(ibd12(index,:));

if( ~isempty(neighbor1) )
    [p_clique_forced, m_clique_forced error] ...
        = kinship_partition(index, neighbor1, kinship2, posteriorIBD1);
    if( error ~= 0 )
        disp('impossible chromosome sharing given this family structure');
        return;
    end
else
    p_clique_forced = [];
    m_clique_forced = [];
end

% modify ibd, given reference point at individual index
triple_viewpair = triples.triple_viewpair;
for i = 1:nGENO
    for j = 1:nGENO
        if( ibd1(index,i) && ibd1(index,j) && ibd1(i,j) )
            if( triple_viewpair(index,i,j) > 0 )
                ibd12(i,j) = true;
            end
            if( triple_viewpair(index,i,j) < 0 )
                ibd12(i,j) = false;
            end
        end
    end
end

% modify ibd, given forced ancestral partition
ibd12(p_clique_forced,p_clique_forced) = true;
ibd12(m_clique_forced,m_clique_forced) = true;
ibd12(p_clique_forced,m_clique_forced) = false;
ibd12(m_clique_forced,p_clique_forced) = false;
% does not exhaust all pedigree information

temp = reshape(singleLINEAL(index,:,1:2), [nGENO,2]);
transmission = unique(temp, 'rows');
for i = 1:size(transmission,1)
    if( all(transmission(i,:) > 0) )
        % self coding, ambiguous
        group = ismember(temp, transmission(i,:), 'rows') & ibd12(index, :)';
        ibd12(group, group) = true;
    end
end

if( ~isempty(neighbor1) )
    % if no prio partition, eigen_partition has one degree of freedom
    % set by random
    posteriorIBD12 = posterior(1:nGENO,1:nGENO,2) + posterior(1:nGENO,1:nGENO,3);
    [p_clique_expanded, m_clique_expanded error] = eigen_partition(p_clique_forced, m_clique_forced, neighbor1, ibd12, posteriorIBD12);
    if( error ~= 0 )
        disp('error in DNAphoresis');
        return;
    end
    if( isempty(m_clique_forced) && isempty(m_clique_forced) && isempty(p_clique_expanded) )
        % ambiguous assigning, only one side is not empty
        % let the paternal side be non-empty
        % pre-order, supplement reassign, founder allele
        % reassign() skips very short recombination block
        p_clique_expanded = m_clique_expanded;
        m_clique_expanded = [];
    end
else
    p_clique_expanded = p_clique_forced;
    m_clique_expanded = m_clique_forced;
end


final_p_clique = p_clique_expanded;
final_m_clique = m_clique_expanded;
final_pm_clique = neighbor2;

if( isempty(final_pm_clique) )
    error = 1;
    disp('error in clique partition seeding');
    return;
end

end






