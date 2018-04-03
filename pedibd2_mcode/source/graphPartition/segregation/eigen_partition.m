
function [output_p_clique output_m_clique error] = eigen_partition(p_clique, m_clique, neighbor1, ibd12, posteriorIBD12)

% will be equivalent to eigenvector, after first iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p_clique and m_clique will never shrink
% all predisposed cliques are preserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ibd2 nodes will receive equal vote from all nodes, because they are neutral
% ibd0 nodes will receive no vote, if included
% ibd1 nodes will receive vote from its own clique

% let deterministic nodes vote first

[num, c] = size(ibd12);
if( num <= 0 || c ~= num )
    error = 1;
    disp('error in neighborhood affinity');
    return;
end

error = 0;
output_p_clique = [];
output_m_clique = [];

if( isempty(neighbor1) )
    error = 1;
    disp('neighborhood empty');
    return;
end

if( ~isempty(p_clique) )
    if( ~all(ismember(p_clique, neighbor1)) )
        error = 1;
        disp('invalid kinship partition');
        return;
    end
end

if( ~isempty(m_clique) )
    if( ~all(ismember(m_clique, neighbor1)) )
        error = 1;
        disp('invalid kinship partition');
        return;
    end
end

neighbor_bit(1:num) = false;
neighbor_bit(neighbor1) = true;
p_bit(1:num) = false;
p_bit(p_clique) = true;
m_bit(1:num) = false;
m_bit(m_clique) = true;
cluster_vote(1:num) = 0;

if( any(p_bit) )
    for i = 1:num
        if( p_bit(i) == 1 )
            temp = ibd12(i,1:num);
            affinity = neighbor_bit & temp;
            repellent = neighbor_bit & (~temp);
            cluster_vote(affinity) = cluster_vote(affinity) + posteriorIBD12(i, affinity);
            cluster_vote(repellent) = cluster_vote(repellent) - 1;
        end
    end
end

% originally +/- 1, to create more discrimination
if( any(m_bit) )
    for i = 1:num
        if( m_bit(i) == 1 )
            temp = ibd12(i,1:num) == 1;
            affinity = neighbor_bit & temp;
            repellent = neighbor_bit & (~temp);
            cluster_vote(affinity) = cluster_vote(affinity) - posteriorIBD12(i, affinity);
            cluster_vote(repellent) = cluster_vote(repellent) + 1;            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bypass situation if no predisposed seeds
% ibd12 is logical matrix
nneighbors = nnz(neighbor_bit);
matrix = zeros(nneighbors, nneighbors);
matrix(ibd12(neighbor_bit,neighbor_bit)) = 1;
matrix = posteriorIBD12(neighbor_bit,neighbor_bit);
matrix(~ibd12(neighbor_bit,neighbor_bit)) = -1;
if( ~any(p_bit) && ~any(m_bit) )
    [V,D] = eig(matrix);
    [d,i] = max(diag(D));
    output_p_clique = neighbor1(V(1:nneighbors, i) >= 0);
    output_m_clique = neighbor1(V(1:nneighbors, i) < 0);
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let the remaining nodes vote
% guess their membership from the beginning

for i = 1:num
    if( neighbor_bit(i) == 1 )
        if( ~p_bit(i) && ~m_bit(i) )
            temp = ibd12(i,1:num);
            affinity = neighbor_bit & temp;
            repellent = neighbor_bit & (~temp); 
            sum1 = 0;
            sum2 = 0;
            if( any(affinity) )
                sum1 = sum(cluster_vote(affinity));
            end
            if( any(repellent) )
                sum2 = sum(cluster_vote(repellent));
            end
            if( sum1 - sum2 >= 0 )
                cluster_vote(affinity) = cluster_vote(affinity) + posteriorIBD12(i, affinity);
                cluster_vote(repellent) = cluster_vote(repellent) - 1;
            else
                cluster_vote(affinity) = cluster_vote(affinity) - posteriorIBD12(i, affinity);
                cluster_vote(repellent) = cluster_vote(repellent) + 1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_p_bit(1:num) = false;
output_p_bit(p_bit) = true;
output_m_bit(1:num) = false;
output_m_bit(m_bit) = true;

for i = 1:num
    if( neighbor_bit(i) == 1 )
        if( ~p_bit(i) && ~m_bit(i) ) 
            temp = ibd12(i,1:num);
            affinity = neighbor_bit & temp;
            repellent = neighbor_bit & (~temp);
            eigenv = 0;
            if( any(affinity) )
                eigenv = eigenv + sum(cluster_vote(affinity));
            end
            if( any(repellent) )
                eigenv = eigenv - sum(cluster_vote(repellent));
            end
            if( eigenv >= 0 ) 
                output_p_bit(i) = true;
            end
            if( eigenv < 0 )
                output_m_bit(i) = true;
            end
        end
    end
end

output_p_clique = find(output_p_bit);
output_m_clique = find(output_m_bit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end











