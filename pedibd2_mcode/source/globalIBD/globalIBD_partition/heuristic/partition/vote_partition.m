
function [output_p_clique output_m_clique vote_both vote_neither] = vote_partition(p_clique, m_clique, neighbor12, ibd12)

% will be equivalent to eigenvalue method, given first iteration


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p_clique and m_clique will never shrink
% all predisposed cliques are preserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ibd2 nodes will receive equal vote from all nodes, because they are neutral
% ibd0 nodes will receive no vote, if included
% ibd1 nodes will receive vote from its own clique


% let deterministic nodes vote first
num = length(ibd12);
cluster_vote(1:num) = 0;

if( ~isempty(p_clique) )
    for i = 1:length(p_clique)
        temp = find(ibd12(p_clique(i),:) == 1);
        temp = intersect(temp, neighbor12);
        c_temp = setdiff(neighbor12, temp);
        if( ~isempty(temp) && ~isempty(c_temp) )
            cluster_vote(temp) = cluster_vote(temp) + 1;
            cluster_vote(c_temp) = cluster_vote(c_temp) - 1;
        end
    end
end

if( ~isempty(m_clique) )
    for i = 1:length(m_clique)
        temp = find(ibd12(m_clique(i),:) == 1);
        temp = intersect(temp, neighbor12);
        c_temp = setdiff(neighbor12, temp);
        if( ~isempty(temp) && ~isempty(c_temp) )
            cluster_vote(temp) = cluster_vote(temp) - 1;
            cluster_vote(c_temp) = cluster_vote(c_temp) + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let the remaining nodes vote

% order nodes according to degree in the maximum neighborhood
% add high degree node first to the vote
% if two cliques are not equal in size
% degree roughly reflects its membership
local_count(1:length(neighbor12)) = 0;
for i = 1:length(neighbor12)
    local_count(i) = nnz(ibd12(neighbor12(i),neighbor12));
end
[X,I] = sort(local_count, 'descend');
neighbor12 = neighbor12(I);

ibd12_temp = ibd12(neighbor12, neighbor12);

for i = 1:length(neighbor12)
    if( ~ismember(neighbor12(i), p_clique) && ~ismember(neighbor12(i), m_clique) )
        temp = find(ibd12(neighbor12(i),:) == 1);
        temp = intersect(temp, neighbor12);
        c_temp = setdiff(neighbor12, temp);
        if( ~isempty(temp) && ~isempty(c_temp) )
            % ibd0, ibd2 to index not voting
            mean1 = mean(cluster_vote(temp));
            mean2 = mean(cluster_vote(c_temp));
            if( mean1 >= mean2 )
                cluster_vote(temp) = cluster_vote(temp) + 1;
                cluster_vote(c_temp) = cluster_vote(c_temp) - 1;
            else
                cluster_vote(temp) = cluster_vote(temp) - 1;
                cluster_vote(c_temp) = cluster_vote(c_temp) + 1;
            end
        end
    end
end

[X,I] = sort(cluster_vote, 'descend');
ibd12_temp = ibd12(I,I);
% spy(ibd12_temp);
temp_vote = cluster_vote(neighbor12);
[X,I] = sort(temp_vote, 'descend');
ibd12_temp = ibd12(neighbor12(I),neighbor12(I));
% spy(ibd12_temp);
neighbor12 = neighbor12(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% partition according to vote mean bias
pivot = mean(cluster_vote(neighbor12));
% random sampling will yield a mean at pivot
% mean1 == mean2 == pivot will also happen, if ibd2/ibd0 to index
output_pm_clique = [];
vote_both = [];
vote_neither = [];

for i = 1:length(neighbor12)
    if( ~ismember(neighbor12(i), p_clique) && ~ismember(neighbor12(i), m_clique) ) 
        temp = find(ibd12(neighbor12(i),:) == 1);
        temp = intersect(temp, neighbor12);
        c_temp = setdiff(neighbor12, temp);
        if( ~isempty(temp) && ~isempty(c_temp) )
            mean1 = mean(cluster_vote(temp));
            mean2 = mean(cluster_vote(c_temp));            
            if( mean1 >= mean2 ) % or mean1 > pivot
                p_clique = [p_clique, neighbor12(i)];
            end
            if( mean1 < mean2 )
                m_clique = [m_clique, neighbor12(i)];
            end
        end
        if( isempty(c_temp) )
            vote_both = [vote_both, neighbor12(i)];
        end
        if( isempty(temp) )
            % this will never happen in the neighborhood of index
            % any node is at least connected to index
            vote_neither = [vote_neither, neighbor12(i)];
        end
    end
end



[X,I] = sort(cluster_vote(p_clique), 'descend');
output_p_clique = p_clique(I);
[X,I] = sort(cluster_vote(m_clique), 'descend');
output_m_clique = m_clique(I);
% spy(ibd12([p_clique,m_clique],[p_clique,m_clique]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end