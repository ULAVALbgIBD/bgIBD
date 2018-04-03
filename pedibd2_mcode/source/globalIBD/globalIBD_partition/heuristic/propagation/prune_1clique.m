
% only works for large cliques, small cliques will fail

function [final_clique] = prune_1clique(pm_clique, clique, ibd12, posterior, tolerance)



final_clique = [pm_clique];
% must contain pm_clique;

neighbor = [pm_clique, clique];
len = double(length(neighbor));

tolerance = tolerance;


for i = 1:length(clique)
    node = clique(i);
    if( nnz(ibd12(node,neighbor)) >= (len - 1) * tolerance + 1 )
        temp = neighbor(ibd12(node,neighbor)>0);
        temp = setdiff(temp, node);
        if( ~isempty(temp) )
            temp_sum = sum(sum(posterior(node, temp, 2:3)));
            temp_mean = temp_sum/length(temp);
            if( temp_mean > 0 )
                final_clique = [final_clique, node];
            end
        end
    else
        if( sum(posterior(pm_clique, node, 2)) == length(pm_clique) )
            % parent-child relation must be added
            final_clique = [final_clique, node];
        end
    end
end


end

