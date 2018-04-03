
% only works for large cliques, small cliques will fail

function [final_p_clique, final_m_clique, final_pm_clique] = reorganize_2chr(index, neighbor12, ibd12, p_clique_forced, p_clique_expanded, m_clique_forced, m_clique_expanded, vote_both, posterior, tolerance)


% make a list of individual can not share one allele IBD
% if all individuals are related, does it mean they can share one single
% allele?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do not add pm_clique to the pruning process,
% pm_clique is not informative, might disrupt the tolerance level

% pruning according to evidence total vote
% remove nodes which is not fully connected to the clique

final_pm_clique = [];
lenpm = double(length(vote_both));
% ibd2 clique
final_p_clique = [];
lenp = double(length(p_clique_expanded));
final_m_clique = [];
lenm = double(length(m_clique_expanded));


% tolerance for false negative, in most cases false negative is 0
% most errors are false positive
tolerance = tolerance;

% add consideration of one clique only
% if one clique is empty, all nodes will be added to pm_clique

for i = 1:length(ibd12)
    node = i;
    if( ~ismember(i, neighbor12) )
        if( nnz(ibd12(node, neighbor12)) < 10 )
            continue;
        end
        if( posterior(index, node, 1) > 0.5 )
            continue;
        end
        % only test non-neighborhood nodes, of very strong evidence
        
    end
    p_cohesion = nnz(ibd12(node,p_clique_expanded));
    m_cohesion = nnz(ibd12(node,m_clique_expanded));
    
    
    if( p_cohesion >= lenp * tolerance || m_cohesion >= lenm * tolerance )
        if( p_cohesion >= lenp * tolerance && m_cohesion >= lenm * tolerance && posterior(index, node, 3) > 0.5 )
            final_pm_clique = [node, final_pm_clique];
        else
            if( lenp == 0 && lenm == 0 )
                % cannot work in this situation
                % arbitrarily assign to p
                if( ismember(node, vote_both) )
                    % if vote_both cannot be identified as ibd2, then
                    % assign it to an arbirary clique
                    final_p_clique = [node, final_p_clique];
                end
            end
            if( lenp == 0 && lenm ~= 0 )
                if( m_cohesion >= lenm * tolerance )
                    final_m_clique = [node, final_m_clique];
                end
            end
            if( lenp ~= 0 && lenm == 0 )
                if( p_cohesion >= lenp * tolerance )
                    final_p_clique = [node, final_p_clique];
                end
            end
            if( lenp ~= 0 && lenm ~= 0 )
                if( ismember(node, p_clique_forced) )
                    if( p_cohesion >= lenp * tolerance )
                        final_p_clique = [node, final_p_clique];
                    else
                        if( posterior(index, node, 2) == 1 )
                            final_p_clique = [node, final_p_clique];
                            % parent child must be added
                        end
                    end
                end
                if( ismember(node, m_clique_forced) )
                    if( m_cohesion >= lenm * tolerance )
                        final_m_clique = [node, final_m_clique];
                    else
                        if( posterior(index, node, 2) == 1 )
                            final_m_clique = [node, final_m_clique];
                            % parent child must be added
                        end
                    end
                end
                if( ~ismember(node, p_clique_forced) && ~ismember(node, m_clique_forced) )
                    if( p_cohesion/lenp >= m_cohesion/lenm )
                        final_p_clique = [node, final_p_clique];
                    else
                        final_m_clique = [node, final_m_clique];
                    end
                end
            end
        end
    end
    if( p_cohesion <= lenp * tolerance && m_cohesion <= lenm * tolerance )
        % errors
        % remove these neighbors
    end
end

% add nodes, if even if it is not connected to vote_both

end

