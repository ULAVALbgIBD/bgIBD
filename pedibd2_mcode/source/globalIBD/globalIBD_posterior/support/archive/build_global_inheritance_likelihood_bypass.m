% build global relationship from pairwise relationship

function [ stat_assignment_all, stat_assignment_result ] = build_global_inheritance_likelihood_bypass( markers, likelihood, map)


n_ind = length(map.list);

num_m = length(markers(:,1));

stat_relation_prob(1:n_ind,1:n_ind,1:3) = 0;
stat_relation(1:n_ind,1:n_ind) = 0;
stat_assignment(1:n_ind,1:2) = 0;



time = cputime;
for l = 1:num_m

    m_l = markers(l,3);
    s_l = markers(l,1);
    t_l = markers(l,2);

    
    for i = 1:length(map.list)
        for j = i+1:length(map.list)
            for num = 1:3
                stat_relation_prob(i,j,num) = mean(likelihood{map.reverse_pair(i,j)}(s_l:t_l,num));
                stat_relation_prob(j,i,num) = stat_relation_prob(i,j,num);
            end       
            [temp stat_relation(i,j)] = max( stat_relation_prob(i,j,:) );
            stat_relation(j,i) = stat_relation(i,j);
        end
    end
    
    % time likelihood together to optimize, branch and bound search

    [stat_assignment, recombination, lk, diff, complexity] = assign_allele_nonrecur_likelihood_bypass(stat_relation, stat_relation_prob, map, stat_assignment);
    complexity
    recombination
    
    

    temp_result = 0;
    for i = 1:length(stat_relation)
        for j = i+1:length(stat_relation)
            temp = count_ibd( stat_assignment(i,1:2), stat_assignment(j,1:2) );
            if( temp + 1 ~= stat_relation(i,j) )
                temp_result = temp_result + 1;
%                 [i, j, l, temp, stat_relation(i,j)-1]
            end
        end
    end
    

    stat_assignment_all{l} = stat_assignment;
    stat_assignment_result(l,1:3) = [lk, temp_result, complexity];

end

cputime - time
end

