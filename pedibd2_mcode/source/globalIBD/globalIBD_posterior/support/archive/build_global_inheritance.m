% build global relationship from pairwise relationship

function [ stat_assignment_all, stat_assignment_result ] = build_global_inheritance( markers, output_dec_state, output_dec_state_sc, map)

count = 0;


n_ind = length(map.list);

num_m = length(markers);
stat_relation(1:n_ind,1:n_ind) = 0;
stat_relation_prob(1:n_ind,1:n_ind) = 0;
stat_relation_sc(1:n_ind) = 0;

time = cputime;
for l = 1:num_m

    locus = markers(l);

    for i = 1:length(map.list)
        for j = i+1:length(map.list)
            stat_relation(i,j) = output_dec_state{map.reverse_pair(i,j)}(locus);
            stat_relation(j,i) = stat_relation(i,j);
            stat_relation_prob(i,j) = output_dec_state{map.reverse_pair(i,j)}(locus,2);
            stat_relation_prob(j,i) = stat_relation_prob(i,j);           
        end
        stat_relation_sc(i) = output_dec_state_sc{map.list(i)}(locus);
    end

    %start individual
    %end individual
    %number of alleles

    stat_assignment = [];
    stat_assignment(1:length(stat_relation),1:2) = 0;
    
%     [stat_assignment, success, complexity] = assign_allele(stat_assignment, 1, length(stat_assignment), stat_relation, stat_relation_sc, map, 0);
%     complexity
%     success
    
    % time likelihood together to optimize, branch and bound search
    success = 0;
    tolerance = 0;
    while( success == 0 )
        [stat_assignment, success, complexity] = assign_allele_nonrecur(stat_relation, stat_relation_sc, map, tolerance);
        tolerance = tolerance + 1;
    end
    l
    
    
%     stat_assignment = assign_allele_heuristic(stat_relation, stat_relation_sc);
    
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
%     temp_result = 0;

    

    stat_assignment_all{l} = stat_assignment;
    stat_assignment_result(l) = temp_result;

end

cputime - time
end

