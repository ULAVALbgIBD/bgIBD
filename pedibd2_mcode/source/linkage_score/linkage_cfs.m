
%%



%%

stat_interval =  generate_consistent_interval( input_family_range, output_dec_state, output_dec_state_sc );

%%

% generate global inheritance only for sampled interval midpoint

[stat_assignment_all, stat_assignment_result] = build_global_inheritance(stat_interval(:,3), input_family_range, output_dec_state, output_dec_state_sc);

stat_interval(:,4) = stat_assignment_result;


%%

stat_score = [];
stat_score(1:length(stat_interval(:,1))) = 0;

for l = 1:length(stat_interval(:,1))
    stat_assignment = stat_assignment_all{l};
    allele_times(1:2*length(input_family_range)) = 0;
    if( stat_interval(l,4) == 1 )
        stat_score(l) = 1;
        for i = 1:length(assigned_list)
            if( pedigree_all_missing(input_family_range(assigned_list(i)),6) == 2 )
                for j = 1:2
                    allele_times(stat_assignment(i,j)) = allele_times(stat_assignment(i,j)) + 1;
                end
            end
        end
        for k = 1:length(allele_times)
            if( allele_times(k) ~= 0 )
                stat_score(l) = stat_score(l) * factorial(allele_times(k));
            end
        end
    end
end

stat_interval(:,5) = stat_score(1:end);

clear l allele_times i j k stat_assignment stat_score;

%%

% score permutation

stat_assignment_perm = [];


f = 43;
count = 0;

temp = find(pedigree_all_missing(:,1) == f);
[temp1, temp2] = sort(pedigree_all_missing(temp,9));

order = temp(temp2);

temp_assignment(1:length(order),1:2) = 0;

for t = 1:100000

    for i = 1:length(order)
        ind = order(i);
        inner = temp2(i);
        if( pedigree_all_missing(ind, 3) == 0 )
            temp_assignment(inner,1) = inner;
        else
            father = pedigree_all_missing(ind, 3);
            if( rand() > 0.5 )
                temp_assignment(inner,1) = temp_assignment(father,1);
            else
                temp_assignment(inner,1) = temp_assignment(father,2);
            end
        end
        if( pedigree_all_missing(ind, 4) == 0 )
            temp_assignment(inner,2) = -inner;
        else
            mother = pedigree_all_missing(ind, 4);
            if( rand() > 0.5 )
                temp_assignment(inner,2) = temp_assignment(mother,1);
            else
                temp_assignment(inner,2) = temp_assignment(mother,2);
            end
        end
    end

    stat_assignment_perm = temp_assignment(pedigree_all_missing(input_family_range(assigned_list),2), 1:2);

    stat_assignment = stat_assignment_perm + length(order) + 1;
    allele_times(1:2*length(order)) = 0;
    
    
    stat_score(t) = 1;
    for i = 1:length(assigned_list)
        if( pedigree_all_missing(input_family_range(assigned_list(i)),6) == 2 )
            for j = 1:2
                allele_times(stat_assignment(i,j)) = allele_times(stat_assignment(i,j)) + 1;
            end
        end
    end
    for k = 1:length(allele_times)
        if( allele_times(k) ~= 0 )
            stat_score(t) = stat_score(t) * factorial(allele_times(k));
        end
    end
    
end
% plot_relationship(input_family_range(assigned_list), pedigree_all_missing, stat_assignment_perm);


stat_score_perm = sort(stat_score);

clear ind inner i f  temp temp1 temp2 father mother temp_assignment stat_assignment_perm stat_assignment k t stat_score;

%%

temp = [];
for i = 1:length(stat_score_perm)
    temp = union(temp, stat_score_perm(i));
end
length(temp)

clear temp i;

%%

file_head = 'D:/IBD/script/data/CFS/results';

plot_multi_locus_inheritance( file_head, stat_interval, sampled_markerlist, assigned_list, input_family_range, pedigree_all_missing, stat_assignment_all );

clear file_head;

%%













