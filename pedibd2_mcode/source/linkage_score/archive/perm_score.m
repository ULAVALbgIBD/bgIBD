
%be careful of two options below, additive and nonadditive

function stat_score_perm = perm_score(range, pedigree)

stat_assignment_perm = [];


family = pedigree(range.pedigree_range_full,:);
family_range = range.pedigree_range_full(range.family_range);

temp_assignment(1:length(family(:,1)),1:2) = 0;

stat_interval(:,1) = 1:100000;

% this assumes pedigree is partially ordered

for t = 1:length(stat_interval)
    
    for i = 1:length(family(:,1))
        if( family(i,3) == 0 )
            temp_assignment(i,1) = i;
        else
            father = family(i,3);
            if( rand() > 0.5 )
                temp_assignment(i,1) = temp_assignment(father,1);
            else
                temp_assignment(i,1) = temp_assignment(father,2);
            end
        end
        if( family(i,4) == 0 )
            temp_assignment(i,2) = -i;
        else
            mother = family(i,4);
            if( rand() > 0.5 )
                temp_assignment(i,2) = temp_assignment(mother,1);
            else
                temp_assignment(i,2) = temp_assignment(mother,2);
            end
        end
    end

    stat_assignment_perm{t} = temp_assignment(pedigree(family_range,2), 1:2);     
    
end
% plot_relationship(input_family_range(assigned_list), pedigree_all_missing, stat_assignment_perm);

stat_score = non_paramatric_linkage(stat_interval, stat_assignment_perm, family_range, pedigree);

% stat_score = non_paramatric_linkage_additive(stat_interval, stat_assignment_perm, family_range, pedigree);
stat_score_perm = sort(stat_score);


end


