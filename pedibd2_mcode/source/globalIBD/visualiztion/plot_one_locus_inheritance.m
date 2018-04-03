
% assigned_list specifies which individual to be included in the global
% inheritance

function plot_one_locus_inheritance( assigned_list, input_family_range, pedigree_all_missing, stat_assignment )

expand_stat_assignment(1:length(input_family_range),1:2) = 0;

for i = 1:length(assigned_list)
    expand_stat_assignment(assigned_list(i),1:2) = stat_assignment(i,1:2);
end

plot_relationship(input_family_range, pedigree_all_missing, expand_stat_assignment);

end

