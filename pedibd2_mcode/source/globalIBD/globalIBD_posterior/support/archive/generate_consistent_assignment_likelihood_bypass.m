function [ result ] = generate_consistent_assignment_likelihood_bypass( likelihood, ia, range, pedigree )

% quick method for simple families, not complete enumeration of all
% tranmission patterns


n_ind = ceil(sqrt(2*length(likelihood)));
n_ind = length(range.family_range);

assigned_list = 1:n_ind;

% generate global map for selected markers

stat_interval =  generate_consistent_interval_likelihood( likelihood );
stat_interval(:,4:5) = 0; %reserved for global map status and linkage score

map.list = assigned_list;

% the immediate neighborhood method works only when assigned_list ==
% family_range, otherwise regenerate the immediate list again

map.list_range = range.family_range(map.list);
map.reverse_list(1:max(range.family_range)) = 0; %not included families members maps to 0
for i = 1:length(assigned_list)
    map.reverse_list(range.family_range(assigned_list(i))) = i;
end

map.reverse_pair(1:length(assigned_list),1:length(assigned_list)) = 0;

for i = 1:length(assigned_list)
    for j = 1:length(assigned_list)
        map.reverse_pair(i,j) = range.reverse_pair(assigned_list(i),assigned_list(j));
    end
end


for i = 1:length(range.pedigree_range_full)
    map.immediate_ancestor{i,1} = ia{range.pedigree_range_full(i),1};
    map.immediate_ancestor{i,2} = ia{range.pedigree_range_full(i),2};
end

[stat_assignment_all, stat_assignment_result] = build_global_inheritance_likelihood_bypass(stat_interval(:,1:3), likelihood, map);
stat_interval(:,4:4+length(stat_assignment_result(1,:))-1) = stat_assignment_result;

input_family_range = range.family_range_outer;

stat_interval(:,4+length(stat_assignment_result(1,:))) = non_paramatric_linkage(stat_interval, stat_assignment_all, input_family_range, pedigree);

result.intervals = stat_interval;
result.alleles = stat_assignment_all;

end

