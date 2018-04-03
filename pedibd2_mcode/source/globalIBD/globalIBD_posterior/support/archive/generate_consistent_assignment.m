function [ stat_interval stat_assignment_all ] = generate_consistent_assignment( output_dec_state, output_dec_state_sc, ia, range )




n_ind = ceil(sqrt(2*length(output_dec_state)));
n_ind = length(range.family_range);

assigned_list = 1:n_ind;

% generate global map for selected markers

pair_result = output_dec_state;
single_result = output_dec_state_sc;

% pair_result = random_ref_pair{1};
% single_result = random_ref_single{1};

stat_interval =  generate_consistent_interval( assigned_list, n_ind, pair_result, single_result );
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

[stat_assignment_all, stat_assignment_result] = build_global_inheritance(stat_interval(:,3), pair_result, single_result, map);
stat_interval(:,4) = stat_assignment_result;



end

