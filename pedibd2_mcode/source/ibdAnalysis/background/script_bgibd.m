%%


bg_epl_merge = dc_condensed_generate_emission_ibs_all_loci(sampled_af, pm_missing_rate, pm_typing_error, sampled_markerlist);

%%

bg_ibd_sc_all(1:2, 1:length(pm_member_list)) = 0;

for i = 1:length(pm_member_list)
    bg_ibd_sc_all(:,i) = ibd_expected_sc(all_data(pm_member_list(i), 7:end), input_epl_s);
end

clear i;

%%

bg_ibd_all(1:3, 1:length(pm_sibling_list)) = 0;

for i = 1:1:length(pm_sibling_list)
    ibs_data = ibs(all_data(pm_sibling_list(i,1), 7:end), all_data(pm_sibling_list(i,2), 7:end));
    bg_ibd_all(:,i) = ibd_expected(ibs_data, bg_epl_merge);
    
end
clear i ibs_data;

%%

bg_pi(1:1:length(pm_sibling_list)) = 0;
for i = 1:1:length(pm_sibling_list)
    bg_pi(i) = bg_ibd_all(2,i)/2 + bg_ibd_all(3,i);
end

clear i;
%%
