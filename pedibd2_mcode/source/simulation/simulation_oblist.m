

function [random_ibd_result random_ibs_result] = simulation_oblist(pedigree_all_missing, s_id, input_pairlist, sm, pm)


% generate genotypes for a specific pair of individuals

f = input_pairlist(s_id,3);

%input parameters
%f
[family, random_ibd, random_genotype] =  simulation_oblist_family(f, pedigree_all_missing, sm, pm);
 

[random_ibd_result, random_ibs_result] = simulation_oblist_pair(s_id, input_pairlist, pedigree_all_missing, random_ibd, random_genotype);


end

