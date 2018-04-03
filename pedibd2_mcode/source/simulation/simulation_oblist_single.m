function [ random_ibd_result, random_ibs_result ] = simulation_oblist_single( s_id_sc, input_family_range, pedigree_all_missing, random_ibd, random_genotype )

% transform to inner id
sp = pedigree_all_missing(input_family_range(s_id_sc),2);

random_ibs_result = genotype_sc(random_genotype(sp, 1:end));

random_ibd_result = genotype_sc(random_ibd(sp, 1:end));


end

