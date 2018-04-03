
% obtain pairwise result from simulated families

function [ ref_ibd_pair, ref_ibd_sc, ibs, ibs_sc ] = simulation_all_pairs( input_family_range, input_pairlist, pedigree_all_missing, family_all )


for iteration = 1:length(family_all)
    
    random_pedigree = family_all{iteration};
    f = input_pairlist(1,3);
    
    for s_id = 1:length(input_pairlist)        

        [random_ibd_result, random_ibs_result] = simulation_oblist_pair(s_id, input_pairlist, pedigree_all_missing, random_pedigree{f,2}, random_pedigree{f,3}); 
        ref_ibd_pair{iteration}{s_id} = random_ibd_result + 1;
        % refers to ibd number 0, 1, 2
        ibs{iteration}{s_id} = random_ibs_result;
    end

    for s_id_sc = 1:length(input_family_range)

        [random_ibd_result_sc, random_ibs_result_sc] = simulation_oblist_single( s_id_sc, input_family_range, pedigree_all_missing, random_pedigree{f,2}, random_pedigree{f,3}); 
        ref_ibd_sc{iteration}{s_id_sc} = random_ibd_result_sc + 1;               
        ibs_sc{iteration}{s_id_sc} = random_ibs_result_sc;
    end

end

end

