
function [input pedigree expanded_genotype error] = process_input(pedigree_info, all_data, marker_list)

    global debug_mode;

    input = [];
    pedigree = [];
    expanded_genotype = [];
    parameters = [];
    error = 0;

    [input.family_range pedigree error] = process_pedigree(pedigree_info, 1);

    if( error == 1 )
        disp('error in pedigree file');
        return;
    end
    
    % column 1 ~ 6
    % family_id, individual_id, father_id, mother_id, sex, disease_status
    % (0 for missing)

    % column 7
    % if_genotyped (1 yes, 0 no)
    
    
    % column 1 ~ 6
    % family_id, individual_id, father_id, mother_id, sex, disease_status 
    % (0 for missing)
    % column 7 ~
    % bi-allelic genotype, each column one allele: 1 2 1 1 ....
    

    [expanded_genotype error] = pedigree2genotype(pedigree.structure, all_data);
    if( error == 1 )
        disp(['error in genotype file']);
        return;
    end
    
    
end


