function [output_family error] = extract1family(pedigree_structure, family_id)

    error = 0;
    
    [family_range, pedigree_range_full, child_range error] = generate_list_1family(pedigree_structure, family_id);
    if( error ~= 0 )
        disp('error in parsing family structures');
        return;
    end
    

    [pairs reverse_pairs] = generate_pairlist(family_range);

    structure = pedigree_structure(pedigree_range_full,:);
    
    family_id = structure(1,12);
    
    output_family.family_range = family_range;
    output_family.pedigree_range_full = pedigree_range_full;
    output_family.child_range = child_range;
    output_family.pairs = pairs;
    output_family.reverse_pairs = reverse_pairs;
    output_family.family_id = family_id;
    output_family.structure = structure;

end