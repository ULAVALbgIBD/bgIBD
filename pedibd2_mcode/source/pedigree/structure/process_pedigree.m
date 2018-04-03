
function [family_range pedigree error] = process_pedigree(reordered, split)

    error = 0;
    family_range = [];
    pedigree = [];
    
    [pedigree, error] = pedigree_structure(reordered);
    if( error ~= 0 )
        disp('error in pedigree structure');
        return;
    end
    
    % split families
    if( split == 1 )
        [pedigree.family_size pedigree.structure error] = family_size(pedigree.structure);
    end    
    if( error ~= 0 )
        disp('error in processing pedigrees');
        return;
    end     
    
    [family_range error] = generate_families(pedigree);    
    if( error ~= 0 )
        disp('error in family structures');
        return;
    end
    
    [error] = check_family_structure(pedigree.structure, family_range);
    if( error ~= 0 )
        disp('error in family structures');
        return;
    end
    
    % load information back to pedigree data
    % index are given in absolute pedigree order
    [pedigree.founder_list pedigree.sibling_list error] = generate_list_allfamilies(family_range);
    if( error ~= 0 )
        disp('error in processing family structures');
        return;
    end
    [pedigree.allele_source error] = allele_source(family_range);
    if( error ~= 0 )
        disp('error in processing family structures');
        return;
    end
    [pedigree.founder_source error] = founder_source(family_range);
    if( error ~= 0 )
        disp('error in processing family structures');
        return;
    end
    
end







