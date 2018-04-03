function [family_range error] = generate_families(input_pedigree)

    error = 0;
    family_range = [];
    
    pedigree = input_pedigree.structure;

    % check whether pedigree/family ids are correct
    error = check_family_id(pedigree);   
    if( error ~= 0 )
        disp('error in processing pedigrees');
        return;
    end      
    f_id = unique(pedigree(:,1));
    nfam = length(f_id);
    if( nfam <= 0 )
        error = 1;
        disp('error in processing family structures');
        return;
    end

    family_range = cell(nfam,1);
    
    % generate input pairs for each family  
    for f = 1:nfam
        [family_range{f} error] = extract1family(pedigree, f_id(f));
        if( error ~= 0 )
            disp('error in processing pedigrees');
            return;
        end
        family_range{f}.coding = input_pedigree.coding;
    end
    
    [error] = check_family_structure(pedigree, family_range);
    if( error ~= 0 )
        disp('error in processing family structures');
        return;
    end

    for f = 1:nfam
        [family_range{f}.allele_source.relevance, error] = allele_source1family(family_range{f}.structure); 
        if( error ~= 0 )
            disp('error in processing family structures');
            return;
        end
    end
    
    for f = 1:nfam
        [family_range{f}.founder_source, error] = founder_source1family(family_range{f}.structure);
        if( error ~= 0 )
            disp('error in processing family structures');
            return;
        end
    end


end


