
function [all_inheritance, ...
    kinship, kinship2, kinship2ex, ...
    singleLINEAL, ...
    error] ...
    = generate_inheritance_1family(family_range)

    error = 0;
    all_inheritance = [];
    kinship = [];
    kinship2 = [];
    kinship2ex = [];
    singleLINEAL = [];
    
    if( isempty(family_range) )
        error = 1;
        disp('error in family structures');
        return;
    end
    
    if( length(family_range.pedigree_range_full) < 1 )
        error = 1;
        disp('empty family');
        return;
    end
    
    if( length(family_range.family_range) < 2 )
        disp(['family ', num2str(family_range.family_id), ': < 2 genotyped individuals: kinship skipped']);
        return;
    end
    
    display(' ');    
    display(['family ',num2str((family_range.structure(1,12))), ': calculating kinship, total ', num2str((length(family_range.pedigree_range_full))^2), ' relative pairs ...']);
    time = cputime;
    
    
    [singleLINEAL, error] = generate_allpaths(family_range.structure);
    if( error ~= 0 )
        disp('error in generating kinship');
        return;
    end
    
    [all_inheritance, error] = descent_graph_2alleles_all(family_range);
    if( error ~= 0 )
        disp('error in generating kinship');
        return;
    end

    [kinship, kinship2, kinship2ex, error] = ...
        all_alleles_kinship(all_inheritance, family_range.family_range, family_range.pedigree_range_full);
    if( error ~= 0 )
        disp('error in generating kinship');
        return;
    end
    [error] = check_kinship1family(family_range, kinship);
    if( error ~= 0 )
        disp('error in generating kinship');
        return;
    end

    
    
    display(['kinship computing time: ', num2str(cputime - time), ' seconds']);    
        
end
        









