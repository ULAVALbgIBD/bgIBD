function error = output_kinship(fid, all_families, kinship)

    error = 0;
    nfam = length(all_families);
    if( nfam <= 0 )
        error = 1;
        return;
    end
    if( length(kinship) ~= nfam )
        error = 1;
        return;
    end
    
    for i = 1:nfam
        if( isempty(all_families{i}) )
            error = 1;
            return;
        end
        [error] = check_kinship1family(all_families{i}, kinship{i});
        if( error ~= 0 )
            return;
        end        
    end

    fprintf(fid, 'family,      id1    -    id2 \t\tkinship-coefficient\n');
    for f = 1:nfam
        family = all_families{f}.structure;
        family_range = all_families{f}.family_range;
        family_id = all_families{f}.family_id;
        nind = length(family_range);
        if( nind < 2 )
            continue;
        end      
        for i = 1:nind
            for j = 1:nind
                id1 = family(i,8);
                id2 = family(j,8);
                k = kinship{f}(i,j);
                fprintf(fid, '%6d:   %6d    - %6d \t\t', family_id, id1, id2);
                fprintf(fid, '%f\t', k);
                fprintf(fid, '\n');                
            end
        end        
    end
    
end







