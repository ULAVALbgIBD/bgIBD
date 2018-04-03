function [error] = output_posterior(fid, all_families, posterior, parameters)
    
    error = 0;
    
    nfam = length(all_families);
    if( nfam <= 0 )
        error = 1;
        return;
    end
    
    if( length(posterior) ~= nfam )
        error = 1;
        return;
    end
    
    sampled_markerlist = parameters.sampled_markerlist;
    [nloc, ~] = size(sampled_markerlist);
    
    for i = 1:nfam
    
        [error] = check_posterior1family(all_families{i}, posterior{i}, parameters.sampled_markerlist);
        if( error ~= 0 )
            return;
        end
    
    end
    
    fprintf(fid, '*****          ibd probability           *****\n\n\n\n');
    
    fprintf(fid, 'family, \t      id1    -     id2\t\t    \t\teach marker\n');
    for i = 1:length(all_families)
        pairs = all_families{i}.pairs;
        if( isempty(pairs) )
            continue;
        end
        [num, c] = size(pairs);
        if( num <= 0 || c ~= 2 )
            error = 1;
            return;
        end
        family = all_families{i}.structure;
        [nind, fields] = size(family);
        if( nind <= 0 || fields < 12 )
            error = 1;
            return;
        end
        if( any(any(pairs(:,1:2) > nind)) || any(any(pairs(:,1:2) < 1)) )
            error = 1;
            return;
        end
        family_id = all_families{i}.family_id;
        for j = 1:num
            id1 = family(pairs(j,1),8);
            id2 = family(pairs(j,2),8);
            temp = reshape(posterior{i}(j, 1:nloc, 1:3), [nloc,3]);    %flattened format
            for k = 1:3            
                fprintf(fid, '%6d:   \t   %6d    -  %6d\t\tIBD%d\t\t', family_id, id1, id2, k-1);
                fprintf(fid, '%f\t', temp(:,k));
                fprintf(fid, '\n');
            end
        end
    end
    
    
end












