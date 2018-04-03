function [error] = check_kinship1family(family, kinship)


    error = 0;
    if( isempty(family) )
        error = 1;
        return;
    end
    family_range = family.family_range;
    nind = length(family_range);
    if( nind < 2 && ~isempty(kinship) )
        error = 1;
        return;
    end
    if( nind >= 2 )
        [r, c] = size(kinship);
        if( r ~= nind || c ~= nind )
            error = 1;
            return;
        end
    end
    pairs = family.pairs;
    if( nind < 2 )
        if( ~isempty(pairs) )
            error = 1;
            return;
        end
    else
        [num, cols] = size(pairs);
        if( num <= 0 || cols ~= 2 )
            error = 1;
            return;
        end
    end
    

end






