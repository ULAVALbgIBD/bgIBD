function [error] = check_posterior1family(family, posterior, markers)

error = 0;

if( isempty(family) )
    error = 1;
    return;
end

pairs = family.pairs;

if( isempty(pairs) )
    if( ~isempty(posterior) )
        error = 1;
        return;
    else
        return;
    end
end

[num, cols] = size(pairs);

if( num < 1 || cols ~= 2 )
    error = 1;
    return;
end



[nmarkers, cols] = size(markers);
if( nmarkers < 1 || cols ~= 2 )
    error = 1;
    return;
end

if( ndims(posterior) ~= 3 )
    error = 1;
    return;
end
[d1 d2 d3] = size(posterior);
if( d1 ~= num || d2 ~= nmarkers || d3 ~= 3 )
    error = 1;
    return;
end

end

