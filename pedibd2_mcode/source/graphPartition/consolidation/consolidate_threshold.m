function [T, error] = consolidate_threshold(T, dist, cutoff)

error = 0;

% use maximum distance

num = length(T);
if( num <= 0 )
    error = 1;
    return;
end

[d1 d2] = size(dist);
if( d1 <= 0 || d1 ~= num || d2 ~= d1 )
    error = 1;
    return;
end

if( cutoff < 0 )
    error = 1;
    return;
end

if( any(T <= 0 | T > num) )
    error = 1;
    disp('invalid cluster labeling');
    return;
end





while(true)
    mindist = inf;
    imin = 0;
    jmin = 0;
    for i = 1:num
        groupi = (T == i);
        if( ~any(groupi) )
            continue;
        end
        for j = i+1:num        
            groupj = (T == j);
            if( ~any(groupj) )
                continue;
            end
            Tdist = max(max(dist(groupi,groupj)));
            if( Tdist < mindist )
                mindist = Tdist;
                imin = i;
                jmin = j;
            end
        end
    end
    if( imin ~= 0 && jmin ~= 0 && mindist <= cutoff )
        if( imin < jmin )
            T( T == jmin ) = imin;
        else
            T( T == imin ) = jmin;
        end
    else
        break;
    end
end




end
