

function [smoothed_vit error] = resist_fluctuation(oblist, viterbi, posteriorIBD1, posterior, map)


smoothed_vit = [];
error = 0;

[nmarkers, c] = size(map);
if( nmarkers <= 0 || c ~= 1 )
    error = 1;
    disp('error in marker list');
    return;
end

[r, c] = size(viterbi);
if( r <= 0 || r ~= nmarkers || c ~= 2 )
    error = 1;
    disp('error in viterbi decoding');
    return;
end

segments =  generate_consistent_interval_viterbi( viterbi );
[nseg, c] = size(segments);
if( nseg <= 0 || nseg > nmarkers || c ~= 3 || segments(nseg,2) ~= nmarkers )
    error = 1;
    disp('error in segmenting chromosomes');
    return;
end

[d1, d2, d3] = size(posteriorIBD1);
if( d1 <= 0 || d1 ~= nmarkers || d2 ~= 2 || d3 ~= 2)
    error = 1;
    disp('error in posterior decoding');
    return;
end

[r, c] = size(posterior);
if( r <= 0 || r ~= nmarkers || c ~= 3 )
    error = 1;
    disp('error in posterior decoding');
    return;
end

len = length(oblist);
if( len <= 0 || len ~= nmarkers )
    error = 1;
    disp('error in observation list');
    return;
end

smoothed_vit(1:nmarkers, 1:2) = viterbi(1:nmarkers, 1:2);


segments =  generate_consistent_interval_viterbi( smoothed_vit );
[nseg, c] = size(segments);
if( nseg <= 0 || c ~= 3 )
    error = 1;
    disp('error in segmentation');
    return;
end
foreground = smoothed_vit(1:nmarkers, 1);
% long range ibd 0 region, mis-suppressed
cutoff = 10 * 10^6;
for i = 1:nseg
    a = segments(i,1);
    b = segments(i,2);
    if( ~all(foreground(a:b) == foreground(a)) )
        error = 1;
        disp('error in segmenting chromosomes');
        return;
    end
    if( a <= 0 || b <= 0 || a > nmarkers || b > nmarkers || a > b )
        error = 1;
        disp('error in segmenting chromosomes');
        return;
    end
    
    if( map(b) - map(a) < cutoff )
        % look around a neighboring area of same length
        nearby = (map(b) - map(a));
        left = 1;   % to left side marker of distance nearby bps
        right = nmarkers;   % to right side marker of distance nearby bps
        close_left = a;     % marker preserve the same ibd 1 status to left
        % left and right boundaries with ibd 1
        for j = a:-1:1
            if( foreground(j) == 2 )
                close_left = j;
            end
            if( map(a) - map(j) > nearby )
                left = j;
                break;
            end
        end
        close_right = b;    % marker preseve the same ibd 1 status to right
        for j = b:1:nmarkers
           if( foreground(j) == 2 )
               close_right = j;
           end
           if( map(j) - map(b) > nearby )
               right = j;
               break;
           end
        end
        
       if( (b - a + 1) > 500 )
           % bypass very confident dense pikes, many markers
           if( mean( posterior(a:b, smoothed_vit(a,1)) ) > 0.99 )
               continue;
           end
       end         
        
        % average ibd of nearby two portions, left and right
        average_vit = round(mean(foreground([left:a-1,b+1:right])));
        averageL = round(mean(foreground(left:a-1)));
        averageR = round(mean(foreground(b+1:right)));
        if( left < (a-1) && (b+1) < right && averageL ~= averageR )
           continue;
        end            
        if( average_vit == 2 )
            % the region is not on chromosome tips
            % or valid change boundaries
            if( close_left < a && close_right > b )
                m1 = mean(posteriorIBD1(close_left:a-1,1,1)) + mean(posteriorIBD1(close_left:a-1,1,2)) - mean(posteriorIBD1(close_left:a-1,2,1)) - mean(posteriorIBD1(close_left:a-1,2,2));
                m2 = mean(posteriorIBD1(b+1:close_right,1,1)) + mean(posteriorIBD1(b+1:close_right,1,2)) - mean(posteriorIBD1(b+1:close_right,2,1)) - mean(posteriorIBD1(b+1:close_right,2,2));
                if( abs(m1-m2) > 1 )
                    continue;
                end
            end
        end       
        if( smoothed_vit(a,1) < average_vit && smoothed_vit(a,1) == 1 )
            % the region is not on chromosome tips
            % or valid change boundaries
            if( close_left < a && close_right > b )
                if( nnz(oblist(a:b) == 1) <= 2 )
                   smoothed_vit(a:b,1) = average_vit;
                end
            end
        end
    end
end



end

