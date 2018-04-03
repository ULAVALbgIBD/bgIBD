
function [smoothed_vit error] = remove_pikes(oblist, viterbi, posteriorIBD1, posterior, map)


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
change = true;
while(change)
    change = false;
    foreground = smoothed_vit(1:nmarkers, 1);
    cutoff = 10^6;
    % very short LIFTED regions
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
           nearby = 2 * (map(b) - map(a));
           left = 1;
           right = nmarkers;
           close_left = a;
           for j = a:-1:1
               if( foreground(j) == 2 )
                   close_left = j;
               end
               if( map(a) - map(j) > nearby )
                   left = j;
                   break;
               end
           end
           close_right = b;
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
           average_vit = round(mean(foreground([left:a-1,b+1:right])));
           averageL = round(mean(foreground(left:a-1)));
           averageR = round(mean(foreground(b+1:right)));

           if( left < (a-1) && (b+1) < right && averageL ~= averageR && abs(averageL - smoothed_vit(a,1)) <= 1 && abs(averageR - smoothed_vit(a,1)) <= 1 )
               continue;
           end
           if( average_vit == 2 )
               % exclude ibd 1 maternal/paternal shift
               % popping up ibd 2 is valid
               if( close_left < a && close_right > b )
                   % 1 w/ 1 vs 2 w/ 2
                   % large interval 1, for difference maybe inferenced by
                   % noise at pike
                   % ideal difference is 2, 1 - 0 versus 0 - 1
                   m1 = mean(posteriorIBD1(close_left:a-1,1,1)) + mean(posteriorIBD1(close_left:a-1,1,2)) - mean(posteriorIBD1(close_left:a-1,2,1)) - mean(posteriorIBD1(close_left:a-1,2,2));
                   m2 = mean(posteriorIBD1(b+1:close_right,1,1)) + mean(posteriorIBD1(b+1:close_right,1,2)) - mean(posteriorIBD1(b+1:close_right,2,1)) - mean(posteriorIBD1(b+1:close_right,2,2));
                   if( abs(m1-m2) > 1 )
                       continue;
                   end
               end
           end 
           if( a > 1 && b < nmarkers )
               if( smoothed_vit(a,1) > average_vit )
                   smoothed_vit(a:b,1) = average_vit;
                   change = true;
               end             
           else
               if( smoothed_vit(a,1) > average_vit && map(b) - map(a) < 0.5 * 10^6 )
                   smoothed_vit(a:b,1) = average_vit;
                   change = true;
               end           
           end
        end
    end
end

segments =  generate_consistent_interval_viterbi( smoothed_vit );
[nseg, c] = size(segments);
if( nseg <= 0 || c ~= 3 )
    error = 1;
    disp('error in segmentation');
    return;
end
change = true;
while(change)
    change = false;
    foreground = smoothed_vit(1:nmarkers, 1);
    cutoff = 0.5 * 10^6;
    % extra very short DEPRESSED segments
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
           nearby = 2 * (map(b) - map(a));
           left = 1;
           right = nmarkers;
           close_left = a;
           for j = a:-1:1
               if( foreground(j) == 2 )
                   close_left = j;
               end
               if( map(a) - map(j) > nearby )
                   left = j;
                   break;
               end
           end
           close_right = b;
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
           % average does not include the region itself
           average_vit = round(mean(foreground([left:a-1,b+1:right])));
           averageL = round(mean(foreground(left:a-1)));
           averageR = round(mean(foreground(b+1:right)));
           if( left < (a-1) && (b+1) < right && averageL ~= averageR )
               continue;
           end           
           if( average_vit == 2 )
               if( close_left < a && close_right > b )
                   m1 = mean(posteriorIBD1(close_left:a-1,1,1)) + mean(posteriorIBD1(close_left:a-1,1,2)) - mean(posteriorIBD1(close_left:a-1,2,1)) - mean(posteriorIBD1(close_left:a-1,2,2));
                   m2 = mean(posteriorIBD1(b+1:close_right,1,1)) + mean(posteriorIBD1(b+1:close_right,1,2)) - mean(posteriorIBD1(b+1:close_right,2,1)) - mean(posteriorIBD1(b+1:close_right,2,2));
                   if( abs(m1-m2) > 1 )
                       continue;
                   end
               end
           end           
           if( smoothed_vit(a,1) == 1 && smoothed_vit(a,1) < average_vit && nnz(oblist(a:b) == 1) < 5 )
               smoothed_vit(a:b,1) = average_vit;
               change = true;
           end
           if( smoothed_vit(a,1) == 2 && smoothed_vit(a,1) < average_vit && nnz(oblist(a:b) <= 2) < 5 )
               smoothed_vit(a:b,1) = average_vit;
               change = true;
           end           
        end
    end
end

end
