function [smoothed_vit error] = suppress_ld(oblist, viterbi, posterior, map)


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
    cutoff = 2 * 10^6;
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
           for j = a:-1:1
               if( map(a) - map(j) > nearby )
                   left = j;
                   break;
               end
           end
           for j = b:1:nmarkers
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
           average_vit = round(mean(foreground(left:right)));
           if( a > 1 && b < nmarkers )
               if( smoothed_vit(a,1) == 2 && smoothed_vit(a,1) > average_vit )
                   smoothed_vit(a:b,1) = average_vit;
                   change = true;
               end            
           else
               if( smoothed_vit(a,1) == 2 && smoothed_vit(a,1) > average_vit && map(b) - map(a) < 0.5 * 10^6 )
                   smoothed_vit(a:b,1) = average_vit;
                   change = true;
               end           
           end
        end
    end
end


segments =  generate_consistent_interval_viterbi( smoothed_vit(1:nmarkers, 1:2) );
[nseg, c] = size(segments);
if( nseg <= 0 || c ~= 3 )
    error = 1;
    disp('error in segmentation');
    return;
end
foreground(1:nmarkers) = smoothed_vit(1:nmarkers, 1);
% very long regions, check for ibd0 spots, reassign if exceed threshold
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
        nearby = map(b) - map(a);
        % translate nearby distance from bp to markers
        % according to cutoff
        left = 1;
        right = nmarkers;
        for j = a:-1:1
           if( map(a) - map(j) > nearby )
               left = j;
               break;
           end
        end
        for j = b:1:nmarkers
           if( map(j) - map(b) > nearby )
               right = j;
               break;
           end
        end
        average_vit = round(mean(foreground(left:right)));
        if( left < a )
            average1 = round(mean(foreground(left:a-1)));
            if( average1 == 3 )
                continue;
            end
        end
        if( right > b )
            average2 = round(mean(foreground(b+1:right)));
            if( average2 == 3 )
                continue;
            end
        end
        % only consider ibd regions poping up from neighboring ones
        if( smoothed_vit(a) < average_vit  )
            continue;
        end

        ibd = sum(posterior(a:b,1:3),1);
        [~, index] = max(ibd);
        % use posterior max decoding, to replace viterbi, first
        if( a > 1 && b < nmarkers && map(b) - map(a) <= 5 * 10^6 )
            % not chromosome tips
            smoothed_vit(a:b, 1) = index;
        end
        % apply to suppress ibd 1 only
        if( smoothed_vit(a,1) == 2 )
            if( a > 1 && b < nmarkers )
                % not on chromosome tips
                % ibs 0 threshold
                if( map(b) - map(a) < 3 * 10^6 )
                    if( nnz(oblist(a:b)==1) > round(0.001 * (b-a+1)) )
                        smoothed_vit(a:b, 1) = 1;
                    end
                end
                if( map(b) - map(a) >= 3 * 10^6 && map(b) - map(a) <= 5 * 10^6 )
                    if( nnz(oblist(a:b)==1) > round(0.002 * (b-a+1)) )
                        smoothed_vit(a:b, 1) = 1;
                    end
                end
                if( map(b) - map(a) > 5 * 10^6 )
                    if( nnz(oblist(a:b)==1) > round(0.005 * (b-a+1)) )
                        smoothed_vit(a:b, 1) = 1;
                    end                    
                end
            else
                if( nnz(oblist(a:b)==1) > round(0.008 * (b-a+1)) )
                    smoothed_vit(a:b, 1) = 1;
                end
            end
        end
    end    
end



end
