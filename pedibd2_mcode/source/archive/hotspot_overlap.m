function result = hotspot_overlap(rec, hotspots)

% temp = mean(hotspots(:,1:2), 2);
% hotspots(:,1) = temp;
% hotspots(:,2) = temp;
result(1:length(rec)) = 0;

max_length = max(hotspots(:,1));
move = round(max_length * rand(1,length(rec)));
for j = 1:length(rec)
    
    a = rec(j, 1);
    b = rec(j, 2);

    if( b < a )
        temp = b;
        b = a;
        a = temp;
    end
    
%     a = mean(rec(j,1:2)) - 0;
%     b = mean(rec(j,1:2)) + 0;
    
    


    min = inf;
    for k = 1:length(hotspots)
        if( a > hotspots(k, 2) )
            if( a - hotspots(k, 2) < abs(min) )
                min = hotspots(k, 1) - a;
            end
            continue;
        end
        if( b < hotspots(k, 1) )
            if( hotspots(k, 1) - b < abs(min) )
                min = hotspots(k, 1) - b;
            end
            continue;
        end
        min = 0;
        break;               
    end
    result(j) = min;
end



