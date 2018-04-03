function result = hotspot_overlap(rec, hotspots, map)


result(1:length(rec)) = 0;

max_length = max(hotspots(:,1));
move = round(max_length * rand(1,length(rec)));

%permutation over all snp locations, instead of all base pair positions

move2 = 1 + floor(length(map) * rand(1,length(rec)));
move3 = 1 + floor(length(map) * rand(1,length(hotspots)));

for i = 1:length(hotspots)
    temp = hotspots(i,2) - hotspots(i,1);
    hotspots(i,1) = map(move3(i));
    hotspots(i,2) = map(move3(i)) + temp;
end

for j = 1:length(rec)
    

    a = rec(j, 1); 
    b = rec(j, 2); 

    if( b < a )
        temp = b;
        b = a;
        a = temp;
    end

    %a = map(move2(j));
    %b = map(move2(j)) + rec(j, 2) - rec(j, 1);    
    
    min = inf;
    for k = 1:length(hotspots)
        if( a > hotspots(k, 1) )
            if( a - hotspots(k, 1) < min )
                min = a - hotspots(k, 1);
            end
            continue;
        end
        if( b < hotspots(k, 2) )
            if( hotspots(k, 2) - b < min )
                min = hotspots(k, 2) - b;
            end
            continue;
        end
        result(j) = 1;
        min = 0;
        break;               
    end
end



