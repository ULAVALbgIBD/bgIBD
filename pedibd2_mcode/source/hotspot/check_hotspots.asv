function [error] = check_hotspots(hotspots)

error = 0;

if( isempty(hotspots) )
    return;
end

[num c] = size(hotspots);
if( num <= 0 || c ~= 2 )
    error = 1;
    disp('error in hotspots file format');
    return;
end

if( any(hotspots(1:num,2) < hotspots(1:num,1)) )
    error = 1;
    disp('error in hotspots file format');
    return;
end

if( any(any(hotspots(1:num,1:2) <= 0)) )
    error = 1;
    disp('error in hotspots file format');
    return;
end

approx_range = max(hotspots(1:num,2)) - min(hotspots(1:num,1)) + 1;
hot_range = sum(hotspots(1:num,2) - hotspots(1:num,1));

if( approx_range <= 0 )
    error = 1;
    disp('error in hotspots file format');
    return;
end



end