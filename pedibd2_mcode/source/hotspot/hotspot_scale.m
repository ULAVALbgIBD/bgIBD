function [newmap error] = hotspot_scale(map, hotspots)

error = 0;
newmap = [];

[nloc, c] = size(map);
if( nloc <= 1 || c ~= 2 )
    error = 1;
    disp('error in physical map');
    return;
end

phy_map = map(1:nloc,2);
intervals = (phy_map(2:nloc) - phy_map(1:nloc-1));
if( any(intervals < 0) )
    error = 1;
    disp('error in physical map');
    return;
end

if( isempty(hotspots) )
    newmap = map;
    return;
end

[error] = check_hotspots(hotspots);
if( error ~= 0 )
    disp('error in hotspots input');
    return;
end

[num, c] = size(hotspots);

if( num <= 0 || c ~= 2 )
    error = 1;
    disp('error in hotspots list');
    return;
end



hot(1:nloc-1) = false;
hot_count(1:nloc-1) = 0;
cast_hot(1:num) = false;

newmap = zeros(nloc,2);
newmap(1:nloc,1) = map(1:nloc,1);

% overlay with hotspots
% mark all intervals overlapping with, containing or contained in hotspots
for i = 1:num
    spotL = hotspots(i,1);
    spotR = hotspots(i,2);
    left = (phy_map(1:nloc-1) <= spotR);
    right = (phy_map(2:nloc) >= spotL);
    hot( left & right ) = true;
    % how many hotspots cast effect in the interval
    hot_count( left & right ) = hot_count( left & right ) + 1/nnz( left & right );
    % do not count in hotspots out of range of current snp array
    if( any( left & right ) )
        cast_hot(i) = true;
    end
end

% enlarge hot area to 80% 
% shrink cold area to 20%

hot_area = sum(intervals(hot));
cold_area = sum(intervals(~hot));
total_area = hot_area + cold_area;
if( hot_area > 0 && cold_area > 0 && nnz(hot) > 0 && hot_area < 0.8 * total_area )
    increment = (0.8*total_area - hot_area)/nnz(hot);
    shrink = 0.2*total_area/cold_area;
else
    newmap = map;
    return;
end

% rescale intervals
% rounding

newintervals(hot) = intervals(hot) + increment;
newintervals(~hot) = intervals(~hot) * shrink;

newintervals = ceil(newintervals);

newmap(1,2) = map(1,2);
for i = 2:nloc
    newmap(i,2) = newmap(i-1,2) + newintervals(i-1);
end


% plot snp intervals overlapping with hotspots
% figure;
% hold on;
% for i = 1:nloc-1
%     line([phy_map(i),phy_map(i+1)], [hot(i),hot(i)], 'color', 'b', 'linewidth', 2);
% end
% ylim([-0.5,1.5]);
% for i = 1:num
%     line(hotspots(i,1:2), [0.9,0.9], 'color', 'r', 'linewidth', 2);
% end




end













