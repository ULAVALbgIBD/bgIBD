function [error unordered] = check_marker(marker_list)
    
error = 0;
unordered = false;
[nrows ncols] = size(marker_list);
if( nrows <= 0 )
    error = 1;
    disp('empty marker list');
    return;
end
    
if( ncols ~= 2 && ncols ~= 3 )
    error = 1;
    disp('marker list data wrong format, require two or three columns');
    return;
end

error = check_integer(marker_list(:,1:2), [0,inf]);
if( error ~= 0 )
    disp('marker list data wrong format, all fields must be non-negative integers');
    return;
end

error = check_integer(marker_list(:,1), [0,2^15]);
if( error ~= 0 )
    disp('marker list data wrong format, invalid chromosome number');
    return;
end

chr = unique(marker_list(:,1));
if( length(chr) ~= 1 )
    error = 1;
    disp('marker list data wrong format, no chromosome or more than 1 chromosomes specified');
    return;
end


if( nrows > 1 && any( marker_list(2:nrows,2) < marker_list(1:nrows-1,2) ) ) 
    disp('warning: marker location is not in non-decreasing order');
    disp('assuming the same order in genotype data');
    disp('output result will be ordered non-decreasingly by marker');
    unordered = true;
end

if( nrows > 1 && ncols == 3 && any( marker_list(2:nrows,3) < marker_list(1:nrows-1,3) ) ) 
    disp('warning: marker location is not in non-decreasing order');
    disp('assuming the same order in genotype data');
    disp('output result will be ordered non-decreasingly by marker');
    unordered = true;
end
    
end