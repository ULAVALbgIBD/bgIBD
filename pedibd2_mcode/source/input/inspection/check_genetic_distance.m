
function [dense_marker, error] = check_genetic_distance(list)

global echo;
error = 0;
dense_marker = false;

if( isempty(list) )
    error = 1;
    return;
end

[nLOC, ~] = size(list);

if( nLOC <= 1 )
    disp('<2 markers found');
    return;
end

if( any( list(2:nLOC,3) < list(1:nLOC-1,3) ) )
    disp('error in marker ordering');
    error = 1;
    return;    
end


average_marker_interval = (list(nLOC,3)-list(1,3))/(nLOC-1);
if( echo == 1 )
    disp(['Average genetic distance between markers: ', num2str(average_marker_interval, '%.2f'), ' cM']);
end


% 0.1cM is the threshold for dense marker check
if( average_marker_interval <= 0.1 )
    dense_marker = true;
else
    dense_marker = false;
    disp('warning: low marker density, average interval > 0.1cM');
end



end