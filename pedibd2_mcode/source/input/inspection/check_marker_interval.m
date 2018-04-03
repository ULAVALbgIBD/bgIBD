
function [dense_marker, error] = check_marker_interval(list)

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

if( nLOC > 1 && any( list(2:nLOC,2) < list(1:nLOC-1,2) ) )
    disp('error in marker ordering');
    error = 1;
    return;    
end


average_marker_interval = (list(nLOC,2)-list(1,2))/(nLOC-1);
if( average_marker_interval < 100 )
    if( average_marker_interval < 1 )
        error = 1;
        disp(['average marker interval too small, ~', num2str(average_marker_interval, '%f'), ' b']);
        return;
    else
        disp(['average marker interval, ~', num2str(average_marker_interval, '%.1f'), ' b']);       
    end
else
    if( echo == 1 )
        disp(['average marker interval distance: ', num2str(average_marker_interval/1000, '%.2f'), ' kb']);
    end
end

% 100kb is the threshold for dense marker check
if( average_marker_interval <= 100 * 10^3 )
    dense_marker = true;
else
    dense_marker = false;
    disp('warning: low marker density, average interval > 100kb');
end


end