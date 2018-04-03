

function [data_re, list_re, error] = reorder_genotype(all_data, marker_list)

error = 0;
data_re = [];
list_re = [];
unordered = 0;

[error unordered] = check_marker(marker_list);
if( error == 1 )
    return;
end

[error] = check_genotype(all_data);
if( error == 1 )
    return;
end

[nLOC, nFIELDS] = size(marker_list);
if( nFIELDS < 2 || nFIELDS > 3 || nLOC <= 0 )
    error = 1;
    disp('error in reading marker data');
    return;
end

temp = all_data(:,7:end);
[nIND, len2] = size(temp);
if( nIND <= 0 )
    error = 1;
    disp('error in reading genotype data');
    return
end
if( nLOC*2 ~= len2 )
    error = 1;
    disp('genotype data and marker list do not match');
    return;
end

header = all_data(1:nIND,1:6);
genotype = int8(cat(3, temp(1:nIND,1:2:len2)', temp(1:nIND,2:2:len2)'));

if( unordered )

    [~, index] = sort(marker_list(:,2), 'ascend');

    list_re = marker_list(index,1:nFIELDS);
    
    genotype = genotype(index, 1:nIND, 1:2);

else 
    list_re = marker_list;
end

data_re.header = header;
data_re.genotype = genotype;

end



