function [all_data_reorder, marker_list_reorder, dense_marker, error] = check_input(pedigree_info, all_data, marker_list)

error = 0;
all_data_reorder = [];
marker_list_reorder = [];
dense_marker = false;


[error] = check_pedigree(pedigree_info);
if( error ~= 0 )
    return;
end

[all_data_reorder, marker_list_reorder, error] = reorder_genotype(all_data, marker_list);
if( error ~= 0 )
    return;
end

[dense_marker, error] = check_marker_interval(marker_list_reorder);
if( error ~= 0 )
    return;
end

[~, nFIELDS] = size(marker_list_reorder);
if( nFIELDS == 3 )
    [dense_marker, error] = check_genetic_distance(marker_list_reorder);
    if( error ~= 0 )
        return;
    end
end

end