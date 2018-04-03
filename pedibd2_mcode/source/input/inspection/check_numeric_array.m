function error = check_numeric_array(pedigree_info, all_data, marker_list, hotspots)

error = 0;

if( ~isnumeric(pedigree_info) )
    error = 1;
    disp('error importing pedigree file: all fields must be integers');
    return;
end

if( ~isnumeric(all_data) )
    error = 1;
    disp('error importing genotype file: all fields must be integers');
    return;
end

if( ~isnumeric(marker_list) )
    error = 1;
    disp('error importing marker list: all fields must be numeric');
    return;
end

if( ~isempty(hotspots) )
    if( ~isnumeric(hotspots) )
        error = 1;
        disp('error importing hotspots');
        return;
    end
end

end