    
function [input, output, pedigree, parameters, error] ...
    = dbFamily_secure(pedigree_info, ...
    genotype_data, ...
    marker_list, ...
    hotspots, option)

    global echo;
    echo = 0;
    global debug_mode;
    debug_mode = 1; 
    global release;
    
    if( ~isempty(release) && release == 1 )
        debug_mode = 0;
        echo = 0;
    else
        debug_mode = 1;
        echo = 0;
    end   
    
    input = [];
    output = [];
    pedigree = [];
    parameters = [];
    error = 0;
    
    nIn = nargin;
    if( nIn ~= 5 )
        error = 1;
        disp('input arguments wrong numbers');
        return;
    end
    if( isempty(option) || length(option) < 3 )
        error = 1;
        disp('internal error');
        return;
    end
    
    error = check_numeric_array(pedigree_info, genotype_data, marker_list, hotspots);
    if( error ~= 0 )
        disp('input data format wrong');
        error = 1;
        return;
    end
    
    if( error == 0 )
        [genotype_data, marker_list, dense, error] = check_input(pedigree_info, genotype_data, marker_list);
    end    
    if( error ~= 0 )
        return;
    else
        if( dense )
            option(4) = 1;
        else
            option(4) = 0;
        end
    end
    
    
    if( error == 0 )        
        map.physical_map = marker_list(:,1:2);
        if( size(marker_list,2) == 3 )
            map.genetic_map = marker_list(:,3);
            disp('genetic distance specified');
            disp(' ');
        else
            if( isempty(hotspots) )        
                disp('using a flat genetic map 1cM/1Mb');
                disp(' ');
            end
        end
        [map.scaled_map, error] = hotspot_scale(map.physical_map, hotspots);
    end
    if( error ~= 0 )
        return;
    end
    
    
    if( error == 0 )    
        [input, output, pedigree, parameters, error] = ...
            dbFamily(pedigree_info, ...
            genotype_data, ...
            map, option);
    end
    
    if( error ~= 0 )
        return;
    end
    
    clear global echo;
    clear global debug_mode;
    
end



