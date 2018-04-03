function [founder_list sibling_list error] = generate_list_allfamilies(all_families)
    
    error = 0;
    
    founder_list = [];
    sibling_list = [];
    fcount = 0;
    scount = 0;
    
    if( isempty(all_families) )
        error = 1;
        disp('no family structure found');
        return;
    end
    
    
    for f = 1:length(all_families)
        %temp is part of input_structure
        family = all_families{f}.structure;
        global_order = all_families{f}.pedigree_range_full;
        nind = length(global_order);
        
        if( nind <= 0 || any(global_order <= 0) )
            error = 1;
            disp('error in family partitioning');
            return;
        end
        [r, c] = size(family);
        if( r <= 0 || r ~= nind || c < 12 )
            error = 1;
            disp('error in family structures');
            return;
        end
        
        for j = 1:nind
            if( family(j,7) ~= 0 )
                if( family(j,3) == 0 )
                    if( family(j,4) == 0 )
                        fcount = fcount + 1;
                        founder_list(fcount) = global_order(j);
                    end
                end
            end
        end
        for j1 = 1:nind
            for j2 = j1+1:nind
                if( family(j1,7) ~= 0 && family(j2,7) ~= 0 )
                    if( ( family(j1,3) == family(j2,3) && family(j1,3) ~= 0 ) || ( family(j1,4) == family(j2,4) && family(j1,4) ~= 0 ) )
                        scount = scount + 1;
                        sibling_list(scount,1) = global_order(j1);
                        sibling_list(scount,2) = global_order(j2); 
                    end
                end
            end
        end
        
    end
    
end








