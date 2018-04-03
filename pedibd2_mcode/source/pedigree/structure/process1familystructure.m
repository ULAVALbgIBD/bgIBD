function [reorderF, error] = process1familystructure(family)

        error = 0;
        reorderF = [];

        [rows, cols] = size(family);
        if( rows <= 0 )
            error = 1;
            disp('empty family');
            return;
        end

        error = check_family(family);
        if( error ~= 0 )
            error = 1;
            disp('error in family structures');
            return;
        end
        
        f_id = unique(family(:,1));
        if( length(f_id) ~= 1 )
            error = 1;
            disp('error in processing families');
            return;
        end
        
        maxid = max(max(family(:,2:4)));        
        defined(1:maxid) = 0;
        assigned(1:maxid) = 0;
        sex(1:maxid) = 0;
        original_id(1:maxid) = 0;
        
        % maintain the partial order that parents proceeds children
        % add all seen ids to the list
        
        % sequential order
        nind = 0;
        for j = 1:rows
            id = family(j,2);
            defined(id) = j;
            for k = [2,3,4]
                id = family(j,k);
                if( id ~= 0 )
                    if( assigned(id) == 0 )
                        nind = nind + 1;
                        order = nind;
                        assigned(id) = order;
                        original_id(order) = id;
                        if( k == 2 )
                            sex(order) = family(j,5);
                        end
                        if( k == 3 )
                            sex(order) = 1;
                        end
                        if( k == 4 )
                            sex(order) = 2;
                        end
                    else
                        order = assigned(id);
                        if( k == 2 )
                            if( sex(order) ~= family(j,5) )
                                error = 1;
                            end
                        end
                        if( k == 3 )
                            if( sex(order) == 0 )
                                sex(order) = 1;
                            end
                            if( sex(order) ~= 1 )
                                error = 1;
                            end
                        end
                        if( k == 4 )
                            if( sex(order) == 0 )
                                sex(order) = 2;
                            end
                            if( sex(order) ~= 2 )
                                error = 1;
                            end                            
                        end
                        if( error ~= 0 )
                            disp('error in family structure');
                            disp('conflicts in parents sex');
                            return;
                        end
                    end
                end
            end
        end
                
        % recoded family
        % column 1, 2, 5, 8, 12
        % fields for both defined and undefined individuals
        recodeF = [];
        recodeF(1:nind,1:8) = 0;
        recodeF(1:nind,1) = f_id;
        recodeF(1:nind,2) = 1:nind;   
        recodeF(1:nind,5) = sex(1:nind);       
        recodeF(1:nind,8) = original_id(1:nind);
        recodeF(1:nind,12) = f_id;
        
        

        % fields for defined individuals only
        % update disease status column 6, genotyped status column 7
        % update parents ids column 3 and 4
        % update sequential order column 11
        % double-check column 1, 2, 5, 8, 12
        
        
        for j = 1:rows % assign all 2-column individuals
            order = assigned(family(j,2));
            if( family(j,3) ~= 0 )
                father_id = assigned(family(j,3));
            else
                father_id = 0;
            end
            if( family(j,4) ~= 0 )
                mother_id = assigned(family(j,4));
            else
                mother_id = 0;
            end
            recodeF(order,3:4) = [father_id, mother_id];
            if( recodeF(order,5) ~= family(j,5) )
                error = 1;
                disp('error in processing families');
                return;
            end
            recodeF(order,6) = family(j,6);
            recodeF(order,7) = family(j,7);
            if( recodeF(order,8) ~= family(j,2) )
                error = 1;
                disp('error in processing families');
                return;
            end         
            recodeF(order,11) = j;
            if(recodeF(order,12) ~= family(j,1))
                error = 1;
                disp('error in processing pedgree');
                return;
            end
        end

        
        % update column 9
        generation = 9;
        recodeF(1:order,generation) = 0;
        change = 1;
        while( change == 1 )
            change = 0;
            for j = 1:nind
                father_id = recodeF(j,3);
                mother_id = recodeF(j,4);
                if( father_id ~= 0 )
                    if( recodeF(j,generation) < 1 + recodeF(father_id,generation) )
                        recodeF(j,generation) = 1 + recodeF(father_id,generation);
                        change = 1;
                    end
                end
                if( mother_id ~= 0 )
                    if( recodeF(j,generation) < 1 + recodeF(mother_id,generation) )
                        recodeF(j,generation) = 1 + recodeF(mother_id,generation);
                        change = 1;
                    end
                end
            end
        end
        
        
        % reorder according to generation
        [temp,ix] = sort(recodeF(1:nind,generation));
        map(ix) = 1:length(ix);
        reorderF = recodeF(ix,1:12);
        for i = 1:length(ix)
            for j = 2:4
                if( reorderF(i,j) ~= 0 )
                    reorderF(i,j) = map(reorderF(i,j));
                end
            end
        end


end











