
%make sure parents are relative indexed by id

function lk_ancestor = ancestor_matrix( parents )

    lk_ancestor(1:max(parents(:,1)), 1:max(parents(:,1))) = 0;
    for j = 1:length(parents(:,1))
        for temp4 = 2:3
            parent = parents(j,temp4);
            id = parents(j,1);
            if( parent ~= 0 )
                lk_ancestor(id,parent) = 1;
                lk_ancestor(id,:) = max(lk_ancestor(id,:), lk_ancestor(parent,:));
                for k = 1:length(lk_ancestor(:,id))
                    if( lk_ancestor(k,id) == 1 )
                        lk_ancestor(k,:) = max(lk_ancestor(k,:), lk_ancestor(id,:));
                    end
                end
            end
        end
    end

    
end

