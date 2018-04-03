

function [id] = add_inheritance(vec, genF)

    % look up a specific alle combination that is already calculated
    % not to be confused with terminal states

    global all_inheritance;
    
    [id, inheritance] = lookup_inheritance(vec);
    if( id ~= 0 )
        return;
    end
        
    if( isempty(all_inheritance.list) )
        id = 1;
    else
        len = length(all_inheritance.list(:,1));
        id = len + 1;
    end
    
    all_inheritance.list(id,:) = vec;
    all_inheritance.content{id} = genF;
    if( length(vec) == 2 )
        all_inheritance.fast_index(vec(1),vec(2)) = id;
        all_inheritance.fast_index(vec(2),vec(1)) = id;
    end
    
end