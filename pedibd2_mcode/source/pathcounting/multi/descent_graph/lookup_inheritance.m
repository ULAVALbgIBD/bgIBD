

function [id, inheritance] = lookup_inheritance(vec)

    % look up a specific alle combination that is already calculated
    % not to be confused with terminal states

    global all_inheritance;
    
    id = 0;
    inheritance = [];
    list = all_inheritance.list;
    
    if(isempty(list))
        return;
    end
    
    if( length(vec) == 2 )
       id = all_inheritance.fast_index(vec(1),vec(2));
       if( id ~= 0 )
           inheritance = all_inheritance.content{id};
       end
       return;
    end    
    
    % all fields should match, vec is already in all_inheritance
    % allow permutation of alleles, setxor
    
    for i = 1:length(list(:,1))
        if( isempty( setxor(vec, list(i,:)) ) )
            id = i;
            break;
        end
    end
    
    if( id ~= 0 )
        inheritance = all_inheritance.content{id};
    end
    
end