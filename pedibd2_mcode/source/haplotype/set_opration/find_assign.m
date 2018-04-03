function [ assign ] = find_assign( index )

% assign is either 1 or 2

global ancestral;

if( index < 0 )
    assign = ancestral(-index,5);
else
    assign = ancestral(index,6);
end

end

