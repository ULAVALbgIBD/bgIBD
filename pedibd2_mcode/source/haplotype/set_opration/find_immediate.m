function [ ancestor, offset, assign ] = find_immediate( index )

global ancestral;

if( index < 0 )
    ancestor = ancestral(-index,1);
    offset = ancestral(-index,3);
    assign = ancestral(-index,5);
else
    ancestor = ancestral(index,2);
    offset = ancestral(index,4);
    assign = ancestral(index,6);
end

end

