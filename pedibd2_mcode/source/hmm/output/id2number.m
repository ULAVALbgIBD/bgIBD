

function [map] = id2num(is)



for i = 1:length(is(:,1))
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    if( is(i,1) == is(i,3) )
        temp1 = temp1 + 1;
    end
    if( is(i,2) == is(i,4) )
        temp1 = temp1 + 1;
    end
    if( is(i,1) == is(i,4) )
        temp2 = temp2 + 1;
    end
    if( is(i,2) == is(i,3) )
        temp2 = temp2 + 1;
    end
    
    if( is(i,1) == is(i,2) )
        temp3 = temp3 + 1;
    end
    
    if( is(i,3) == is(i,4) )
        temp4 = temp4 + 1;
    end
    
    map(i,1) = max(temp1,temp2);
    map(i,2) = temp3;
    map(i,3) = temp4;
end




%%







