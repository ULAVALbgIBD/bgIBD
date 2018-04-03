function [i1, i2] = count_ibd1(a,b)
    i1 = 0;
    i2 = 0;
    for i = 1:2
        for j = 1:2
            if( a(i) == b(j) )
                i1 = i;
                i2 = j;
            end
        end
    end
end