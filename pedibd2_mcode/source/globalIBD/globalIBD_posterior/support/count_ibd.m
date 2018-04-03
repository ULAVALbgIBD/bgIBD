function result = count_ibd(a,b)
    sum1 = 0;
    sum2 = 0;
    if( a(1) == b(1) )
        sum1 = sum1 + 1;
    end
    if( a(2) == b(2) )
        sum1 = sum1 + 1;
    end
    if( a(1) == b(2) )
        sum2 = sum2 + 1;
    end
    if( a(2) == b(1) )
        sum2 = sum2 + 1;
    end
    result = max( sum1, sum2 );
end