function result = get_freq(seq)

count = 0;
result(1:3) = 0;
for i = 1:length(seq)
    count = count + 1;
    if( seq(i) == 0 )
        result(3) = result(3) + 1;
    end
    if( seq(i) == 1 )
        result(1) = result(1) + 1;
    end
    if( seq(i) == 2 )
        result(2) = result(2) + 1;
    end
end
