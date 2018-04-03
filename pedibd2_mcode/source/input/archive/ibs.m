function result = ibs(ind1, ind2)

count = 0;
result(1:length(ind1)/2) = 0;

% 0: ibs0
% 1: ibs1
% 2: ibs2
% 3: missing genotype

for i = 1:2:length(ind1)
    count = count + 1;
    if( ind1(i) == 0 || ind1(i+1) == 0 )
        result(count) = 3;
        continue;
    end
    if( ind2(i) == 0 || ind2(i+1) == 0 )
        result(count) = 3;
        continue;
    end
    sum1 = 0;
    sum2 = 0;
    if( ind1(i) == ind2(i) )
        sum1 = sum1 + 1;
    end
    if( ind1(i+1) == ind2(i+1) )
        sum1 = sum1 + 1;
    end
    if( ind1(i) == ind2(i+1) )
        sum2 = sum2 + 1;
    end
    if( ind1(i+1) == ind2(i) )
        sum2 = sum2 + 1;
    end
    if( sum1 < sum2 )
        result(count) = sum2;
    else
        result(count) = sum1;
    end
end
