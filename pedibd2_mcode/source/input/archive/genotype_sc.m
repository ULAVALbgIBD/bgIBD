function result = genotype_sc(ind1)

count = 0;
result(1:length(ind1)/2) = 0;

g11 = 1;
g12 = 2;
g22 = 3;
g00 = 4;

for i = 1:2:length(ind1)
    count = count + 1;
    if( ind1(i) == 0 && ind1(i+1) == 0 )
        result(count) = g00;
    end
    if( ind1(i) == 1 && ind1(i+1) == 2 )
        result(count) = g12;
    end
    if( ind1(i) == 2 && ind1(i+1) == 1 )
        result(count) = g12;
    end
    if( ind1(i) == 1 && ind1(i+1) == 1 )
        result(count) = g11;
    end
    if( ind1(i) == 2 && ind1(i+1) == 2 )
        result(count) = g22;
    end
end
