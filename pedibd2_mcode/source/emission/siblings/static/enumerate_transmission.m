function [trans genotype] = enumerate_transmission(genotype1, genotype2)
    
count = 0;
for i1 = 1:2
    for i2 = 1:2
        for j1 = 1:2
            for j2 = 1:2
                count = count + 1;
                trans(count,1:4) = [i1,i2,j1,j2];
                genotype(count,1:4) = [genotype1(i1),genotype2(i2),genotype1(j1),genotype2(j2)];
            end
        end
    end
end
            
end