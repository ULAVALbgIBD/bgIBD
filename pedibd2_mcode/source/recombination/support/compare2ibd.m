
function [result] = compare2ibd(allele1, allele2)

result(1:2,1:2) = 0;

for i = 1:2
    for j = 1:2
        if( allele1(i) == allele2(j) )
            result(i,j) = 1;
        end
    end
end

end
