function [result result2] = generate_kinship(outer_range, f, pedigree)


for i = 1:length(outer_range)
    for j = 1:i-1
        result(i,j) = kinship(pedigree, f, outer_range(i), outer_range(j));
        result(j,i) = result(i,j);
    end
    result(i,i) = 1;
end

for i = 1:length(outer_range)
    for j = 1:length(outer_range)
        result2(i,j,1:2,1:2) = kinship2(pedigree, f, outer_range(i), outer_range(j));
    end
end

end