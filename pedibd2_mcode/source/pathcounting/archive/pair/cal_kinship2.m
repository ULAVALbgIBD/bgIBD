function [ allele_kinship ] = cal_kinship2( pedigree, family, ind1, ind2 )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

map1(2) = 1;
map2(2) = 1;
map1(3) = 1;
map2(3) = 2;
map1(4) = 2;
map2(4) = 1;
map1(5) = 2;
map2(5) = 2;


input_allpaths = linkage_compute_allele(pedigree, family, ind1, ind2);

% inbreeding part 1 and 4 of temp shall not be included actually
temp = sum(input_allpaths);

output_args = 0;

for i = 1:length(temp)
    output_args = output_args + 0.25 * temp(i) * (0.5)^i;
end

for i = [2,3,4,5]
    temp2 = input_allpaths(i,:);
    sum_temp2 = 0;
    for j = 1:length(temp2)
        sum_temp2 = sum_temp2 + 0.25 * temp2(j) * (0.5)^j;
    end
    allele_kinship(map1(i),map2(i)) = sum_temp2;
end

end

