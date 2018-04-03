function [genotype2pair pair_missingcode] = genotype_code(allele2genotype)

g11 = allele2genotype(1,1);
g12 = allele2genotype(1,2);
g21 = allele2genotype(2,1);
g22 = allele2genotype(2,2);

range = unique([g11,g12,g21,g22]);
len = length(range);

if( max(range) - min(range) ~= len - 1 )
    disp('allele coding not compact');
end

if( isempty(range) )
    disp('error in allele coding');
end

genotype2pair(range,range) = 0;

for i = 1:len
    for j = 1:len
        genotype2pair(range(i),range(j)) = len*(i-1)+j;
    end
end

% generate consecutive coding for genotype pair

pair_missingcode = len*len + 1;

end