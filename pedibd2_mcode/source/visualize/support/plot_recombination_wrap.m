function plot_recombination_wrap(input, output, parameters, f, id1, id2, limit)


temp = input.family_range{f}.structure(input.family_range{f}.pedigree_range_full,8);
in1 = find(temp == id1);
in2 = find(temp == id2);
if( length(in1) ~= 1 || length(in2) ~= 1 || id1 == id2 )
    return;
end


coor = parameters.sampled_markerlist(:,2);
if( length(limit) ~= 2 )
    limit = [1,length(parameters.sampled_markerlist(:,1))];
else
    if( limit(1) > limit(2) || limit(1) < 1 || limit(2) > length(parameters.sampled_markerlist(:,1)) )
        return;
    end
end

intervals = output.assignment{f}.intervals;
[nseg c] = size(intervals);
if( nseg <= 0 || c <= 3 )
    return;
end
alleles_all = output.assignment{f}.alleles_all;
nloc = length(coor);
haplotype1 = zeros(nloc, 2);
haplotype2 = zeros(nloc, 2);

for i = 1:nseg
    haplotype1(intervals(i,1):intervals(i,2),1) = alleles_all{i}(in1,1);
    haplotype1(intervals(i,1):intervals(i,2),2) = alleles_all{i}(in1,2);
    haplotype2(intervals(i,1):intervals(i,2),1) = alleles_all{i}(in2,1);
    haplotype2(intervals(i,1):intervals(i,2),2) = alleles_all{i}(in2,2);    
end

plot_ibd(haplotype1, haplotype2, coor, id1, id2, limit);

end