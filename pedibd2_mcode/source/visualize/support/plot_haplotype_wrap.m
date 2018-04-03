function plot_haplotype_wrap(input, output, parameters, f, id1, id2, limit)


temp = input.family_range{f}.structure(input.family_range{f}.pedigree_range_full,8);
in1 = find(temp == id1);
in2 = find(temp == id2);
if( length(in1) ~= 1 || length(in2) ~= 1 || id1 == id2 )
    return;
end


coor = parameters.sampled_markerlist(:,1);
if( length(limit) ~= 2 )
    limit = [1,length(parameters.sampled_markerlist(:,1))];
else
    if( limit(1) > limit(2) || limit(1) < 1 || limit(2) > length(parameters.sampled_markerlist(:,1)) )
        return;
    end
end

haplotype1(:,1) = output.haplotype{f}(in1,1:2:end);
haplotype1(:,2) = output.haplotype{f}(in1,2:2:end);
haplotype2(:,1) = output.haplotype{f}(in2,1:2:end);
haplotype2(:,2) = output.haplotype{f}(in2,2:2:end);

plot_haplotype(haplotype1, haplotype2, coor, id1, id2, limit)

end