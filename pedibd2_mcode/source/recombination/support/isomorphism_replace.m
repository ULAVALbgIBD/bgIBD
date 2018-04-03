function new_config = isomorphism_replace(config)

new_config = [];

[rows, cols] = size(config);
if( rows <= 0 || cols ~= 2 )
    return;
end

all_alleles = unique(config);
alleles_dist = histc(config(:,1), all_alleles) + histc(config(:,2), all_alleles);

len1 = length(alleles_dist);
if( len1 <= 0 )
    return;
end

len2 = length(all_alleles);
if( len2 ~= len1 )
    return;
end


new_config(1:rows, 1:2) = 0;
for i = 1:rows
    for j = 1:2
        allele = config(i,j);
        if( nnz(all_alleles == allele) ~= 1 )
            return;
        end
        new_config(i,j) = alleles_dist(all_alleles == allele);
    end
end
    

end

