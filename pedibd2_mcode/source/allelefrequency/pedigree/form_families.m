function [family_genotype, error] = form_families(genotype_all, all_families)

error = 0;
family_genotype = [];

[d1, d2, d3] = size(genotype_all);

if( d3 ~= 2 )
    error = 1;
    disp('error in genotype data');
    return;
end

if( d1 <= 0 || d2 <= 0 )
    error = 1;
    disp('error in genotype data');
    return;
end

nind = d2;
nloc = d1;
nfam = length(all_families);
if( nfam <= 0 )
    error = 1;
    disp('no families');
    return;
end


checkmap(1:nind) = 0;
family_genotype = cell(nfam,1);

for i = 1:nfam
    pedigree_range_full = all_families{i}.pedigree_range_full;
    if( any(pedigree_range_full > nind) || any(pedigree_range_full <= 0) )
        error = 1;
        disp('error in family structures');
        return;
    end
    if( any(checkmap(pedigree_range_full)>0) )
        error = 1;
        disp('error in family structures');
        return;
    else
        checkmap(pedigree_range_full) = 1;
    end
    family_genotype{i} = genotype_all(1:nloc, pedigree_range_full, 1:2);
    family = all_families{i}.structure;
    num = length(pedigree_range_full);
    [r, c] = size( family );
    if( r <= 0 || r ~= num || c < 12 )
        error = 1;
        disp('error in family structures');
        return;
    end

end



end






