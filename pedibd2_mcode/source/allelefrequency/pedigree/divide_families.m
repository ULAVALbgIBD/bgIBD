function [ family_genotype mendel_error error ] = divide_families(pedigree, expanded_genotype, all_families)



error = 0;
family_genotype = [];
mendel_error = [];

[d1, nind, d3] = size( expanded_genotype );
[r, c] = size( pedigree );

if( nind <= 0 || nind ~= r )
    disp('error in genotype data');
    error = 1;
    return;
end

if( c < 12 )
    error = 1;
    disp('error in pedigree structure');
    return;
end

if( d1 <= 0 || d3 ~= 2 )
    error = 1;
    disp('error in genotype data');
    return;
end

loci = d1;

% all id mapped to the correct location in the genotype file
% untyped individuals are paddled with all 0's

[family_genotype error] = form_families(expanded_genotype, all_families);
if( error ~= 0 )
    disp('error in genotype data');
    return;
end


nfam = length(family_genotype);
if( nfam <= 0 || nfam ~= length(all_families) )
    error = 1;
    disp('error in obtaining genotype for each family');
    return;
end
for i = 1:nfam
    [d1, d2, d3] = size(family_genotype{i});
    if( d2 ~= length(all_families{i}.pedigree_range_full) || d1 ~= loci || d3 ~= 2 )
        error = 1;
        disp('error in obtaining genotype for each family');
    end
end



end

