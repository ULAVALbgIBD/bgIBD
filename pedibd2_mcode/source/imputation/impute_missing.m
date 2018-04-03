function [imputed_alleles error] = impute_missing(assignment, family_range, kinship2ex)

global debug_mode;
imputed_alleles = [];
error = 0;

[alleles_all, intervals, error] = check_allsegments(assignment, family_range);
if( error ~= 0 )
    disp('error in processing global IBD');
    return;
end


genotyped = family_range.family_range;
nIND = length(family_range.pedigree_range_full);
if( nIND <= 0 || any(genotyped > nIND) || any(genotyped <= 0) )
    error = 1;
    disp('error in family structures');
    return;
end

family = family_range.structure;
[r c] = size(family);
if( r <= 0 || r ~= nIND || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

[r c] = size(intervals);
if( r <= 0 || c < 2 )
    error = 1;
    disp('error in segmentation');
    return;
end

if( ndims(alleles_all) ~= 3 )
    error = 1;
    disp('error in global IBD');
    return;
else
    [nSEG, d2, d3] = size(alleles_all);
    if( nSEG <= 0 || d2 ~= nIND || d3 ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
end


if( debug_mode == 1 )
    disp(' ');
    disp('imputing missing individuals');
end

imputed_alleles = alleles_all;
for i = 1:nSEG
    if( debug_mode == 1 )
        disp(['          segment ', num2str(i)]);
    end
    [geno_alleles error] = erase_ungenotyped(reshape(alleles_all(i, 1:nIND, 1:2), [nIND,2]), family);
    if( error ~= 0 )
        disp('error in global IBD');
        return;
    end
    if( i == 127 )
    end
    [imputed_alleles(i, 1:nIND, 1:2) error] ...
        = imputation(geno_alleles, family_range, kinship2ex);
    if( error ~= 0 )
        disp('error in imputing missing family members');
        return;
    end
    % one more iteration
    nITERATION = 1;
    for j = 1:nITERATION
        [imputed_alleles(i, 1:nIND, 1:2) error] ...
            = imputation(reshape(imputed_alleles(i, 1:nIND, 1:2), [nIND,2]), ...
            family_range, kinship2ex);
    end
    if( error ~= 0 )
        disp('error in imputing missing family members');
        return;
    end    
end


end

