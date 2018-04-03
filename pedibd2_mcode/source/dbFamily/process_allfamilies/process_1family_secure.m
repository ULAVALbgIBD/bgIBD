

function [output1family, state_code, error] = ...
    process_1family_secure...
    (family_range, family_genotype, parameters, option)

error = 0;
output1family = [];
state_code = [];

[error] = check_1processedfamily(family_range);
if( error ~= 0 )
    error = 1;
    disp('error in family structures');
    return;
end

full_range = family_range.pedigree_range_full;
nind = length(full_range);

if( nind <= 0 )
    error = 1;
    disp('error in family structures');
    return;
end

if( isempty(parameters) )
    error = 1;
    disp('error in marker list');
    return;
end

[nmarkers, fields] = size( parameters.sampled_markerlist );

if( nmarkers < 1 || fields ~= 2 )
    error = 1;
    disp('error in marker list');
    return;
end

if( isempty(family_genotype) )
    error = 1;
    disp('error in genotype data');
    return;
end

if( ndims(family_genotype) ~= 3 )
    error = 1;
    disp('error in genotype data');
    return;
end
[d1, d2, d3] = size(family_genotype);
if( d2 <= 0 || nind ~= d2 || d3 ~= 2 )
    error = 1;
    disp('error in genotype data');
end

if( d1 <= 0 || d1 ~= nmarkers )
    error = 1;
    disp('error in marker list');
    return;
end


genotyped_range = family_range.family_range;
if( length(genotyped_range) < 2 )
    disp(['warning: family ', num2str(family_range.family_id), ': < 2 genotyped individuals']);
end

[output1family, state_code, error] = ...
    process_1family ...
    (family_range, family_genotype, parameters, option);



end





