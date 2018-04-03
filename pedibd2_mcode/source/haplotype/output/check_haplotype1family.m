function [error] = check_haplotype1family(family, haplotype, markers)

error = 0;

if( isempty(family) )
    error = 1;
    disp('family structure not defined');
    return;
end

structure = family.structure;
genotyped_range = family.family_range;

[num_all, fields] = size(structure);
num_geno = length(genotyped_range);


if( num_all <= 0 )
    error = 1;
    disp('error in family structure');
    return;
end

if( fields < 7 )
    error = 1;
    disp('error in family structure');
    return;
end

[num_markers, loc] = size(markers);
if( num_markers <= 0 || loc ~= 2 )
    error = 1;
    disp('error in marker format');
    return;
end

if( ndims(haplotype) ~= 3 )
    error = 1;
    disp('error in haplotype result');
    return;
end

[rows, cols, d3] = size(haplotype);

if( d3 ~= 2 )
    error = 1;
    disp('error in haplotype result');
    return;
end

if( cols <= 0 )
    error = 1;
    disp('empty haplotype result');
    return;
end

if( cols ~= num_all )
    error = 1;
    disp('error in haplotype result');
    return;
end

if( rows <= 0 || rows ~= num_markers )
    error = 1;
    disp('error in haplotype result');
    return;
end




end


