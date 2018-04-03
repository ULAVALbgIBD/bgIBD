function [error] = check_1processedfamily(family_range)

error = 0;
if( isempty(family_range) )
    error = 1;
    return;
end

genotyped_range = family_range.family_range;
full_range = family_range.pedigree_range_full;
family = family_range.structure;

nind = length(full_range);
if( nind <= 0 )
    error = 1;
    disp('empty family');
    return;
end

[rows, cols] = size( family );
if( rows <= 0 || rows ~= nind || cols < 12 )
    error = 1;
    disp('error in family structure');
    return;
end

% check order in the pedigree
for i = 2:nind
    if( full_range(i) ~= full_range(i-1) + 1 )
        error = 1;
        disp('error in family indexing');
        return;
    end
end

if( any(genotyped_range < 1) || any(genotyped_range > nind ) )
    error = 1;
    disp('error in family structures');
    return;
end

if( length(unique(genotyped_range)) ~= length(genotyped_range) )
    error = 1;
    disp('error in family structures');
    return;
end

for i = 1:nind
    if( family(i,2) ~= i )
        error = 1;
        disp('error in family indexing');
        return;
    end
end

for i = 1:nind
    father = family(i,3);
    mother = family(i,4);
    if( father >= i || mother >= i )
        error = 1;
        disp('error in family indexing');
        return;
    end
end


end



















