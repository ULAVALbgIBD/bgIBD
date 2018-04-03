function [geno, map, error] = ...
    generate_markermap(boundaries, ...
    family, family_genotype, ...
    markers1, markers2, ...
    listIND2)

error = 0;
geno = [];
map = [];

if( isempty(family) )
    error = 1;
    disp('error in family structure');
    return;
end

full_range = family.pedigree_range_full;
nIND = length(full_range);

nLOC1 = size(markers1, 1);
nLOC2 = size(markers2, 1);
[d1, d2, d3] = size(family_genotype);

nIND2 = size(listIND2, 1);

if( d2 ~= nIND2 )
    error = 1;
    disp('genotype do not match family structure');
    return;
end

if( d3 ~= 2 || d1 ~= nLOC2 )
    error = 1;
    disp('error in genotype data');
    return;
end

[nSEG, d2] = size(boundaries);
if( nSEG <= 0 || d2 ~= 2 )
    error = 1;
    disp('error in segment boundaries');
end

for i = 1:nSEG
    if( boundaries(i,1) > boundaries(i,2) )
        error = 1;
        disp('error in segment boundaries');
        return;
    end
    if( i >= 2 )
        if( boundaries(i,1) ~= boundaries(i-1,2) + 1 )
            error = 1;
            disp('error in segment boundaries');
            return;
        end
    end
end

start_marker = boundaries(1,1);
end_marker = boundaries(nSEG,2);

if( nLOC2 <= 0 )
    error = 1;
    disp('no genotype data');
    return;
end

if( start_marker ~= 1 || end_marker ~= nLOC1 )
    error = 1;
    disp('warning: segments do not cover all markers');
    return;
end


map = zeros(nLOC2,1);

for i = 1:nSEG
    a = boundaries(i,1);
    b = boundaries(i,2);
    if( a > b )
        error = 1;
        disp('error in segment boundaries');
        return;
    end
    bit = markers2(:,2) <= markers1(b,2) ...
        & markers2(:,2) >= markers1(a,2);
    if( any(map(bit) ~= 0) )
        error = 1;
        disp('segmentation error, overlapping');
        return;
    end
    map(bit) = i;
end

if( any(map(1:nLOC2) == 0) )
    disp([num2str(nnz(map > 0)), ' loci covered out of ', num2str(nLOC2), ' loci']);
end
    
[bit12, list12] = ismember(family.structure(:,8), listIND2); 
geno = int8(zeros(nLOC2, nIND, 2));
geno(1:nLOC2, bit12, 1:2) = family_genotype(1:nLOC2, list12(bit12), 1:2);




end














