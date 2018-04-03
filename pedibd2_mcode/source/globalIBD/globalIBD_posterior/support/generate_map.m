function [map error] = generate_map(range, full)

error = 0;
map = [];

if( isempty(range) )
    error = 1;
    disp('error in family structures');
    return;
end

family = range.structure;
nind = length(range.pedigree_range_full);
[r, c] = size(range.structure);

if( nind <= 0 || r ~= nind || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

error = check_1processedfamily(range);
if( error ~= 0 )
    disp('error in famliy structures');
    return;
end

map.id = family(1:nind, 2);
map.id_genotyped = range.family_range;

if( length(unique(map.id)) ~= nind || length(map.id) ~= nind || any(map.id > nind) || any(map.id < 1))
    error = 1;
    disp('error in family structures');
    return;
end

if( isempty(map.id_genotyped) )
    error = 1;
    disp('error: no genotyped individuals in this family');
    return;
end

if( any(map.id_genotyped < 1) || any(map.id_genotyped > nind) )
    error = 1;
    disp('error in family structures');
    return;
end

for i = 1:nind
    if( map.id(i) ~= i )
        error = 1;
        disp('error in family structures');
        return;
    end
end

map.reverse_list(1:nind) = 0; %not included families members maps to 0
for i = 1:length(map.id_genotyped)
    map.reverse_list(map.id_genotyped(i)) = i;
end


map.parents = family(:, 3:4);
founder_firstchild(1:nind) = -1;
map.paternal_alleles(1:nind, 1:2) = 0;
map.maternal_alleles(1:nind, 1:2) = 0;
for i = 1:nind
    father = map.parents(i,1);
    mother = map.parents(i,2);
    if( father >= i || mother >= i )
        error = 1;
        disp('error in family structures');
        return;
    end
    if( father == 0 && mother == 0 )
        founder_firstchild(i) = 0;
    end
    if( father > 0 )
        if( founder_firstchild(father) == 0 )
            founder_firstchild(father) = i;
            if( full == 0 )
                map.paternal_alleles(i,1:2) = -father;
            else
                map.paternal_alleles(i,1:2) = [-father, father];
            end
        else
            map.paternal_alleles(i,1:2) = [-father, father];
        end
    else
        map.paternal_alleles(i,1:2) = -i;
    end
    if( mother > 0 )
        if( founder_firstchild(mother) == 0 )
            founder_firstchild(mother) = i;
            if( full == 0 )
                map.maternal_alleles(i,1:2) = -mother;
            else
                map.maternal_alleles(i,1:2) = [-mother, mother];
            end
        else
            map.maternal_alleles(i,1:2) = [-mother, mother];
        end
    else
        map.maternal_alleles(i,1:2) = i;
    end
end

end

