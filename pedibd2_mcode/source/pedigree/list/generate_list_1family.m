function [ family_range pedigree_range_full child_range error ] = generate_list_1family( pedigree, f )

error = 0;

count = 0;
count_full = 0;
count_children = 0;
family_range = [];
pedigree_range_full = [];
child_range = [];


[rows, cols] = size(pedigree);
if( rows <= 0 || cols < 12 )
    error = 1;
    disp('error in pedigree structures');
    return;
end

nind = nnz(pedigree(:,1) == f);
pedigree_range_full(1:nind) = 0;

if( nind <= 0 )
    error = 1;
    disp('error in pedigree structures');
    return;
end

for i = 1:length(pedigree(:,1))
    if( pedigree(i,1) == f )
        count_full = count_full + 1;
        pedigree_range_full(count_full) = i;
        if( pedigree(i,7) == 1 )
            count = count + 1;
            family_range(count) = count_full;
            if( pedigree(i,3) ~= 0 && pedigree(i,4) ~= 0 )
                count_children = count_children + 1;
                child_range(count_children) = count_full;
            end
        end
    end
end

if( any(pedigree_range_full <= 0) )
    error = 1;
    disp('error in processing family structures');
    return;
end


end

