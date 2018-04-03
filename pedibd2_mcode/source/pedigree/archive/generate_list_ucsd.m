function [ family_range pedigree_range_full child_range ] = generate_list_ucsd( pedigree_all_missing, f )

count = 0;
count_full = 0;
count_children = 0;
child_range = [];

width = length(pedigree_all_missing(1,:));

for i = 1:length(pedigree_all_missing(:,1))
    if( pedigree_all_missing(i,1) == f )
        count_full = count_full + 1;
        pedigree_range_full(count_full) = i;
        if( pedigree_all_missing(i,7) == 1 )
            count = count + 1;
            family_range(count) = count_full;
            if( pedigree_all_missing(i,3) ~= 0 && pedigree_all_missing(i,4) ~= 0 )
                count_children = count_children + 1;
                child_range(count_children) = count_full;
            end
        end

    end
end




end

