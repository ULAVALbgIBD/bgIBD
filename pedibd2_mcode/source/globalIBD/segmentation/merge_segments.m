function [merged_assignment error] = merge_segments(assignment)

merged_assignment = [];
merged_alleles = [];
merged_intervals = [];
error = 0;

if( isempty(assignment) )
    disp('global IBD not processed');
    return;
end

alleles_all = assignment.alleles_all;
intervals = assignment.intervals;
if( ndims(alleles_all) ~= 3 )
    disp('error in allele assignment');
    error = 1;
    return;
else
    [num_seg, nind, c] = size(alleles_all);
    if( c ~= 2 )
        error = 1;
        disp('error in allele assignment');
        return;
    end
end

if( num_seg <= 0 )
    disp('some segments not assigned');
    error = 1;
    return;
end



sbp = intervals(1,1);
tbp = intervals(1,2);
if( tbp < sbp )
    disp('error in segmentation');
    error = 1;
    return;
end
for i = 2:num_seg
    new_sbp = intervals(i,1);
    new_tbp = intervals(i,2);    
    if( new_tbp < new_sbp || new_sbp <= tbp )
        disp('error in segmentation');
        error = 1;
        return;
    end
    sbp = new_sbp;
    tbp = new_tbp;
end


% merge intervals of the same allele assignment
merged_seg = 0;
segment_range = [];

pre_allele_all = reshape(alleles_all(1, 1:nind, 1:2), [nind,2]);
pre_index = 1;
segment_range = zeros(num_seg,2);
for i = 1:num_seg
    cur_allele_all = reshape(alleles_all(i, 1:nind, 1:2), [nind,2]);
    cur_index = i;
    if( nnz(pre_allele_all - cur_allele_all) > 0 )
        ext1 = allele2pair(pre_allele_all);
        ext2 = allele2pair(cur_allele_all);
        if( nnz(ext1 - ext2) > 0 )
            merged_seg = merged_seg + 1;
            segment_range(merged_seg,1:2) = [pre_index, cur_index-1];
            pre_index = cur_index;
            pre_allele_all = cur_allele_all;
        end
    end
end

merged_seg = merged_seg + 1;
segment_range(merged_seg,1:2) = [pre_index, cur_index];


merged_intervals = zeros(merged_seg, 2);
merged_alleles_all = zeros(merged_seg, nind, 2);
for i = 1:merged_seg
    s_index = segment_range(i,1);
    t_index = segment_range(i,2);
    merged_intervals(i,1:2) = [intervals(s_index,1), intervals(t_index,2)];
    merged_alleles_all(i, 1:nind, 1:2) = assignment.alleles_all(s_index, 1:nind, 1:2);
    if( nnz( (alleles_all(s_index,1:nind,1:2)-alleles_all(t_index, 1:nind, 1:2))) > 0 )
        ext1 = allele2pair(reshape(alleles_all(s_index,1:nind,1:2), [nind,2]));
        ext2 = allele2pair(reshape(alleles_all(t_index,1:nind,1:2), [nind,2]));
        if( nnz(ext1 - ext2) > 0 )
            disp('discrepancy in merging segments');
        end
    end
end

merged_assignment.intervals = merged_intervals;
merged_assignment.alleles_all = merged_alleles_all;

disp(['total consolidated chromosomal regions: ', num2str(merged_seg)]);
        
end











