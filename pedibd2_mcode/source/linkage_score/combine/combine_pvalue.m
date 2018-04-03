function [combined_score error] = combine_pvalue(assignment)

error = 0;
combined_score = [];


n_fam = length(assignment);
if( n_fam <= 0 )
    error = 1;
    disp('error in family strutures');
    return;
end


for i = 1:n_fam
    if( isempty(assignment{i}) )
        error = 1;
        return;
    end
    p_value = assignment{i}.score.p_value;
    intervals = assignment{i}.intervals;
    z_score = assignment{i}.score.z_score;
    len1 = length(p_value);
    len2 = length(z_score);
    [r, c] = size(intervals);
    if( len1 ~= len2 || len2 ~= r || c < 2 )
        error = 1;
        disp('error in linkage score');
        return;
    end
    if( any( p_value > 1 | p_value < 0 ) )
        error = 1;
        disp('error in linkage score');
        return;
    end
    for j = 1:len1
        if( intervals(j,2) < intervals(j,1) )
            error = 1;
            disp('error in chromosome segmentation');
            return;
        end
    end
    for j = 2:len1
        if( intervals(j,1) ~= intervals(j-1,2) + 1 )
            error = 1;
            disp('error in chromosome segmentation');
            return;
        end
    end
    int{i} = intervals;
    p{i} = p_value;
    z{i} = z_score;
end





all_intervals = int{1};
for i = 2:n_fam
    [all_intervals, ~, ~, error] = merge2intervals(all_intervals, int{i});
    if( error ~= 0 )
        error = 1;
        disp('error in merging intervals');
        return;
    end
end

[total_len, cols] = size(all_intervals);
if( total_len <= 0 || cols ~= 2 )
    error = 1;
    disp('error in merging intervals');
    return;
end

all_p_value(1:total_len, 1:n_fam) = 0;
all_z_score(1:total_len, 1:n_fam) = 0;

for i = 1:n_fam
    [all_intervals, index1, index2, error] = merge2intervals(all_intervals, int{i});
    if( error ~= 0 )
        error = 1;
        disp('error in merging intervals');
        return;
    end
    if( length(index1) ~= total_len || length(unique(index1)) ~= total_len || length(index2) ~= total_len)
        error = 1;
        disp('error in merging intervals');
        return;
    end
    for j = 1:total_len
        if( index2(j) > length(p{i}) || index2(j) <= 0 || index2(j) > length(z{i}) )
            error = 1;
            disp('error in merging intervals');
            return;
        end
        all_p_value(j,i) = p{i}(index2(j));
        all_z_score(j,i) = z{i}(index2(j));
    end
end

combined_score.all_intervals = all_intervals;
combined_score.all_p_value = all_p_value;
combined_score.all_z_score = all_z_score;




end



























