function [geno_alleles error] = erase_ungenotyped(seg_assignment, family)

geno_alleles = [];
error = 0;


[num_ind c] = size(family);
if( num_ind <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

seg_assignment = int32(seg_assignment);

[nrows, ncols] = size(seg_assignment);
if( nrows ~= num_ind || ncols ~= 2 )
    error = 1;
    disp('error in global IBD');
    return;
end

[error] = check_1segment(seg_assignment, family);
if( error ~= 0 )
    error = 1;
    disp('error in global IBD');
    return;
end

isgeno(1:num_ind) = (family(1:num_ind,7) == 1);



if( any(any(seg_assignment(1:num_ind,1:2) > num_ind)) || any(any(seg_assignment(1:num_ind,1:2) < -num_ind)) )
    error = 1;
    disp('error in global IBD');
    return;
end


temp(1:num_ind,1:2) = 0;
for i = 1:num_ind
    for j = 1:2
        if( seg_assignment(i,j) >= 0 )
            temp(i,j) = seg_assignment(i,j);
        else
            temp(i,j) = seg_assignment(i,j) + num_ind * 2 + 1;
        end
    end
end

map(1:2*num_ind) = 0;
for i = 1:num_ind
    if( ~isgeno(i) )
        continue;
    end
    if( map(temp(i,1)) == 0 )
        map(temp(i,1)) = -i;
    end
    if( map(temp(i,2)) == 0 )
        map(temp(i,2)) = i;
    end
end

geno_alleles = zeros(num_ind,2);
for i = 1:num_ind
    if( ~isgeno(i) )
        continue;
    end
    for j = 1:2
        geno_alleles(i,j) = map(temp(i,j));
    end
end

[error] = check_1segment(geno_alleles, family);
if( error ~= 0 )
    error = 1;
    disp('error in erasing ungenotyped alleles');
    return;
end



end














